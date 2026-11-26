#!/bin/bash
#SBATCH --job-name=ridl_s05_profiles
#SBATCH --array=1-73
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=16:00:00
#SBATCH --output=logs/%x_%A_%a.log
# Tasks  1-12  main profiles, 6 groups x (all genes, genes with exonic and intronic copies)   -> work/prof_v5
#       13-16  within-gene comparison                                                      -> work/prof_gene
#       17-20  diagnostics: no CpG sites, intron vs intergenic                             -> work/diag6
#       21-32  tertiles of GTEx expression                                                 -> work/dose_tpm
#       33-68  tertiles of stability, nuclear retention, splicing efficiency               -> work/dose_rbp
#       69-73  established lincRNAs (GENCODE v29 class), 4 groups + ordinary sequence      -> work/lincrna_v29
source "$PIPE/config.sh"
t=$SLURM_ARRAY_TASK_ID
P="python $S/06_profiles.py"

if [ $t -le 12 ]; then
    G6=("AluS*" "MIR*" "AluJ*" "L2abc" "AluY*" "MIR3"); i=$((t - 1)); G=${G6[$((i % 6))]}
    if [ $i -lt 6 ]; then X=""; TAG=""; else X="--paired-genes"; TAG="_paired"; fi
    $P --annot $W/annot --subfamilies "$G" --bin-size 20 $X --query "$AUTO" $STRATA_TE \
        --n-perm 1000 --n-boot 200 --out-dir "$W/prof_v5/${G//\*/}_b20$TAG"
elif [ $t -le 16 ]; then
    G=${MAIN_GROUPS[$((t - 13))]}
    $P --annot $W/annot --subfamilies "$G" --bin-size 20 --paired-genes --query "$AUTO" $STRATA_GENE \
        --n-perm 1000 --n-boot 200 --out-dir "$W/prof_gene/${G//\*/}_b20"
elif [ $t -le 20 ]; then
    C="--annot $W/annot --bin-size 20 $STRATA_TE --n-perm 300 --n-boot 100 --no-plots"
    case $t in
    17) $P $C --subfamilies 'AluS*' --query "$AUTO and cpg == False" --out-dir $W/diag6/AluS_noCpG ;;
    18) $P $C --subfamilies 'MIR*' --query "$AUTO and cpg == False" --out-dir $W/diag6/MIR_noCpG ;;
    19) $P $C --subfamilies 'AluS*' --query "$AUTO" --fg intron --bg intergenic --out-dir $W/diag6/AluS_intron_vs_intergenic ;;
    20) $P $C --subfamilies 'MIR*' --query "$AUTO" --fg intron --bg intergenic --out-dir $W/diag6/MIR_intron_vs_intergenic ;;
    esac
elif [ $t -le 32 ]; then
    i=$((t - 21)); G=${MAIN_GROUPS[$((i % 4))]}; T=$((i / 4))
    read LO HI < <(python - "$T" "$W/annot" << 'PY'
import glob, os, sys
import pandas as pd
t = int(sys.argv[1])
d = pd.concat([pd.read_parquet(f, columns=["gene_id", "max_tpm"]).drop_duplicates("gene_id")
               for f in glob.glob(os.path.join(sys.argv[2], "*.parquet"))]).drop_duplicates("gene_id").dropna()
q = d["max_tpm"].quantile([0, 1/3, 2/3, 1]).values
print(q[t] - 1e-9, q[t + 1] + 1e-9)
PY
)
    stamp "group $G, expression tertile $((T + 1)): max_tpm in [$LO, $HI)"
    $P --annot $W/annot --subfamilies "$G" --bin-size 50 --paired-genes \
        --query "$AUTO and max_tpm >= $LO and max_tpm < $HI" $STRATA_TE \
        --n-perm 300 --n-boot 100 --no-plots --out-dir "$W/dose_tpm/${G//\*/}_t$((T + 1))"
elif [ $t -le 68 ]; then
    i=$((t - 33)); G4=("AluS*" "MIR*" "L2abc" "AluJ*"); PROPS=(stab nuc splice)
    G=${G4[$((i % 4))]}; PR=${PROPS[$(((i / 4) % 3))]}; T=$((i / 12 + 1))
    $P --annot $W/annot --subfamilies "$G" --bin-size 50 --paired-genes --gene-table $W/gene_props.tsv \
        --query "$AUTO and ${PR}_t == $T" $STRATA_TE \
        --n-perm 300 --n-boot 100 --no-plots --out-dir "$W/dose_rbp/${G//\*/}_${PR}_t$T"
else
    SETS=("AluS*" "MIR*" "L2abc" "AluJ*" "nonTE"); G=${SETS[$((t - 69))]}
    if [ "$G" = "nonTE" ]; then A=$W/annot_nonte; ST="--pop-strata tri,cpg"; else A=$W/annot; ST=$STRATA_TE; fi
    $P --annot $A --subfamilies "$G" --bin-size 50 --paired-genes --gene-table $W/gene_class.tsv \
        --query "$AUTO and project_class == 'lincRNA'" $ST \
        --n-perm 300 --n-boot 100 --no-plots --out-dir "$W/lincrna_v29/${G//\*/}"
fi
stamp "task $t done"
