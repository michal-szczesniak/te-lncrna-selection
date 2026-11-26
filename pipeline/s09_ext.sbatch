#!/bin/bash
#SBATCH --job-name=ridl_s09_ext
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%A_%a.log
# Usage: sbatch s09_ext.sbatch prepare|annotate|summary
#        sbatch --array=1-8 s09_ext.sbatch profiles   (4 groups x main/within-gene)
source "$PIPE/config.sh"
case "$1" in
prepare)
    python $S/01_build_compartments.py $STEP1_COMMON --keep-unmatched-gtex --other-genes exons --out-dir $W/step1_ext
    python $S/02_select_te_copies.py $TE_SELECT --compartments $W/step1_ext/compartments.tsv.gz --out-dir $W/step2_ext
    map_and_group $W/step2_ext $W/step3_ext
    python $S/09_nonte_windows.py --compartments $W/step1_ext/compartments.tsv.gz $NONTE_WINDOWS --out-dir $W/nonte_ext ;;
annotate)
    annotate $W/step3_ext/map_grouped.parquet $W/step3_ext/copies_grouped.tsv.gz $W/step1_ext $W/gnomad_ext \
        $W/annot_ext $W/step3_ext/consensus_grouped.fa
    python $S/03_map_to_consensus.py --copies $W/nonte_ext/te_copies.tsv.gz --mode linear --out $W/nonte_ext/map.parquet
    annotate $W/nonte_ext/map.parquet $W/nonte_ext/te_copies.tsv.gz $W/step1_ext $W/gnomad_nonte_ext $W/annot_nonte_ext ;;
profiles)
    i=$((SLURM_ARRAY_TASK_ID - 1)); G=${MAIN_GROUPS[$((i % 4))]}; T=${G//\*/}
    if [ $i -lt 4 ]; then X=$STRATA_TE; OUT=$W/prof_ext/${T}_b20_paired; else X=$STRATA_GENE; OUT=$W/prof_ext/${T}_b20_gene; fi
    python $S/06_profiles.py --annot $W/annot_ext --subfamilies "$G" --bin-size 20 --paired-genes --query "$AUTO" $X \
        --n-perm 1000 --n-boot 200 --out-dir "$OUT" ;;
summary)
    mkdir -p $W/comp_ext $W/trend_ext
    C="--annot $W/annot_nonte_ext --subfamilies nonTE --bin-size 100 --pop-strata tri,cpg --n-perm 1000 --n-boot 200 --no-plots"
    python $S/06_profiles.py $C --out-dir $W/prof_nonte_ext/exon_vs_intron
    python $S/06_profiles.py $C --paired-genes --out-dir $W/prof_nonte_ext/exon_vs_intron_paired
    python $S/06_profiles.py $C --fg intron --bg intergenic --out-dir $W/prof_nonte_ext/intron_vs_intergenic
    python $S/14_composition.py --annot $W/annot_nonte_ext --paired-genes --query "$AUTO" --out $W/comp_ext/lnc_paired
    for G in "${MAIN_GROUPS[@]}"; do
        python $S/14_composition.py --annot $W/annot_ext --subfamilies "$G" --paired-genes --query "$AUTO" \
            --out "$W/comp_ext/te_${G//\*/}_paired"
    done
    python $S/12_gene_trend.py --annot $W/annot_nonte_ext --subfamilies nonTE --gene-table $W/gene_props.tsv \
        --query "$AUTO" --n-boot 500 --out $W/trend_ext/nonte
    python $S/12_gene_trend.py --annot $W/annot_ext --subfamilies 'AluS*,MIR*,L2abc,AluJ*' --gene-table $W/gene_props.tsv \
        --query "$AUTO" --pool-orientations --n-boot 500 --baseline $W/trend_ext/nonte_genes.tsv --out $W/trend_ext/te_vs_nonte ;;
*) echo "unknown part: $1"; exit 1 ;;
esac
stamp "done"
