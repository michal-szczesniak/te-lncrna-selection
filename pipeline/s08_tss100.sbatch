#!/bin/bash
#SBATCH --job-name=ridl_s08_tss100
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%j.log
# Usage: sbatch s08_tss100.sbatch prepare|analyze
source "$PIPE/config.sh"
case "$1" in
prepare)
    python $S/01_build_compartments.py $STEP1_COMMON --keep-unmatched-gtex --exon-tss-pad 100 --out-dir $W/step1_tss100
    python $S/02_select_te_copies.py $TE_SELECT --compartments $W/step1_tss100/compartments.tsv.gz --out-dir $W/step2_tss100
    map_and_group $W/step2_tss100 $W/step3_tss100 ;;
analyze)
    annotate $W/step3_tss100/map_grouped.parquet $W/step3_tss100/copies_grouped.tsv.gz $W/step1_tss100 \
        $W/gnomad_tss100 $W/annot_tss100 $W/step3_tss100/consensus_grouped.fa
    mkdir -p $W/comp_tss100 $W/trend_tss100
    for G in "${MAIN_GROUPS[@]}"; do
        T=${G//\*/}
        python $S/06_profiles.py --annot $W/annot_tss100 --subfamilies "$G" --bin-size 20 --paired-genes \
            --query "$AUTO" $STRATA_TE --n-perm 300 --n-boot 200 --no-plots --out-dir "$W/prof_tss100/${T}_b20_paired"
        python $S/14_composition.py --annot $W/annot_tss100 --subfamilies "$G" --paired-genes --query "$AUTO" \
            --out "$W/comp_tss100/te_${T}_paired"
    done
    python $S/12_gene_trend.py --annot $W/annot_tss100 --subfamilies 'AluS*,MIR*,L2abc,AluJ*' \
        --gene-table $W/gene_props.tsv --query "$AUTO" --pool-orientations --n-boot 500 \
        --baseline $W/trend/nonte_genes.tsv --out $W/trend_tss100/te_vs_nonte ;;
*) echo "unknown part: $1"; exit 1 ;;
esac
stamp "done"
