#!/bin/bash
#SBATCH --job-name=ridl_s04_annotate
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=16:00:00
#SBATCH --output=logs/%x_%j.log
# Usage: sbatch s04_annotate.sbatch te|nonte|pc
source "$PIPE/config.sh"
case "$1" in
te)
    n=$(ls $W/gnomad/gn.chr*.vcf.gz.tbi | wc -l); [ "$n" -eq 23 ] || { echo "gnomAD subsets: $n of 23"; exit 1; }
    annotate $W/step3/map_grouped.parquet $W/step3/copies_grouped.tsv.gz $W/step1_all $W/gnomad $W/annot \
        $W/step3/consensus_grouped.fa
    stamp "step 5: power (Fig. 1c)"
    python $S/05_power.py --annot $W/annot --query "$AUTO" --out-dir $W/power
    for b in 20 30; do
        python $S/05_power.py --annot $W/annot --bin-size $b --query "$AUTO" --out-dir $W/power_b$b
    done ;;
nonte)
    python $S/03_map_to_consensus.py --copies $W/nonte/te_copies.tsv.gz --mode linear --out $W/nonte/map.parquet
    annotate $W/nonte/map.parquet $W/nonte/te_copies.tsv.gz $W/step1_all $W/gnomad_nonte $W/annot_nonte ;;
pc)
    python $S/03_map_to_consensus.py --copies $W/pc_windows/te_copies.tsv.gz --mode linear --out $W/pc_windows/map.parquet
    annotate $W/pc_windows/map.parquet $W/pc_windows/te_copies.tsv.gz $W/step1_pc $W/gnomad_pc $W/annot_pc ;;
*) echo "unknown set: $1"; exit 1 ;;
esac
stamp "done"
