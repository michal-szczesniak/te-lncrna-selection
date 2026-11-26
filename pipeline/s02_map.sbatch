#!/bin/bash
#SBATCH --job-name=ridl_s02_map
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x_%j.log
# Steps 3 and 3b for the main set: linear mapping (only for the comparison report), exact mapping from the
# RepeatMasker alignments, then subfamilies placed on a shared reference consensus.
source "$PIPE/config.sh"
mkdir -p $W/step3
python $S/03_map_to_consensus.py --copies $W/step2_v2/te_copies.tsv.gz --mode linear --out $W/step3/map_linear.parquet
map_and_group $W/step2_v2 $W/step3 "--compare-to $W/step3/map_linear.parquet"
stamp "done"
