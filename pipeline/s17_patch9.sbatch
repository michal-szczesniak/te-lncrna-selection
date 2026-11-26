#!/bin/bash
#SBATCH --job-name=ridl_s17_patch9
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x_%A_%a.log
# Patch 9: candidate loci for the illustrative figure (step 44). Outputs in work/rev/p9/ and results/locus/.
# Usage: sbatch s17_patch9.sbatch locus
source "$PIPE/config_rev.sh"
P9=$RW/p9
mkdir -p $P9
case "$1" in
locus)
    python $S/44_locus_candidates.py --sites $RW/splice_age_sites.tsv.gz --usage $RW/junction_usage.tsv.gz \
        --copies $W/step3/copies_grouped.tsv.gz --map-align $W/step3/map_align.parquet \
        --consensus $W/step3/consensus_rec.fa --gtf $GTF --genome-fasta $D/hg38.fa $(species_args) \
        --flank 40 --top 12 --out $P9/locus
    O=results/locus; mkdir -p $O
    cp $P9/locus/candidates.tsv $P9/locus/bundle.json $O/
    (cd $O && md5sum * > MD5SUMS)
    tar czf results/locus_candidates.tgz -C results locus
    ls -l results/locus_candidates.tgz ;;
*) echo "unknown part: $1"; exit 1 ;;
esac
