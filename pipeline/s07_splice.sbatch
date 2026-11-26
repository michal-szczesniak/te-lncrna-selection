#!/bin/bash
#SBATCH --job-name=ridl_s07_splice
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x_%j.log
# Usage: sbatch s07_splice.sbatch prepare|analyze
#   prepare: positions of TE-derived and other splice sites with matched unused positions; writes the regions
#            whose gnomAD data are not yet in work/gnomad (fetched by s03 with set "ss")
#   analyze: variant and conservation comparisons (Fig. 6)
source "$PIPE/config.sh"
case "$1" in
prepare)
    python $S/18_splice_prepare.py --gtf $GTF --lnc-genes $W/step1_all/lnc_genes.tsv.gz \
        --pc-genes $W/step1_pc/lnc_genes.tsv.gz --rmsk-out $D/hg38.fa.out.gz --genome-fasta $D/hg38.fa \
        --copies $W/step3/copies_grouped.tsv.gz --map $W/step3/map_grouped.parquet \
        --te-regions $W/step2_v2/gnomad_regions.bed --out-dir $W/splice ;;
analyze)
    python $S/19_splice_analyze.py --positions $W/splice/positions.tsv.gz \
        $(gnomad_args $W/gnomad) $(gnomad_args $W/gnomad_ss) --an-tsv $D/gnomad_AN_v4.1.tsv.bgz \
        $PHYLOP --genome-fasta $D/hg38.fa --out $W/splice/res ;;
*) echo "unknown part: $1"; exit 1 ;;
esac
stamp "done"
