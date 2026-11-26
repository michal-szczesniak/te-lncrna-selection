#!/bin/bash
#SBATCH --job-name=ridl_s03r_gnomad
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%A_%a.log
# Usage: sbatch --array=1-23 s03r_gnomad.sbatch ltr|l1   (task 23 = chrX)
# As s03_gnomad.sbatch, for the TE sets of the revision (LTR and L1 copies, stage s12).
source "$PIPE/config_rev.sh"
case "$1" in
    ltr) BED=$W/step2_ltr/gnomad_regions.bed; OUT=$W/gnomad_ltr ;;
    l1)  BED=$W/step2_l1/gnomad_regions.bed;  OUT=$W/gnomad_l1 ;;
    *) echo "unknown set: $1"; exit 1 ;;
esac
if [ "$SLURM_ARRAY_TASK_ID" -eq 23 ]; then CHR=chrX; else CHR=chr$SLURM_ARRAY_TASK_ID; fi
URL=$GNOMAD_URL/gnomad.genomes.v4.1.sites.$CHR.vcf.bgz
mkdir -p $OUT
awk -v c=$CHR '$1==c' $BED > $OUT/regions_$CHR.bed
stamp "$1 $CHR: $(wc -l < $OUT/regions_$CHR.bed) regions"
if [ ! -s $OUT/regions_$CHR.bed ]; then
    bcftools view -h $URL | bcftools annotate -x ^INFO/AC,INFO/AN -Oz -o $OUT/gn.$CHR.vcf.gz
else
    for attempt in 1 2 3; do
        if bcftools view -R $OUT/regions_$CHR.bed -f PASS -v snps $URL \
            | bcftools annotate -x ^INFO/AC,INFO/AN -Oz -o $OUT/gn.$CHR.vcf.gz; then
            break
        fi
        stamp "attempt $attempt failed, retrying in 60 s"
        sleep 60
        [ "$attempt" -eq 3 ] && exit 1
    done
fi
bcftools index -f -t $OUT/gn.$CHR.vcf.gz
stamp "$1 $CHR: $(bcftools view -H $OUT/gn.$CHR.vcf.gz | wc -l) variants"
