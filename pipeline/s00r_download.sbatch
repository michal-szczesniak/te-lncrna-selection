#!/bin/bash
#SBATCH --job-name=ridl_s00r_download
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%j.log
# Revision downloads into data/rev/: UCSC liftOver chains hg38 -> six species and their 2bit genomes (age of
# splice sites, step 25), GTEx v10 junction read counts and sample attributes (junction usage, step 26).
source "$PIPE/config_rev.sh"
cd $R
get () { local f=$(basename "$1"); [ -s "$f" ] && { echo "have $f"; return 0; }; curl -sSf -o "$f.part" "$1" && mv "$f.part" "$f"; }
UCSC=https://hgdownload.soe.ucsc.edu/goldenPath

stamp "liftOver chains (hg38 -> species) and md5 list"
get $UCSC/hg38/liftOver/md5sum.txt
mv md5sum.txt md5sum_liftOver.txt
for sp in "${SPECIES[@]}"; do
    cap="$(tr '[:lower:]' '[:upper:]' <<< ${sp:0:1})${sp:1}"
    get $UCSC/hg38/liftOver/hg38To$cap.over.chain.gz
    grep " hg38To$cap.over.chain.gz$" md5sum_liftOver.txt | md5sum -c -
done

stamp "2bit genomes"
for sp in "${SPECIES[@]}"; do
    get $UCSC/$sp/bigZips/$sp.2bit
    if get $UCSC/$sp/bigZips/md5sum.txt; then
        grep " $sp.2bit$" md5sum.txt | md5sum -c - || { echo "md5 mismatch for $sp.2bit"; exit 1; }
        mv md5sum.txt md5sum_$sp.txt
    fi
done

stamp "GTEx v10 junction read counts (4.1 GB) and sample attributes"
get https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/GTEx_Analysis_v10_STARv2.7.10a_junctions.gct.gz
echo "51c6ba950588e21f6b4e73d85c77ed83  GTEx_Analysis_v10_STARv2.7.10a_junctions.gct.gz" | md5sum -c - \
    || echo "WARNING: md5 of the junction file differs from the CFDE record (check the download)"
get https://storage.googleapis.com/adult-gtex/annotations/v10/metadata-files/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt \
    || echo "WARNING: sample attributes not found at the expected URL; step 26 then runs without tissues"
ls -l --time-style=long-iso
date -u +%F > download_date.txt
stamp "done"
