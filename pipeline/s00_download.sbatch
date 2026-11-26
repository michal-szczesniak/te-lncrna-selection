#!/bin/bash
#SBATCH --job-name=ridl_s00_download
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%j.log
# Downloads the reference data into data/ (skip this stage if data/ links to an existing download).
# gnomAD v4.1 genome VCFs are not downloaded; stage s03 streams only the needed regions.
source "$PIPE/config.sh"
mkdir -p $D
cd $D
get () { curl -sSf -C - -O "$1"; }

stamp "GENCODE 50, GTEx v10"
get https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.annotation.gtf.gz
get https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz

stamp "UCSC hg38: chromosome sizes, genome, segmental duplications"
UCSC=https://hgdownload.soe.ucsc.edu/goldenPath/hg38
get $UCSC/bigZips/hg38.chrom.sizes
get $UCSC/bigZips/hg38.fa.gz
gunzip -c hg38.fa.gz > hg38.fa && samtools faidx hg38.fa
get $UCSC/database/genomicSuperDups.txt.gz
zcat genomicSuperDups.txt.gz | cut -f2-4 > sd.bed

stamp "RepeatMasker 4.0.6 / Dfam 2.0 annotation and alignments (same run)"
RM=https://www.repeatmasker.org/genomes/hg38/rmsk4.0.6_dfam2.0
get $RM/hg38.fa.out.gz
get $RM/hg38.fa.align.gz

stamp "phyloP 447-way and primate subset"
PP=$UCSC/phyloP447way
get $PP/md5sum.txt
get $PP/hg38.phyloP447way.bw
get $PP/hg38.phyloP447wayPrimates.bw
grep -E " hg38.phyloP447way.bw$| hg38.phyloP447wayPrimates.bw$" md5sum.txt | md5sum -c -

stamp "gnomAD v4.1 allele numbers at all genome sites -> tabix-indexed chrom/pos/AN table"
F=gnomad.genomes.v4.1.allele_number_all_sites.tsv.bgz
curl -sSf -C - -O https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.1/tsv/genomes/$F \
    || curl -sSf -C - -O https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/tsv/genomes/$F
zcat $F | tail -n +2 | sed 's/:/\t/' | bgzip -@ 4 -c > gnomad_AN_v4.1.tsv.bgz
tabix -s1 -b2 -e2 gnomad_AN_v4.1.tsv.bgz

ls -l --time-style=long-iso
stamp "done"
