# Settings shared by all stages. Sourced by every stage script after PIPE (the pipeline/ directory) is set.
# PROJECT must contain scripts/, figures/, pipeline/ and data/ (data/ may be a symlink to an existing download).
# Machine-specific settings (conda installation and environment, other data locations) go in pipeline/local.sh,
# which is not tracked by git; see pipeline/local.sh.example.
[ -f "$PIPE/local.sh" ] && source "$PIPE/local.sh"
PROJECT=${PROJECT:-$(cd "$PIPE/.." && pwd)}
# per-gene RNA properties of lncRNA genes from ENCODE HepG2 and K562 data, GENCODE v29 gene models
# (Supplementary Table S8; deposited at Zenodo with the data of this study)
RNA_PROPS=${RNA_PROPS:-data/gene_properties_controls_v29.tsv.gz}

set -eo pipefail
if [ -n "$CONDA_SH" ]; then
    source "$CONDA_SH"
    conda activate "$CONDA_ENV"
fi
cd "$PROJECT"
mkdir -p logs work

D=data
W=work
S=scripts
GNOMAD_URL=https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.1/vcf/genomes
GTF=$D/gencode.v50.annotation.gtf.gz
GCT=$D/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz
STEP1_COMMON="--gtf $GTF --chrom-sizes $D/hg38.chrom.sizes --gtex-gct $GCT --mask-bed $D/sd.bed"
PHYLOP="--phylop p447=$D/hg38.phyloP447way.bw --phylop prim=$D/hg38.phyloP447wayPrimates.bw"
AN="--an-tsv $D/gnomad_AN_v4.1.tsv.bgz --an-cols 1,2"
TE_SELECT="--rmsk $D/hg38.fa.out.gz --rmsk-format rmout --subfamily-regex ^(AluS|AluJ|AluY|MIR|L2) --assign perbase --max-intronic 2000 --max-intergenic 1000"
NONTE_WINDOWS="--rmsk-out $D/hg38.fa.out.gz --window 100 --max-per-comp 20000"
SF_GROUPS=(--group 'AluS*:^AluS:AluSx' --group 'AluJ*:^AluJ:AluJb' --group 'AluY*:^AluY:AluY'
        --group 'MIR*:^MIR[bc]?$:MIRb' --group 'L2abc:^L2[abc]$:L2b' --group 'L2d*:^L2d2?$:L2d')
AUTO="chrom != 'chrX'"
STRATA_TE="--pop-strata tri,cpg,subfamily_orig --cont-strata subfamily_orig"
STRATA_GENE="--pop-strata tri,cpg,subfamily_orig,gene_id --cont-strata subfamily_orig,gene_id --perm-scheme within_gene"
MAIN_GROUPS=("AluS*" "MIR*" "AluJ*" "L2abc")

gnomad_args () {  # all per-chromosome subsets of one gnomAD directory as --gnomad-vcf arguments
    for f in "$1"/gn.chr*.vcf.gz; do echo -n "--gnomad-vcf $f "; done
}

map_and_group () {  # $1 step2 dir, $2 step3 dir: exact mapping to consensus, then subfamily grouping
    mkdir -p "$2"
    python $S/03_map_to_consensus.py --copies "$1/te_copies.tsv.gz" --mode align --align-file $D/hg38.fa.align.gz \
        --consensus-out "$2/consensus_rec.fa" --out "$2/map_align.parquet" $3
    python $S/03b_group_subfamilies.py --bases "$2/map_align.parquet" --copies "$1/te_copies.tsv.gz" \
        --consensus "$2/consensus_rec.fa" "${SF_GROUPS[@]}" --out-bases "$2/map_grouped.parquet" \
        --out-copies "$2/copies_grouped.tsv.gz" --out-consensus "$2/consensus_grouped.fa" --report "$2/groups.tsv"
}

annotate () {  # $1 bases, $2 copies, $3 step1 dir, $4 gnomAD dir, $5 out dir, $6 optional consensus fasta
    local cons=""
    [ -n "$6" ] && cons="--consensus-fasta $6"
    python $S/04_annotate_bases.py --bases "$1" --copies "$2" --compartments "$3/compartments.tsv.gz" \
        --genome-fasta $D/hg38.fa --splice-sites "$3/splice_sites.tsv.gz" --ref-tx "$3/ref_tx.tsv.gz" $cons \
        $PHYLOP $(gnomad_args "$4") $AN --out-dir "$5"
}

stamp () { echo "[$(date '+%F %T')] $*"; }
