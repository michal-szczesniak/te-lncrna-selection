#!/usr/bin/env bash
# End-to-end test on synthetic data with planted signals (per-base assignment + subfamily grouping).
# Expected: peaks ONLY for AluS*/sense at 101-130 (AluSx + AluSz grouped on AluSx) and MIRb/sense at
# 51-80 (all three layers); no peaks for pseudoexon-vs-intron and intron-vs-intergenic; reconstructed
# consensus = truth; per-base assignment yields more exonic copies than whole-copy assignment.
set -euo pipefail
HERE=$(cd "$(dirname "$0")/.." && pwd)
T=${1:-$(mktemp -d)}
S=$T/syn; W=$T/work
mkdir -p "$W"
cd "$HERE"
STRATA="--pop-strata tri,cpg,subfamily_orig --cont-strata subfamily_orig"

python3 tests/make_synthetic_data.py "$S"
python3 01_build_compartments.py --gtf $S/gencode.gtf.gz --chrom-sizes $S/chrom.sizes --gtex-gct $S/gtex.gct.gz \
    --mask-bed $S/sd_mask.bed --chroms chr1,chr2 --out-dir $W/step1
python3 02_select_te_copies.py --rmsk $S/rmsk.out.gz --rmsk-format rmout --subfamily-regex '^(AluS|MIR)' \
    --assign whole --compartments $W/step1/compartments.tsv.gz --chroms chr1,chr2 --out-dir $W/step2_whole
python3 02_select_te_copies.py --rmsk $S/rmsk.out.gz --rmsk-format rmout --subfamily-regex '^(AluS|MIR)' \
    --assign perbase --compartments $W/step1/compartments.tsv.gz --chroms chr1,chr2 --out-dir $W/step2
python3 03_map_to_consensus.py --copies $W/step2/te_copies.tsv.gz --mode linear --out $W/map_linear.parquet
python3 03_map_to_consensus.py --copies $W/step2/te_copies.tsv.gz --mode align --align-file $S/rmsk.align.gz \
    --consensus-out $W/consensus_rec.fa --out $W/map_align.parquet --compare-to $W/map_linear.parquet
python3 03b_group_subfamilies.py --bases $W/map_align.parquet --copies $W/step2/te_copies.tsv.gz \
    --consensus $W/consensus_rec.fa --group 'AluS*:^AluS:AluSx' --out-bases $W/map_grouped.parquet \
    --out-copies $W/copies_grouped.tsv.gz --out-consensus $W/consensus_grouped.fa --report $W/groups.tsv
python3 04_annotate_bases.py --bases $W/map_grouped.parquet --copies $W/copies_grouped.tsv.gz \
    --compartments $W/step1/compartments.tsv.gz \
    --genome-fasta $S/genome.fa --splice-sites $W/step1/splice_sites.tsv.gz --ref-tx $W/step1/ref_tx.tsv.gz \
    --consensus-fasta $W/consensus_grouped.fa --phylop p447=$S/phylop.bw \
    --gnomad-vcf $S/gnomad_chr1.vcf.gz --gnomad-vcf $S/gnomad_chr2.vcf.gz \
    --an-tsv $S/an.tsv.gz --an-cols 1,2 --out-dir $W/annot
python3 05_power.py --annot $W/annot --out-dir $W/power
python3 06_profiles.py --annot $W/annot --out-dir $W/prof_main --n-perm 300 --n-boot 200 $STRATA
python3 07_pseudoexon_control.py --annot $W/annot --compartments $W/step1/compartments.tsv.gz \
    --out-dir $W/annot_pseudo
python3 06_profiles.py --annot $W/annot_pseudo --fg pseudoexon --bg intron --min-units 10 \
    --out-dir $W/prof_pseudo --n-perm 300 --n-boot 200 --no-plots $STRATA
python3 06_profiles.py --annot $W/annot --fg intron --bg intergenic --out-dir $W/prof_ig \
    --n-perm 300 --n-boot 200 --no-plots $STRATA
python3 08_positive_controls.py --controls $S/controls.bed --annot $W/annot \
    --profiles $W/prof_main/profiles.tsv --peaks $W/prof_main/peaks.tsv --out $W/controls_report.tsv

python3 - "$W" "$S" <<'EOF'
import sys, pandas as pd, pysam
W, S = sys.argv[1], sys.argv[2]
def peaks(d):
    try:
        return pd.read_csv(f"{W}/{d}/peaks.tsv", sep="\t")
    except pd.errors.EmptyDataError:
        return pd.DataFrame()
def fasta(p):
    return {r.name.split("#")[0]: r.sequence for r in pysam.FastxFile(p)}
truth, rec = fasta(f"{S}/consensus.fa"), fasta(f"{W}/consensus_rec.fa")
for k in ("AluSx", "AluSz", "MIRb"):
    assert rec[k] == truth[k], f"reconstructed consensus differs for {k}"
grp = pd.read_csv(f"{W}/groups.tsv", sep="\t")
assert grp.loc[grp.subfamily == "AluSz", "identity"].iloc[0] > 0.9
cw = pd.read_csv(f"{W}/step2_whole/copy_counts.tsv", sep="\t")
nw = cw.loc[cw.compartment == "exon", "copies"].sum()
ann = pd.concat([pd.read_parquet(f"{W}/annot/{f}", columns=["copy_id", "compartment"])
                 for f in ("part-chr1.parquet", "part-chr2.parquet")])
nb = ann.loc[ann.compartment == "exon", "copy_id"].nunique()
assert nb > nw, f"per-base should add exonic copies ({nb} vs {nw})"
for c in ("chr1", "chr2"):
    d = pd.read_parquet(f"{W}/annot/part-{c}.parquet", columns=["gn_poly", "gn_callable"])
    r = d.loc[d.gn_callable, "gn_poly"].mean()
    assert 0.12 < r < 0.25, f"polymorphism rate on {c} is {r:.3f} (expected ~0.18): gnomAD lookup broken?"
pk = peaks("prof_main")
got = set(zip(pk["subfamily"], pk["orientation"], pk["cons_from"], pk["cons_to"]))
exp = {("AluS*", "sense", 101, 130), ("MIRb", "sense", 51, 80)}
assert got == exp, f"unexpected main peaks: {got}"
assert len(pk) == 6, "expected 3 layers x 2 planted regions"
cpk = pd.read_csv(f"{W}/prof_main/peaks_centered.tsv", sep="\t")
cgot = set(zip(cpk["subfamily"], cpk["orientation"], cpk["cons_from"], cpk["cons_to"]))
assert cgot == exp, f"unexpected centered peaks: {cgot}"
assert peaks("prof_pseudo").empty, "pseudo-exon control produced peaks"
assert peaks("prof_ig").empty, "intron-vs-intergenic produced peaks"
ctl = pd.read_csv(f"{W}/controls_report.tsv", sep="\t")
assert ctl.loc[ctl.control == "planted_sense_copy", "in_peaks"].fillna("").str.len().gt(0).all()
print(f"exonic copies: whole {nw}, per-base {nb}")
print("ALL CHECKS PASSED")
EOF
echo "outputs in $W"
