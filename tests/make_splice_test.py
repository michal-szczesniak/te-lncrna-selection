#!/usr/bin/env python3
"""Test data for the revision splice-site scripts (25, 26, 27, 18b), built on the synthetic data set of
scripts/tests/make_synthetic_data.py.

Writes into <out>:
  genome_ss.fa(.fai)   synthetic genome with canonical GT..AG at every annotated intron (both strands)
  sp1.fa(.fai), hg38ToSp1.chain     'species 1': identical coordinates on chr1 (+ chain), chr2 reverse-complemented
                       (chain on the - strand); 2% random substitutions; every TE-derived splice site destroyed
                       with probability 0.5 (so roughly half of them are 'human-specific')
  sp2.fa(.fai), hg38ToSp2.chain     'species 2': like species 1 but every splice site destroyed with p 0.3 and
                       chr2 unaligned (absent from the chain)
  junctions.gct.gz     GTEx-like junction table (Name chrN_start1_end1, Description, 6 samples) for all
                       introns, plus some unannotated junctions sharing donors; reads ~ Poisson(lambda) with
                       lambda drawn per junction
  attributes.txt       sample -> tissue (two tissues)
Usage: make_splice_test.py <synthetic dir> <work dir with step1/> <out>
"""
import gzip
import os
import sys

import numpy as np
import pandas as pd
import pysam

S, W, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(7)
COMP = str.maketrans("ACGTN", "TGCAN")

fa = pysam.FastaFile(f"{S}/genome.fa")
seq = {c: list(fa.fetch(c)) for c in fa.references}
gtf = pd.read_csv(f"{S}/gencode.gtf.gz", sep="\t", comment="#", header=None,
                  names=["chrom", "src", "feat", "start", "end", "score", "strand", "frame", "attr"])
ex = gtf[gtf["feat"] == "exon"].copy()
ex["tx"] = ex["attr"].str.extract(r'transcript_id "([^"]+)"', expand=False)
ex["start"] -= 1
introns = []
for tx, g in ex.sort_values(["tx", "start"]).groupby("tx"):
    st, en = g["start"].to_numpy(), g["end"].to_numpy()
    for i in range(len(g) - 1):
        introns.append((g["chrom"].iloc[0], int(en[i]), int(st[i + 1]), g["strand"].iloc[0]))
introns = pd.DataFrame(introns, columns=["chrom", "start", "end", "strand"]).drop_duplicates()
for r in introns.itertuples(index=False):
    s = seq[r.chrom]
    if r.strand == "+":
        s[r.start], s[r.start + 1], s[r.end - 2], s[r.end - 1] = "G", "T", "A", "G"
    else:  # transcript GT..AG = genome CT..AC on + strand
        s[r.start], s[r.start + 1], s[r.end - 2], s[r.end - 1] = "C", "T", "A", "C"


def write_fa(path, d):
    with open(path, "w") as fh:
        for c, s in d.items():
            fh.write(f">{c}\n")
            x = "".join(s)
            for i in range(0, len(x), 80):
                fh.write(x[i:i + 80] + "\n")
    pysam.faidx(path)


write_fa(f"{OUT}/genome_ss.fa", seq)
sites = []
for r in introns.itertuples(index=False):
    sites += [(r.chrom, r.start), (r.chrom, r.start + 1), (r.chrom, r.end - 2), (r.chrom, r.end - 1)]


def species(name, p_destroy, chr2_aligned):
    d = {}
    for c, s in seq.items():
        x = np.array(s)
        m = rng.random(len(x)) < 0.02
        x[m] = rng.choice(list("ACGT"), m.sum())
        d[c] = x
    for c, p in sites:
        if rng.random() < p_destroy / 2:  # two bases per dinucleotide -> about p_destroy per site
            d[c][p] = {"A": "C", "C": "A", "G": "T", "T": "G"}.get(d[c][p], "A")
    out = {"chr1": "".join(d["chr1"])}
    if chr2_aligned:
        out["chr2r"] = "".join(d["chr2"]).translate(COMP)[::-1]
    write_fa(f"{OUT}/{name}.fa", out)
    L1, L2 = len(seq["chr1"]), len(seq["chr2"])
    with open(f"{OUT}/hg38To{name.capitalize()}.chain", "w") as fh:
        fh.write(f"chain 1000 chr1 {L1} + 0 {L1} chr1 {L1} + 0 {L1} 1\n{L1}\n\n")
        if chr2_aligned:
            # target chr2 (+) aligned to the reverse strand of query chr2r; q coordinates on the - strand
            fh.write(f"chain 1000 chr2 {L2} + 0 {L2} chr2r {L2} - 0 {L2} 2\n{L2}\n\n")


species("sp1", 0.5, True)
species("sp2", 0.3, False)

# GTEx-like junctions
samples = [f"GTEX-S{i}" for i in range(6)]
with open(f"{OUT}/attributes.txt", "w") as fh:
    fh.write("SAMPID\tSMTSD\n")
    for i, s in enumerate(samples):
        fh.write(f"{s}\t{'TissueA' if i < 3 else 'TissueB'}\n")
rows = []
for r in introns.itertuples(index=False):
    lam = rng.lognormal(1.0, 1.5)
    rows.append((f"{r.chrom}_{r.start + 1}_{r.end}", "ENSGX", rng.poisson(lam, len(samples))))
    if rng.random() < 0.3:  # alternative acceptor sharing the donor
        rows.append((f"{r.chrom}_{r.start + 1}_{r.end + 40}", "ENSGX", rng.poisson(lam / 3, len(samples))))
with gzip.open(f"{OUT}/junctions.gct.gz", "wt") as fh:
    fh.write("#1.2\n")
    fh.write(f"{len(rows)}\t{len(samples)}\n")
    fh.write("Name\tDescription\t" + "\t".join(samples) + "\n")
    for n, d, v in rows:
        fh.write(f"{n}\t{d}\t" + "\t".join(map(str, v)) + "\n")
print(f"introns {len(introns)}, junction rows {len(rows)}")
