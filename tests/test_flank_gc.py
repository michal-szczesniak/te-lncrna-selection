#!/usr/bin/env python3
"""Test of 32_flank_gc.py on a planted regional-GC confound.

A synthetic chromosome consists of 3-kb blocks whose GC content varies continuously from 0.30 to 0.65. Every
block holds one 100-nt unit in its middle whose own sequence has GC 0.45 regardless of the block, so own GC
carries no information on the block. phyloP rises with block GC and has no exon effect; exonic units sit more
often in GC-rich blocks. Expected: a positive unmatched exon - intron difference that own-GC matching leaves in
place and flank-GC matching (quintiles) removes almost completely; a residual remains within quintiles.
Usage: python test_flank_gc.py <scripts dir>
"""
import os
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd
import pysam

S = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), "..", "scripts")
rng = np.random.default_rng(7)
d = tempfile.mkdtemp(prefix="flankgc_")
n_blocks, B, U = 600, 3000, 100


def seq(n, gc):
    return "".join(rng.choice(list("GCAT"), n, p=[gc / 2, gc / 2, (1 - gc) / 2, (1 - gc) / 2]))


chrom, rows = [], []
for b in range(n_blocks):
    gc = rng.uniform(0.30, 0.65)
    left = seq((B - U) // 2, gc)
    unit = seq(U, 0.45)
    right = seq(B - U - len(left), gc)
    start = b * B + len(left)
    chrom.append(left + unit + right)
    gene = f"g{b // 2}"                       # two blocks per gene
    exon = rng.random() < 0.1 + 0.8 * (gc - 0.30) / 0.35
    comp = "exon" if exon else "intron"
    for i, base in enumerate(unit):
        rows.append({"unit_id": f"u{b}", "gene_id": gene, "compartment": comp, "subfamily": "nonTE",
                     "chrom": "chrT", "pos": start + i, "ref": base, "cpg": False,
                     "phylop_p447": 2.0 * (gc - 0.30) + rng.normal(0, 0.3),
                     "gn_poly": bool(rng.random() < 0.1), "gn_callable": True})
with open(os.path.join(d, "g.fa"), "w") as fh:
    s = "".join(chrom)
    fh.write(">chrT\n" + "\n".join(s[i:i + 60] for i in range(0, len(s), 60)) + "\n")
pysam.faidx(os.path.join(d, "g.fa"))
os.makedirs(os.path.join(d, "annot"))
pd.DataFrame(rows).to_parquet(os.path.join(d, "annot", "part-chrT.parquet"), index=False)
subprocess.run([sys.executable, os.path.join(S, "32_flank_gc.py"), "--annot", os.path.join(d, "annot"),
                "--genome-fasta", os.path.join(d, "g.fa"), "--n-boot", "100", "--out", os.path.join(d, "res")],
               check=True)
r = pd.read_csv(os.path.join(d, "res_summary.tsv"), sep="\t").set_index(["metric", "comparison"])
v = {k: r.loc[("phylop_p447", k)] for k in ("unmatched", "gc_own", "gc_flank", "gc_both")}
for k, x in v.items():
    print(f"{k:10s} {x.exon_minus_intron:+.3f} [{x.ci_lo:+.3f}, {x.ci_hi:+.3f}]")
assert v["unmatched"].ci_lo > 0.08, "planted confound not visible without matching"
assert v["gc_own"].ci_lo > 0.08, "own GC should not remove a confound it does not carry"
assert abs(v["gc_flank"].exon_minus_intron) < 0.2 * v["unmatched"].exon_minus_intron, \
    "flank-GC matching should remove most of the confound"
print("PASS: flank-GC matching removes most of a planted regional-GC confound that own-GC matching keeps")
