#!/usr/bin/env python3
"""Optional helper. Per-base alignment depth (number of non-reference species with an aligned,
non-gap base) for the regions from step 2, from a MAF(.gz) file or a UCSC bigMaf (.bb).

  python helper_maf_depth.py --maf hg38.cactus447way.bb --regions gnomad_regions.bed \
      --ref-prefix hg38 --subset-species primates.txt --out depth.tsv.gz

--subset-species: file with one species/assembly prefix per line (as in the MAF src field before
the first '.'); a second column depth_subset counts only those (e.g. primates).
bigMaf access uses pyBigWig.entries(); remote URLs need pyBigWig built with curl support.
Streaming a whole-genome MAF is very slow -- prefer bigMaf with region queries.
Output columns: chrom, pos (0-based), depth[, depth_subset]  -> pass to step 4 as --depth-tsv.
"""
import argparse
from collections import defaultdict

import numpy as np
import pandas as pd

from ridl_utils import index_by_chrom, log, merge_intervals, open_text


def parse_block(lines, ref_prefix, subset, want, acc):
    """lines: MAF block lines ('s ...'). want(chrom, pos)->bool. acc: dict (chrom,pos)->[d, ds]."""
    rows = [l.split() for l in lines if l.startswith("s ")]
    ref = [r for r in rows if r[1].split(".")[0] == ref_prefix]
    if not ref:
        return
    r = ref[0]
    chrom = r[1].split(".", 1)[1]
    start, strand, text = int(r[2]), r[4], r[6]
    if strand != "+":
        return
    others = [(x[1].split(".")[0], x[6]) for x in rows if x is not r]
    ref_arr = np.frombuffer(text.encode(), dtype="S1")
    nongap = ref_arr != b"-"
    pos = start + np.cumsum(nongap) - 1
    d = np.zeros(len(ref_arr), dtype=np.int32)
    ds = np.zeros(len(ref_arr), dtype=np.int32)
    for sp, t in others:
        arr = np.frombuffer(t.encode(), dtype="S1")
        ok = (arr != b"-") & (arr != b"N") & (arr != b"n")
        d += ok
        if subset and sp in subset:
            ds += ok
    for j in np.flatnonzero(nongap):
        p = int(pos[j])
        if want(chrom, p):
            acc[(chrom, p)] = (int(d[j]), int(ds[j]))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--maf", required=True)
    ap.add_argument("--regions", required=True)
    ap.add_argument("--ref-prefix", default="hg38")
    ap.add_argument("--subset-species", default=None)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    reg = pd.read_csv(a.regions, sep="\t", header=None, usecols=[0, 1, 2], names=["chrom", "start", "end"],
                      dtype={"chrom": str})
    idx = index_by_chrom(merge_intervals(reg))
    subset = set(l.strip() for l in open(a.subset_species) if l.strip()) if a.subset_species else set()

    def want(chrom, p):
        if chrom not in idx:
            return False
        s, e = idx[chrom]
        i = np.searchsorted(s, p, side="right") - 1
        return i >= 0 and p < e[i]

    acc = {}
    if a.maf.endswith(".bb"):
        import pyBigWig
        bb = pyBigWig.open(a.maf)
        for r in reg.itertuples(index=False):
            for _, _, txt in bb.entries(r.chrom, int(r.start), int(r.end)) or []:
                parse_block(txt.split(";"), a.ref_prefix, subset, want, acc)
    else:
        block = []
        with open_text(a.maf) as fh:
            for line in fh:
                if line.startswith("a"):
                    if block:
                        parse_block(block, a.ref_prefix, subset, want, acc)
                    block = []
                elif line.startswith("s "):
                    block.append(line)
            if block:
                parse_block(block, a.ref_prefix, subset, want, acc)
    out = pd.DataFrame([(c, p, d, ds) for (c, p), (d, ds) in acc.items()],
                       columns=["chrom", "pos", "depth", "depth_subset"])
    if not subset:
        out = out.drop(columns="depth_subset")
    out.sort_values(["chrom", "pos"]).to_csv(a.out, sep="\t", index=False)
    log(f"wrote depth for {len(out)} positions")


if __name__ == "__main__":
    main()
