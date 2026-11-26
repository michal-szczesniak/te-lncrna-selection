#!/usr/bin/env python3
"""Step 7. Pseudo-exon null control.

Random segments with exon-like lengths (drawn from the exonic compartment pieces of the same gene)
are placed uniformly inside the intronic compartment of each gene. Intronic BASES falling inside a
pseudo-exon are relabelled 'pseudoexon' (per base, like real exons with step 2 --assign perbase).
Rerun step 6 with
    --fg pseudoexon --bg intron
Profiles should be flat (no peaks, whole-element effects ~0). Signals here indicate artefacts
(e.g. positional or mutational biases) rather than mature-RNA constraint.
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd

from ridl_utils import index_by_chrom, log, merge_intervals


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annot", required=True, help="annotation directory from step 4")
    ap.add_argument("--compartments", required=True)
    ap.add_argument("--copies", default=None, help="(unused, kept for compatibility)")
    ap.add_argument("--out-dir", required=True, help="new annotation directory")
    ap.add_argument("--factor", type=float, default=3.0, help="pseudo-exons per real exonic piece")
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    rng = np.random.default_rng(a.seed)

    comp = pd.read_csv(a.compartments, sep="\t", dtype={"chrom": str})
    comp["L"] = comp["end"] - comp["start"]
    ex = comp[comp["compartment"] == "exon"]
    intr = comp[comp["compartment"] == "intron"]
    global_lens = ex["L"].to_numpy()
    pseudo = []
    for gid, gi in intr.groupby("gene_id"):
        lens = ex.loc[ex["gene_id"] == gid, "L"].to_numpy()
        if len(lens) == 0:
            lens = global_lens
        n = max(1, int(round(a.factor * max(1, (ex["gene_id"] == gid).sum()))))
        w = gi["L"].to_numpy().astype(float)
        for _ in range(n):
            L = int(rng.choice(lens))
            fits = w >= L
            if not fits.any():
                continue
            k = rng.choice(np.flatnonzero(fits), p=(w[fits] - L + 1) / (w[fits] - L + 1).sum())
            r = gi.iloc[k]
            s = int(r["start"]) + int(rng.integers(0, int(r["L"]) - L + 1))
            pseudo.append((r["chrom"], s, s + L, gid))
    pseudo = pd.DataFrame(pseudo, columns=["chrom", "start", "end", "gene_id"])
    log(f"  placed {len(pseudo)} pseudo-exons")

    pidx = index_by_chrom(merge_intervals(pseudo[["chrom", "start", "end"]]))
    n_bases = 0
    for f in sorted(glob.glob(os.path.join(a.annot, "*.parquet"))):
        d = pd.read_parquet(f)
        m = np.zeros(len(d), dtype=bool)
        for chrom, g in d[d["compartment"] == "intron"].groupby("chrom"):
            if chrom not in pidx:
                continue
            ps, pe = pidx[chrom]
            p = g["pos"].to_numpy()
            k = np.searchsorted(ps, p, side="right") - 1
            ins = (k >= 0) & (p < pe[np.clip(k, 0, len(pe) - 1)])
            m[d.index.get_indexer(g.index[ins])] = True
        d.loc[m, "compartment"] = "pseudoexon"
        n_bases += int(m.sum())
        d.to_parquet(os.path.join(a.out_dir, os.path.basename(f)), index=False)
    log(f"  relabelled {n_bases} intronic bases as 'pseudoexon'")
    pseudo.to_csv(os.path.join(a.out_dir, "pseudoexons.bed"), sep="\t", index=False, header=False)


if __name__ == "__main__":
    main()
