#!/usr/bin/env python3
"""Revision step 34 (patch 5). What do the junction-usage strata of step 27 mean in reads?

For every set of splice sites (lncRNA outside repeats, protein-coding, lncRNA TE-derived) the sites are joined to
the GTEx v10 junction usage of step 26 by their intron, exactly as in 27_splice_strata.py, and summarized:
  sites                 distinct splice sites of the set
  in_gtex_file          share whose intron (junction) appears in the GTEx junction table at all
  zero_reads            share in the table but with no reads in any sample
  no_reads              share in the 'no reads' stratum of step 27 (not in the table, or zero reads)
and, for the usage tertiles of step 27 (computed within each set over sites with reads > 0, ties broken by order
as in step 27): number of sites, minimum / median / maximum of total reads over all samples, median mean reads per
sample, median number of samples with >= 1 read, and median number of tissues with a mean >= 1 read.
Junctions absent from the GTEx table are not necessarily unused: the GTEx table lists junctions of its own
annotation, so junctions annotated later have no count. Both cases form the 'no reads' stratum.

Outputs: <out>_sets.tsv (one row per set), <out>_tertiles.tsv (one row per set x tertile)
"""
import argparse

import numpy as np
import pandas as pd

from ridl_utils import log

SETS = {"lnc_nonTE_site": "lncRNA outside repeats", "pc_site": "protein-coding", "lnc_TE_site": "lncRNA TE-derived"}


def tertiles(v, bins):
    """Labels 'usage q1..qk' for v > 0 (qcut on ranks, as in step 27); 'no reads' otherwise."""
    out = pd.Series("no reads", index=v.index, dtype=object)
    ok = v > 0
    if ok.sum() >= bins * 10:
        q = pd.qcut(v[ok].rank(method="first"), bins, labels=[f"usage q{i + 1}" for i in range(bins)])
        out[ok] = q.astype(str)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--positions", required=True, help="19_splice_analyze.py <out>_positions.tsv.gz")
    ap.add_argument("--usage", required=True, help="26_junction_usage.py output")
    ap.add_argument("--usage-metric", default="total_reads")
    ap.add_argument("--bins", type=int, default=3)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    pos = pd.read_csv(a.positions, sep="\t", usecols=["set", "site_id"], low_memory=False)
    pos = pos[pos["set"].isin(SETS)].drop_duplicates(["set", "site_id"])
    us = pd.read_csv(a.usage, sep="\t").set_index("intron_key")
    pos["_intron"] = pos["site_id"].astype(str).str[:-1]
    # step 26 lists every intron of the positions; 'matched' = the junction appears in the GTEx table
    matched = us["matched"].astype(bool) if "matched" in us.columns else pd.Series(True, index=us.index)
    pos["in_file"] = pos["_intron"].map(matched).fillna(False).astype(bool)
    for c in ("total_reads", "mean_reads", "samples_ge1", "tissues_mean_ge1", a.usage_metric):
        if c in us.columns:
            pos[c] = pos["_intron"].map(us[c])
    rows, trows = [], []
    for st, d in pos.groupby("set"):
        v = d[a.usage_metric].fillna(0).astype(float)
        lab = tertiles(v, a.bins)
        rows.append({"set": st, "label": SETS[st], "sites": len(d), "in_gtex_file": float(d["in_file"].mean()),
                     "zero_reads": float((d["in_file"] & (v == 0)).mean()), "no_reads": float((lab == "no reads").mean())})
        for t in sorted(set(lab) - {"no reads"}):
            x = d[lab == t]
            r = {"set": st, "label": SETS[st], "tertile": t, "sites": len(x),
                 "reads_min": float(x["total_reads"].min()), "reads_median": float(x["total_reads"].median()),
                 "reads_max": float(x["total_reads"].max())}
            for c in ("mean_reads", "samples_ge1", "tissues_mean_ge1"):
                if c in x.columns:
                    r[f"{c}_median"] = float(x[c].median())
            trows.append(r)
    S, T = pd.DataFrame(rows), pd.DataFrame(trows)
    S.to_csv(f"{a.out}_sets.tsv", sep="\t", index=False)
    T.to_csv(f"{a.out}_tertiles.tsv", sep="\t", index=False)
    pd.set_option("display.width", 220)
    log(S.round(4).to_string(index=False))
    log(T.round(2).to_string(index=False))


if __name__ == "__main__":
    main()
