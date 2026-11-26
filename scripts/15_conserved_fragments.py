#!/usr/bin/env python3
"""Conserved TE fragments: how many, which families, and where (exon / intron / intergenic)?

For every unit (one TE insertion within one compartment) of the step-4 annotation, the best LOCAL stretch
is computed: the maximum mean phyloP over --window consecutive bases (as in step 10). A unit is called
conserved when that stretch exceeds a threshold (--thresholds). For each family (--level: subfamily_orig
or subfamily) and compartment the script reports the share of conserved units, and compares exons and
intergenic regions with introns by a Mantel-Haenszel odds ratio stratified by unit length quartiles
(longer units have more windows). Intronic and intergenic copies are random samples (step-2 caps), so
shares are valid estimates; counts are not genome totals.

Processes one chromosome file at a time (memory). Outputs:
  <out>_units.tsv.gz   one row per unit: ids, family, compartment, orientation, gene, coordinates,
                       n_bases, mean and best-stretch phyloP (each track)
  <out>_summary.tsv    shares per family x compartment and length-adjusted odds ratios vs intron
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.stats import chi2

from ridl_utils import log


def unit_table(df, pcols, W):
    keys = ["subfamily", "subfamily_orig", "compartment", "orientation", "unit_id", "gene_id", "chrom"]
    for k in keys:
        df[k] = df[k].astype(str)
    df = df.sort_values(keys + ["pos"]).reset_index(drop=True)
    grp = df.groupby(keys, observed=True, sort=True)
    agg = {c: "mean" for c in pcols}
    agg["pos"] = ["min", "max", "size"]
    u = grp.agg(agg)
    u.columns = [c if isinstance(c, str) else (c[0] if c[1] == "mean" else f"pos_{c[1]}") for c in u.columns]
    u = u.reset_index().rename(columns={"pos_size": "n_bases", "pos_min": "start", "pos_max": "end"})
    u = u.rename(columns={c: f"mean_{c}" for c in pcols})
    gid = grp.ngroup().to_numpy()
    for c in pcols:
        v = df[c].to_numpy(dtype=float)
        ok = np.isfinite(v)
        cs = np.concatenate([[0.0], np.cumsum(np.where(ok, v, 0.0))])
        cn = np.concatenate([[0], np.cumsum(ok)])
        i = np.arange(max(len(df) - W + 1, 0))
        j = i + W
        cnt = cn[j] - cn[i]
        wm = np.where((gid[i] == gid[j - 1]) & (cnt >= 0.8 * W), (cs[j] - cs[i]) / np.maximum(cnt, 1), np.nan)
        u[f"best_{c}"] = pd.Series(wm).groupby(gid[i]).max().reindex(np.arange(len(u))).to_numpy()
    u["end"] = u["end"] + 1
    return u


def mh_vs_intron(g, comp, col, thr):
    """Length-stratified MH odds ratio (comp vs intron) of best stretch > thr, with CMH p."""
    s = g[g["compartment"].isin([comp, "intron"])].dropna(subset=[col]).copy()
    if s["compartment"].nunique() < 2 or len(s) < 20:
        return np.nan, np.nan
    s["lq"] = pd.qcut(s["n_bases"].rank(method="first"), 4, labels=False)
    num = den = sa = se = sv = 0.0
    for _, t in s.groupby("lq"):
        h, x = t[col] > thr, t["compartment"] == comp
        A, B, C, D = (h & x).sum(), (~h & x).sum(), (h & ~x).sum(), (~h & ~x).sum()
        n = A + B + C + D
        if n < 2:
            continue
        num += A * D / n
        den += B * C / n
        sa += A
        se += (A + B) * (A + C) / n
        sv += (A + B) * (C + D) * (A + C) * (B + D) / (n * n * (n - 1))
    orr = num / den if den > 0 else np.nan
    p = chi2.sf((abs(sa - se) - 0.5) ** 2 / sv, 1) if sv > 0 else np.nan
    return orr, p


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annot", required=True)
    ap.add_argument("--exclude-chroms", default="chrX")
    ap.add_argument("--window", type=int, default=20)
    ap.add_argument("--thresholds", default="1.0,1.5,2.0")
    ap.add_argument("--level", default="subfamily_orig", choices=["subfamily_orig", "subfamily"])
    ap.add_argument("--min-bases", type=int, default=30)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    excl = set(c for c in a.exclude_chroms.split(",") if c)
    files = sorted(glob.glob(os.path.join(a.annot, "*.parquet")))
    schema = pq.read_schema(files[0]).names
    pcols = [c for c in schema if c.startswith("phylop_")]
    cols = ["unit_id", "subfamily", "subfamily_orig", "compartment", "orientation", "gene_id", "chrom", "pos"] + pcols
    parts = []
    for f in files:
        df = pd.read_parquet(f, columns=cols)
        if df.empty or df["chrom"].iloc[0] in excl:
            continue
        parts.append(unit_table(df, pcols, a.window))
        log(f"  {os.path.basename(f)}: {len(parts[-1])} units")
    u = pd.concat(parts, ignore_index=True)
    u = u[u["n_bases"] >= a.min_bases]
    u.to_csv(f"{a.out}_units.tsv.gz", sep="\t", index=False, float_format="%.4g")

    thrs = [float(x) for x in a.thresholds.split(",")]
    main_col = "best_phylop_p447" if "best_phylop_p447" in u.columns else f"best_{pcols[0]}"
    rows = []
    for fam, g in u.groupby(a.level):
        for comp in ("exon", "intron", "intergenic"):
            s = g[g["compartment"] == comp]
            if len(s) == 0:
                continue
            r = {"family": fam, "compartment": comp, "units": len(s), "median_bases": int(s["n_bases"].median()),
                 "median_mean_p447": s["mean_phylop_p447"].median() if "mean_phylop_p447" in s else np.nan}
            for t in thrs:
                r[f"share_gt{t}"] = (s[main_col] > t).mean()
            if comp != "intron":
                r["mh_or_vs_intron_gt1.5"], r["p_vs_intron_gt1.5"] = mh_vs_intron(g, comp, main_col, 1.5)
            rows.append(r)
    s = pd.DataFrame(rows)
    s.to_csv(f"{a.out}_summary.tsv", sep="\t", index=False)
    pd.set_option("display.width", 220)
    log(s.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
