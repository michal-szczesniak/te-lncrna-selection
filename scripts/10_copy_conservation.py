#!/usr/bin/env python3
"""Copy-level conservation: is there a SUBSET of exonic TE copies that is conserved (RIDL-like tail),
even if the average exonic copy is not?

For every unit (one TE insertion within one compartment) the script computes, for each phyloP column,
(a) the mean over the unit and (b) the best LOCAL stretch: the maximum mean over --window consecutive
bases (a short conserved segment is diluted in the unit mean). Exonic and intronic units are compared:
medians (Mann-Whitney), the share of units whose mean exceeds --threshold (Fisher), and the share of units
with a local stretch above each --win-thresholds value, as a Mantel-Haenszel odds ratio stratified by unit
length quartiles (longer units have more windows, so length must be matched; CMH test).
With --paired-genes only genes having both exonic and intronic units of that group are used.
The most locally conserved exonic units are listed as candidates.

Outputs: <out>_summary.tsv, <out>_top_exonic.tsv
"""
import argparse

import numpy as np
import pandas as pd
from scipy.stats import chi2, fisher_exact, mannwhitneyu

from ridl_stats import load_annot
from ridl_utils import log


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annot", required=True)
    ap.add_argument("--subfamilies", required=True, help="comma list, e.g. 'MIR*,L2abc,AluS*'")
    ap.add_argument("--query", default=None)
    ap.add_argument("--paired-genes", action="store_true")
    ap.add_argument("--threshold", type=float, default=1.0, help="unit mean phyloP above this = conserved")
    ap.add_argument("--min-bases", type=int, default=30, help="min bases with phyloP per unit")
    ap.add_argument("--top", type=int, default=15)
    ap.add_argument("--window", type=int, default=20, help="bases in the local-stretch window")
    ap.add_argument("--win-thresholds", default="1.0,1.5,2.0")
    ap.add_argument("--out", required=True, help="output prefix")
    a = ap.parse_args()

    sfs = [s for s in a.subfamilies.split(",") if s]
    import pyarrow.parquet as pq, glob, os
    first = sorted(glob.glob(os.path.join(a.annot, "*.parquet")))[0] if os.path.isdir(a.annot) else a.annot
    pcols = [c for c in pq.read_schema(first).names if c.startswith("phylop_")]
    cols = ["unit_id", "subfamily", "compartment", "orientation", "gene_id", "chrom", "pos"] + pcols
    df = load_annot(a.annot, columns=cols, query=a.query, subfamilies=sfs)
    df = df[df["compartment"].isin(["exon", "intron"])]
    keys = ["subfamily", "compartment", "orientation", "unit_id", "gene_id", "chrom"]
    for k in keys:
        df[k] = df[k].astype(str)
    df = df.sort_values(keys + ["pos"]).reset_index(drop=True)
    agg = {c: "mean" for c in pcols}
    agg["pos"] = ["min", "max", "size"]
    grp = df.groupby(keys, observed=True, sort=True)
    u = grp.agg(agg)
    u.columns = [c if isinstance(c, str) else (c[0] if c[1] == "mean" else f"pos_{c[1]}") for c in u.columns]
    u = u.reset_index().rename(columns={"pos_size": "n_bases", "pos_min": "start", "pos_max": "end"})
    gid = grp.ngroup().to_numpy()
    W = a.window
    for c in pcols:
        v = df[c].to_numpy(dtype=float)
        ok = np.isfinite(v)
        cs = np.concatenate([[0.0], np.cumsum(np.where(ok, v, 0.0))])
        cn = np.concatenate([[0], np.cumsum(ok)])
        i = np.arange(max(len(df) - W + 1, 0))
        j = i + W
        cnt = cn[j] - cn[i]
        wm = np.where((gid[i] == gid[j - 1]) & (cnt >= 0.8 * W), (cs[j] - cs[i]) / np.maximum(cnt, 1), np.nan)
        u[f"{c}_maxwin"] = pd.Series(wm).groupby(gid[i]).max().reindex(np.arange(len(u))).to_numpy()
    u = u[u["n_bases"] >= a.min_bases]
    log(f"units: {len(u)}")

    rows, tops = [], []
    for (sf, ori), g in u.groupby(["subfamily", "orientation"]):
        if a.paired_genes:
            ge = set(g.loc[g["compartment"] == "exon", "gene_id"])
            gi = set(g.loc[g["compartment"] == "intron", "gene_id"])
            g = g[g["gene_id"].isin(ge & gi)]
        ex, it = g[g["compartment"] == "exon"], g[g["compartment"] == "intron"]
        if len(ex) < 10 or len(it) < 10:
            continue
        for c in pcols:
            e, i = ex[c].dropna(), it[c].dropna()
            he, hi = int((e > a.threshold).sum()), int((i > a.threshold).sum())
            orr, pf = fisher_exact([[he, len(e) - he], [hi, len(i) - hi]])
            pm = mannwhitneyu(e, i).pvalue if len(e) and len(i) else np.nan
            row = {"subfamily": sf, "orientation": ori, "layer": c, "units_exon": len(e), "units_intron": len(i),
                   "median_exon": e.median(), "median_intron": i.median(), "p_mannwhitney": pm,
                   f"share_exon_gt{a.threshold}": he / max(len(e), 1),
                   f"share_intron_gt{a.threshold}": hi / max(len(i), 1),
                   "odds_ratio": orr, "p_fisher": pf}
            gg = g.dropna(subset=[f"{c}_maxwin"]).copy()
            gg["lq"] = pd.qcut(gg["n_bases"].rank(method="first"), 4, labels=False)
            for thr in [float(x) for x in a.win_thresholds.split(",")]:
                hit = gg[f"{c}_maxwin"] > thr
                isx = gg["compartment"] == "exon"
                row[f"win_share_exon_gt{thr}"] = hit[isx].mean()
                row[f"win_share_intron_gt{thr}"] = hit[~isx].mean()
                num = den = sa = se = sv = 0.0
                for _, s_ in gg.assign(h=hit, x=isx).groupby("lq"):
                    A = ((s_.h) & (s_.x)).sum(); B = ((~s_.h) & (s_.x)).sum()
                    C = ((s_.h) & (~s_.x)).sum(); D = ((~s_.h) & (~s_.x)).sum()
                    n = A + B + C + D
                    if n < 2:
                        continue
                    num += A * D / n; den += B * C / n
                    sa += A; se += (A + B) * (A + C) / n
                    sv += (A + B) * (C + D) * (A + C) * (B + D) / (n * n * (n - 1))
                row[f"win_mh_or_gt{thr}"] = num / den if den > 0 else np.nan
                row[f"win_cmh_p_gt{thr}"] = chi2.sf((abs(sa - se) - 0.5) ** 2 / sv, 1) if sv > 0 else np.nan
            rows.append(row)
        main_col = "phylop_p447" if "phylop_p447" in pcols else pcols[0]
        t = ex.sort_values(f"{main_col}_maxwin", ascending=False).head(a.top).copy()
        t["start"] = t["start"].astype(int)
        t["end"] = t["end"].astype(int) + 1
        tops.append(t)
    s = pd.DataFrame(rows)
    s.to_csv(f"{a.out}_summary.tsv", sep="\t", index=False)
    top = pd.concat(tops) if tops else pd.DataFrame()
    top.to_csv(f"{a.out}_top_exonic.tsv", sep="\t", index=False)
    pd.set_option("display.width", 200)
    log(s.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
