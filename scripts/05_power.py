#!/usr/bin/env python3
"""Step 5. Power check (decision gate 1).

For every subfamily x fg-orientation x consensus bin: number of independent units and callable
sites in foreground (default exon) and background (default intron), background polymorphism
rate, and the minimal detectable RELATIVE reduction of the polymorphism rate in the foreground
(two-proportion test, Bonferroni over bins of that subfamily, given power).

Caveat: sites within one copy are not independent, so site-based MDEs are optimistic; the
unit counts are reported so you can judge this (see also the bootstrap CIs in step 6).

Outputs: power_bins.tsv, power_summary.tsv (fraction of bins with rel. MDE <= --target).
"""
import argparse
import os

import numpy as np
import pandas as pd
from scipy.stats import norm

from ridl_stats import Agg, load_annot, select_comparison
from ridl_utils import log


def rel_mde(p0, n1, n0, alpha, power):
    if not (0 < p0 < 1) or n1 < 1 or n0 < 1:
        return np.nan
    za = norm.ppf(1 - alpha / 2)

    def pw(p1):
        pbar = (n1 * p1 + n0 * p0) / (n1 + n0)
        se0 = np.sqrt(pbar * (1 - pbar) * (1 / n1 + 1 / n0))
        se1 = np.sqrt(p1 * (1 - p1) / n1 + p0 * (1 - p0) / n0)
        return norm.cdf((abs(p1 - p0) - za * se0) / se1)

    if pw(1e-9) < power:
        return np.nan  # even complete depletion is not detectable
    lo, hi = 0.0, p0  # delta range
    for _ in range(60):
        mid = (lo + hi) / 2
        if pw(p0 - mid) >= power:
            hi = mid
        else:
            lo = mid
    return hi / p0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annot", required=True, help="directory (or parquet) from step 4")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--fg", default="exon")
    ap.add_argument("--bg", default="intron")
    ap.add_argument("--bin-size", type=int, default=10)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--power", type=float, default=0.8)
    ap.add_argument("--target", type=float, default=0.10, help="relative MDE considered adequate")
    ap.add_argument("--query", default=None, help="pandas query applied to bases first")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)

    cols = ["unit_id", "subfamily", "compartment", "orientation", "gene_id", "cons_pos",
            "gn_callable", "gn_poly"]
    df = load_annot(a.annot, columns=cols, query=a.query)
    if "gn_poly" not in df.columns:
        raise SystemExit("no gnomAD columns in the annotation -- rerun step 4 with --gnomad-vcf")
    df = df[df["compartment"].isin([a.fg, a.bg]) & df["gn_callable"]][cols].copy()
    df["bin"] = (df["cons_pos"] - 1) // a.bin_size

    rows = []
    for sf in sorted(df["subfamily"].unique()):
        nb = int(df.loc[df["subfamily"] == sf, "bin"].max()) + 1
        for ori in sorted(df.loc[(df["subfamily"] == sf) & (df["compartment"] == a.fg), "orientation"].unique()):
            d = select_comparison(df, sf, ori, a.fg, a.bg)
            if d["_lab"].nunique() < 2:
                continue
            agg = Agg(d, ["gn_poly"], n_bins=nb)
            S = agg.sums(agg.label, np.ones(len(agg.units)))
            a1, a0 = S["gn_poly"]
            n1, n0 = S["_n"]
            u1, u0 = agg.units_per_bin(agg.label)
            alpha_bin = a.alpha / max(int(np.sum((n1 > 0) & (n0 > 0))), 1)
            for b in range(nb):
                p0 = a0[b] / n0[b] if n0[b] > 0 else np.nan
                rows.append({"subfamily": sf, "orientation": ori, "bin": b,
                             "cons_from": b * a.bin_size + 1, "cons_to": (b + 1) * a.bin_size,
                             "units_fg": int(u1[b]), "units_bg": int(u0[b]),
                             "sites_fg": int(n1[b]), "sites_bg": int(n0[b]),
                             "rate_fg": a1[b] / n1[b] if n1[b] > 0 else np.nan, "rate_bg": p0,
                             "rel_mde": rel_mde(p0, n1[b], n0[b], alpha_bin, a.power)})
    bins = pd.DataFrame(rows)
    bins.to_csv(os.path.join(a.out_dir, "power_bins.tsv"), sep="\t", index=False)
    summ = (bins.groupby(["subfamily", "orientation"])
                .agg(bins=("bin", "size"),
                     frac_bins_ok=("rel_mde", lambda x: float(np.mean(x <= a.target))),
                     median_rel_mde=("rel_mde", "median"),
                     median_units_fg=("units_fg", "median"), median_units_bg=("units_bg", "median"),
                     rate_bg=("rate_bg", "median"))
                .reset_index())
    summ["gate"] = np.where(summ["frac_bins_ok"] >= 0.5, "GO", "coarsen/merge")
    summ.to_csv(os.path.join(a.out_dir, "power_summary.tsv"), sep="\t", index=False)
    log("power summary:\n" + summ.to_string(index=False))


if __name__ == "__main__":
    main()
