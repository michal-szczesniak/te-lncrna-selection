#!/usr/bin/env python3
"""Revision step 24. Positive control for TE copies: when part of a TE copy IS conserved across species, do human
variants show it? This tells how large a population signal real constraint on TE sequence produces with these
data and this design, and therefore how sensitive the exon-intron test is.

For each group, every unit (copy part within one compartment) is scanned with a --window-nt sliding mean of a
phyloP track (as in steps 10 and 15); bases inside windows whose mean exceeds --threshold form 'conserved
stretches'. Population data are not used for the selection, so the comparison below is not circular, apart
from local mutation-rate differences, which comparison B removes.
  A  across copies: bases in conserved stretches vs bases at the same consensus positions (bins of --cons-bin
     nt) in copies without any conserved stretch; strata subfamily_orig x consensus bin x trinucleotide x CpG
  B  within copies: bases in conserved stretches vs the other bases of the same copies (strata copy x
     trinucleotide x CpG); this controls for the local mutation rate of the copy
Effects (as in the paper): Mantel-Haenszel log OR of a site being variable (callable sites; negative = fewer
variants) and the difference in the share of singletons; 95% intervals from --n-boot resamples of genes
(intergenic units: of units). Also reported: the share of units with a stretch and of bases inside stretches,
by compartment. Use --track prim for primate-specific groups (Alu).
--threshold-quantile q replaces the fixed threshold by the q-quantile of all window means of the group (e.g.
0.99 = the top 1% of windows); needed for primate phyloP, whose 20-nt means rarely exceed 1 because the
primate tree is short. The threshold actually used is reported in the 'threshold' column. A group without any
stretch gives one row with the stretch shares and no test.

Output: <out>_posctrl.tsv
"""
import argparse

import numpy as np
import pandas as pd

from ridl_stats import load_annot
from ridl_utils import log


def windows(df, col, W):
    """start/end indices and means of all W-nt windows of consecutive bases within one unit (>= 80% scored)."""
    v = df[col].to_numpy(dtype=float)
    ok = np.isfinite(v)
    gid = df["_u"].to_numpy()
    n = len(df)
    cs = np.concatenate([[0.0], np.cumsum(np.where(ok, v, 0.0))])
    cn = np.concatenate([[0], np.cumsum(ok)])
    i = np.arange(max(n - W + 1, 0))
    j = i + W
    cnt = cn[j] - cn[i]
    keep = (gid[i] == gid[j - 1]) & (cnt >= 0.8 * W)
    i, j, cnt = i[keep], j[keep], cnt[keep]
    return i, j, (cs[j] - cs[i]) / np.maximum(cnt, 1)


def stretch_mask(n, win, thr):
    """bases lying in at least one window with mean > thr."""
    i, j, m = win
    good = m > thr
    mark = np.zeros(n + 1, dtype=np.int64)
    np.add.at(mark, i[good], 1)
    np.add.at(mark, j[good], -1)
    return np.cumsum(mark)[:n] > 0


def mh(lab, strat, poly, sing, w=None):
    """MH log OR (lab=1 vs 0) of poly within strata, and stratified singleton-share difference."""
    w = np.ones(len(lab)) if w is None else w
    s_idx = pd.factorize(strat)[0]
    ns = s_idx.max() + 1 if len(s_idx) else 0
    bc = lambda m: np.bincount(s_idx, weights=w * m, minlength=ns)
    L = lab.astype(bool)
    A, B, C, D = bc(L & poly), bc(L & ~poly), bc(~L & poly), bc(~L & ~poly)
    n = A + B + C + D
    with np.errstate(divide="ignore", invalid="ignore"):
        R = np.nansum(np.where(n > 0, A * D / n, 0))
        S = np.nansum(np.where(n > 0, B * C / n, 0))
        lor = np.log(R / S) if R > 0 and S > 0 else np.nan
        s1, p1, s0, p0 = bc(L & sing), bc(L & poly), bc(~L & sing), bc(~L & poly)
        h = np.where((p1 > 0) & (p0 > 0), p1 * p0 / np.maximum(p1 + p0, 1e-9), 0)
        sd = (np.nansum(h * (np.where(p1 > 0, s1 / p1, 0) - np.where(p0 > 0, s0 / p0, 0))) / h.sum()
              if h.sum() > 0 else np.nan)
    return lor, sd


def run(d, lab_col, strat, n_boot, rng):
    lab = d[lab_col].to_numpy()
    poly = d["gn_poly"].to_numpy().astype(bool)
    sing = d["gn_singleton"].to_numpy().astype(bool)
    st = strat.to_numpy()
    lor, sd = mh(lab, st, poly, sing)
    clus = np.where(d["gene_id"].astype(str).to_numpy() == ".", d["unit_id"].astype(str).to_numpy(),
                    d["gene_id"].astype(str).to_numpy())
    ci, inv = np.unique(clus, return_inverse=True)
    bl, bs = [], []
    for _ in range(n_boot):
        cnt = np.bincount(rng.integers(0, len(ci), len(ci)), minlength=len(ci)).astype(float)
        x, y = mh(lab, st, poly, sing, cnt[inv])
        bl.append(x)
        bs.append(y)
    pct = lambda b: (np.nanpercentile(b, 2.5), np.nanpercentile(b, 97.5)) if np.isfinite(b).any() else (np.nan, np.nan)
    lo, hi = pct(np.array(bl, dtype=float))
    slo, shi = pct(np.array(bs, dtype=float))
    rate1 = poly[lab == 1].mean() if (lab == 1).any() else np.nan
    rate0 = poly[lab == 0].mean() if (lab == 0).any() else np.nan
    return {"bases_fg": int((lab == 1).sum()), "bases_bg": int((lab == 0).sum()), "poly_fg": rate1, "poly_bg": rate0,
            "mh_log_or": lor, "ci_lo": lo, "ci_hi": hi, "one_minus_or_pct": 100 * (1 - np.exp(lor)) if np.isfinite(lor) else np.nan,
            "singleton_diff": sd, "singleton_ci_lo": slo, "singleton_ci_hi": shi, "clusters": len(ci)}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annot", required=True)
    ap.add_argument("--subfamilies", required=True, help="comma list of groups, each analysed separately")
    ap.add_argument("--track", default="p447", help="phyloP track used to define conserved stretches")
    ap.add_argument("--window-nt", type=int, default=20)
    ap.add_argument("--threshold", type=float, default=1.5)
    ap.add_argument("--threshold-quantile", type=float, default=None,
                    help="use this quantile of the group's window means as threshold (e.g. 0.99)")
    ap.add_argument("--cons-bin", type=int, default=20)
    ap.add_argument("--query", default=None)
    ap.add_argument("--n-boot", type=int, default=200)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    col = f"phylop_{a.track}"
    cols = ["unit_id", "copy_id", "subfamily", "subfamily_orig", "compartment", "orientation", "gene_id", "chrom",
            "pos", "cons_pos", "tri", "cpg", "gn_callable", "gn_poly", "gn_singleton", col]
    rows = []
    for g in a.subfamilies.split(","):
        df = load_annot(a.annot, columns=cols, query=a.query, subfamilies=[g])
        if df.empty:
            log(f"{g}: no bases")
            continue
        df["_u"] = pd.factorize(df["unit_id"].astype(str) + "|" + df["compartment"])[0]
        df = df.sort_values(["_u", "pos"]).reset_index(drop=True)
        win = windows(df, col, a.window_nt)
        if a.threshold_quantile is not None:
            thr = float(np.quantile(win[2], a.threshold_quantile)) if len(win[2]) else np.nan
            log(f"{g}: {len(win[2])} windows; {a.threshold_quantile} quantile of window means = {thr:.3f}")
        else:
            thr = a.threshold
        df["in_stretch"] = stretch_mask(len(df), win, thr) if np.isfinite(thr) else False
        has = df.groupby("_u")["in_stretch"].transform("any")
        df["has_stretch"] = has
        df["_cbin"] = (df["cons_pos"] - 1) // a.cons_bin
        df["_ctx"] = df["tri"].astype(str) + "|" + df["cpg"].astype(str)
        for comp in ("intergenic", "intron", "exon", "all"):
            d = df if comp == "all" else df[df["compartment"] == comp]
            d = d[d["gn_callable"].astype(bool)]
            if d.empty:
                continue
            units = d.drop_duplicates("_u")
            base = {"group": g, "track": a.track, "threshold": thr, "threshold_quantile": a.threshold_quantile,
                    "compartment": comp,
                    "units": len(units), "share_units_with_stretch": units["has_stretch"].mean(),
                    "share_bases_in_stretch": d["in_stretch"].mean()}
            # A: across copies
            dA = d[d["in_stretch"] | ~d["has_stretch"]].copy()
            dA["_lab"] = dA["in_stretch"].astype(int)
            stA = dA["subfamily_orig"].astype(str) + "|" + dA["_cbin"].astype(str) + "|" + dA["_ctx"]
            if dA["_lab"].nunique() == 2:
                rows.append({**base, "comparison": "A: stretch vs same positions in copies without stretch",
                             **run(dA, "_lab", stA, a.n_boot, rng)})
            # B: within copies
            dB = d[d["has_stretch"]].copy()
            dB["_lab"] = dB["in_stretch"].astype(int)
            stB = dB["_u"].astype(str) + "|" + dB["_ctx"]
            if dB["_lab"].nunique() == 2:
                rows.append({**base, "comparison": "B: stretch vs rest of the same copies",
                             **run(dB, "_lab", stB, a.n_boot, rng)})
            if not d["in_stretch"].any():
                rows.append({**base, "comparison": "none: no conserved stretch"})
            log(f"  {g} {comp}: units {len(units)}, with stretch {base['share_units_with_stretch']:.3f}")
    out = pd.DataFrame(rows)
    out.to_csv(f"{a.out}_posctrl.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    if out.empty:
        log("no rows (no bases in any group)")
    else:
        log(out.drop(columns=["track"], errors="ignore").round(4).to_string(index=False))


if __name__ == "__main__":
    main()
