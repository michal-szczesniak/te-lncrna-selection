#!/usr/bin/env python3
"""Revision step 39 (patch 5b). Conserved stretches in exonic vs intronic units, pooled per TE group, against
repeat-free sequence, for the main compartments and the loose variants (tss0, open).

Input: step-15 unit tables (<out>_units.tsv.gz: one row per unit with subfamily, subfamily_orig, compartment,
gene_id, n_bases, best_phylop_<track>) for TE units and for repeat-free windows of the same compartments.
A unit is conserved when its best --window stretch (step 15) exceeds a threshold. For each TE group (step-3b
'subfamily', e.g. 'AluS*'), design (all genes; same genes = genes with exonic and intronic units of that group)
and threshold, the exon vs intron Mantel-Haenszel odds ratio is computed with strata = original subfamily x
unit-length quartile (within subfamily; longer units have more windows). Confidence intervals come from a gene
bootstrap (genes resampled with replacement; the same draw is applied to the TE group and to repeat-free
sequence, so the ratio of the two odds ratios - the TE-specific excess - has a joint interval).

Output: <out>.tsv  variant, group, design, track, threshold, units/conserved/share per compartment, genes, mh_or,
ci_lo, ci_hi, p_cmh, or_nonte, ratio_vs_nonte, ratio_lo, ratio_hi, n_boot_finite
"""
import argparse

import numpy as np
import pandas as pd
from scipy.stats import chi2

from ridl_utils import log


def mh_parts(st, ex, hit, w, ns):
    """Weighted 2x2 counts per stratum -> MH odds ratio."""
    A = np.bincount(st, w * (ex & hit), ns)
    B = np.bincount(st, w * (ex & ~hit), ns)
    C = np.bincount(st, w * (~ex & hit), ns)
    D = np.bincount(st, w * (~ex & ~hit), ns)
    n = A + B + C + D
    ok = n > 1
    num = np.sum(A[ok] * D[ok] / n[ok])
    den = np.sum(B[ok] * C[ok] / n[ok])
    return num / den if den > 0 else (np.inf if num > 0 else np.nan)


def cmh_p(st, ex, hit, ns):
    A = np.bincount(st, ex & hit, ns).astype(float)
    B = np.bincount(st, ex & ~hit, ns).astype(float)
    C = np.bincount(st, ~ex & hit, ns).astype(float)
    D = np.bincount(st, ~ex & ~hit, ns).astype(float)
    n = A + B + C + D
    ok = n > 1
    A, B, C, D, n = A[ok], B[ok], C[ok], D[ok], n[ok]
    e = (A + B) * (A + C) / n
    v = (A + B) * (C + D) * (A + C) * (B + D) / (n * n * (n - 1))
    if v.sum() <= 0:
        return np.nan
    return float(chi2.sf((abs(A.sum() - e.sum()) - 0.5) ** 2 / v.sum(), 1))


class Set:
    """Units of one group and design, with strata, compartment flags and gene indices."""

    def __init__(self, u, gene_index):
        u = u.copy()
        u["lq"] = 0
        for _, idx in u.groupby("subfamily_orig").groups.items():
            s = u.loc[idx, "n_bases"]
            k = min(4, len(s))
            u.loc[idx, "lq"] = pd.qcut(s.rank(method="first"), k, labels=False) if k > 1 else 0
        self.st = pd.factorize(u["subfamily_orig"].astype(str) + "|" + u["lq"].astype(str))[0]
        self.ns = int(self.st.max()) + 1 if len(u) else 0
        self.ex = (u["compartment"] == "exon").to_numpy()
        self.g = gene_index.reindex(u["_gene"]).to_numpy()
        self.u = u
        self.genes = u.loc[u["compartment"].isin(["exon", "intron"]), "_gene"].nunique()

    def hit(self, col, t):
        return (self.u[col] > t).to_numpy()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--variant", nargs=3, action="append", required=True, metavar=("NAME", "TE_UNITS", "NONTE_UNITS"),
                    help="variant name, step-15 units of TE copies, step-15 units of repeat-free windows")
    ap.add_argument("--groups", default="AluJ*,AluS*,AluY*,L2abc,MIR*,MIR3")
    ap.add_argument("--tracks", default="p447,prim")
    ap.add_argument("--thresholds", default="1.0,1.5")
    ap.add_argument("--n-boot", type=int, default=500)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    groups = a.groups.split(",")
    thrs = [float(x) for x in a.thresholds.split(",")]
    rng = np.random.default_rng(a.seed)
    rows = []
    for name, fte, fnt in a.variant:
        te = pd.read_csv(fte, sep="\t", low_memory=False)
        nt = pd.read_csv(fnt, sep="\t", low_memory=False)
        nt["subfamily"] = "nonTE"
        nt["subfamily_orig"] = "nonTE"
        U = pd.concat([te, nt], ignore_index=True)
        U = U[U["compartment"].isin(["exon", "intron"])].copy()
        gid = U["gene_id"].astype(str)
        U["_gene"] = np.where(gid.isin(["", ".", "nan", "None"]), "unit:" + U["unit_id"].astype(str), gid)
        genes = pd.Index(U["_gene"].unique())
        gene_index = pd.Series(np.arange(len(genes)), index=genes)
        tracks = [t for t in a.tracks.split(",") if f"best_phylop_{t}" in U.columns]
        log(f"{name}: {len(te)} TE units, {len(nt)} repeat-free units, {len(genes)} genes, tracks {tracks}")
        sets = {}
        for grp in groups + ["nonTE"]:
            g = U[U["subfamily"] == grp]
            if g.empty:
                log(f"  {grp}: no units")
                continue
            both = g.groupby("_gene")["compartment"].nunique()
            sets[(grp, "all genes")] = Set(g, gene_index)
            sets[(grp, "same genes")] = Set(g[g["_gene"].isin(both.index[both == 2])], gene_index)
        mult = [np.bincount(rng.integers(0, len(genes), len(genes)), minlength=len(genes)).astype(float)
                for _ in range(a.n_boot)]
        for tr in tracks:
            col = f"best_phylop_{tr}"
            for t in thrs:
                boot = {}
                for key, s in sets.items():
                    h = s.hit(col, t) & s.u[col].notna().to_numpy()
                    ok = s.u[col].notna().to_numpy()
                    st, ex, hh, gi = s.st[ok], s.ex[ok], h[ok], s.g[ok]
                    obs = mh_parts(st, ex, hh, np.ones(len(st)), s.ns) if len(st) else np.nan
                    bs = np.array([mh_parts(st, ex, hh, m[gi], s.ns) for m in mult]) if len(st) else np.full(a.n_boot, np.nan)
                    boot[key] = (obs, bs)
                    lb = np.log(bs[np.isfinite(bs) & (bs > 0)])
                    r = {"variant": name, "group": key[0], "design": key[1], "track": tr, "threshold": t,
                         "units_exon": int(ex.sum()), "units_intron": int((~ex).sum()),
                         "conserved_exon": int((ex & hh).sum()), "conserved_intron": int((~ex & hh).sum()),
                         "genes": s.genes, "mh_or": obs,
                         "ci_lo": float(np.exp(np.percentile(lb, 2.5))) if len(lb) > 20 else np.nan,
                         "ci_hi": float(np.exp(np.percentile(lb, 97.5))) if len(lb) > 20 else np.nan,
                         "p_cmh": cmh_p(st, ex, hh, s.ns) if len(st) else np.nan, "n_boot_finite": len(lb)}
                    r["share_exon"] = r["conserved_exon"] / r["units_exon"] if r["units_exon"] else np.nan
                    r["share_intron"] = r["conserved_intron"] / r["units_intron"] if r["units_intron"] else np.nan
                    rows.append(r)
                for r in rows:
                    if r["variant"] != name or r["track"] != tr or r["threshold"] != t or r["group"] == "nonTE":
                        continue
                    if "ratio_vs_nonte" in r:
                        continue
                    base = boot.get(("nonTE", r["design"]))
                    obs, bs = boot[(r["group"], r["design"])]
                    if base is None:
                        continue
                    with np.errstate(all="ignore"):
                        rat = bs / base[1]
                    lr = np.log(rat[np.isfinite(rat) & (rat > 0)])
                    r["or_nonte"] = base[0]
                    r["ratio_vs_nonte"] = obs / base[0] if base[0] and np.isfinite(base[0]) else np.nan
                    r["ratio_lo"] = float(np.exp(np.percentile(lr, 2.5))) if len(lr) > 20 else np.nan
                    r["ratio_hi"] = float(np.exp(np.percentile(lr, 97.5))) if len(lr) > 20 else np.nan
                log(f"  {name} {tr} >{t}: done")
    R = pd.DataFrame(rows)
    cols = ["variant", "group", "design", "track", "threshold", "units_exon", "units_intron", "conserved_exon",
            "conserved_intron", "share_exon", "share_intron", "genes", "mh_or", "ci_lo", "ci_hi", "p_cmh",
            "or_nonte", "ratio_vs_nonte", "ratio_lo", "ratio_hi", "n_boot_finite"]
    R = R.reindex(columns=cols)
    R.to_csv(f"{a.out}.tsv", sep="\t", index=False, float_format="%.5g")
    pd.set_option("display.width", 250)
    log(R[(R["track"] == "p447") & (R["threshold"] == 1.5)].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
