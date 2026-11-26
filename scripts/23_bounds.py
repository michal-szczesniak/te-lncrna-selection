#!/usr/bin/env python3
"""Revision step 23. Bounds on TE-specific constraint from the gene-level test (Fig. 3b), with more resamples,
weighting, the conservation layers, an orientation contrast and the translation into functional copies.

Input: one or more per-gene tables from 12_gene_trend.py --baseline (te_vs_nonte_genes.tsv): for every gene and
TE group, exon - intron differences of the TE copies (diff_*) and their excess over repeat-free sequence of the
same gene (excess_* = diff_TE - diff_repeat-free). Layers: poly_noncpg (share of variable sites at callable
non-CpG positions), phylop_p447, phylop_prim (and phylop_consmatch if present).

For every table x group x orientation x layer:
  mean excess over genes (as in the paper) and a weighted mean (weights n_exon*n_intron/(n_exon+n_intron)),
  percentile 95% intervals from --n-boot resamples of genes, bootstrap SE and the minimum detectable effect at
  80% power (MDE80 = 2.80 x SE). For poly_noncpg the values are also given relative to the intronic variant
  rate of the TE group (--base-rate GROUP=rate, or --comp-dir with 14_composition.py *_paired_bins.tsv):
  bound_pct = largest TE-specific REDUCTION of variant frequency excluded at 95% (= -ci_lo / rate x 100).
Orientation contrast: sense - antisense, resampling genes jointly (genes may contribute both orientations).
Functional-copy table: if a fraction phi of exonic copies carries a domain covering a fraction d of the element
with a per-site variant reduction r, the whole-element reduction is phi*d*r, so the bound excludes
phi > bound / (d*r) (95%), and the test detects phi > MDE80 / (d*r) with 80% power.

Outputs: <out>_bounds.tsv, <out>_orientation.tsv, <out>_functional.tsv
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd

from ridl_utils import log


def boot_means(x, w, n_boot, rng):
    n = len(x)
    k = rng.integers(0, n, size=(n_boot, n))
    um = x[k].mean(axis=1)
    wm = (x[k] * w[k]).sum(axis=1) / w[k].sum(axis=1)
    return um, wm


def base_rates(args):
    rates = {}
    for item in args.base_rate:
        g, v = item.split("=")
        rates[g] = float(v)
    if args.comp_dir:
        for f in glob.glob(os.path.join(args.comp_dir, "te_*_paired_bins.tsv")):
            b = pd.read_csv(f, sep="\t")
            if "poly_c_noncpg_intron" not in b or "units_intron" not in b:
                continue
            name = os.path.basename(f)[3:-len("_paired_bins.tsv")]
            r = (b["poly_c_noncpg_intron"] * b["units_intron"]).sum() / b["units_intron"].sum()
            for g in (name, name + "*"):
                rates.setdefault(g, r)
    return rates


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes", action="append", required=True, help="label=path to *_genes.tsv (repeatable)")
    ap.add_argument("--comp-dir", default=None, help="directory with te_<group>_paired_bins.tsv (step 14)")
    ap.add_argument("--base-rate", action="append", default=[], help="GROUP=intronic non-CpG variant rate")
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--domain-fracs", default="0.05,0.1,0.2,0.33,0.5,1")
    ap.add_argument("--site-reductions", default="0.25,0.5,1")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    rates = base_rates(a)
    log("intronic variant rates: " + ", ".join(f"{k} {v:.4f}" for k, v in sorted(rates.items())))
    rows, orows = [], []
    for item in a.genes:
        label, path = item.split("=", 1)
        g = pd.read_csv(path, sep="\t")
        layers = [c for c in g.columns if c.startswith("excess_") or c.startswith("diff_")]
        g["_w"] = g["n_exon"] * g["n_intron"] / (g["n_exon"] + g["n_intron"])
        for (sf, ori), s in g.groupby(["subfamily", "orientation"]):
            for layer in layers:
                d = s[np.isfinite(s[layer])]
                if len(d) < 15:
                    continue
                x, w = d[layer].to_numpy(float), d["_w"].to_numpy(float)
                um, wm = boot_means(x, w, a.n_boot, rng)
                r = {"table": label, "subfamily": sf, "orientation": ori, "layer": layer, "genes": len(d),
                     "mean": x.mean(), "ci_lo": np.percentile(um, 2.5), "ci_hi": np.percentile(um, 97.5),
                     "se": um.std(ddof=1), "wmean": (x * w).sum() / w.sum(),
                     "w_ci_lo": np.percentile(wm, 2.5), "w_ci_hi": np.percentile(wm, 97.5), "w_se": wm.std(ddof=1)}
                r["mde80"] = 2.80 * r["se"]
                r["w_mde80"] = 2.80 * r["w_se"]
                rate = rates.get(sf)
                if layer.endswith("poly_noncpg") and rate:
                    r["intron_rate"] = rate
                    for k in ("mean", "ci_lo", "ci_hi", "mde80", "wmean", "w_ci_lo", "w_ci_hi", "w_mde80"):
                        r[f"{k}_pct"] = 100 * r[k] / rate
                    r["bound_pct"] = max(0.0, -r["ci_lo_pct"])
                    r["w_bound_pct"] = max(0.0, -r["w_ci_lo_pct"])
                rows.append(r)
            # orientation contrast (sense - antisense), genes resampled jointly
        for layer in layers:
            for sf, s in g.groupby("subfamily"):
                p = s.pivot_table(index="gene_id", columns="orientation", values=layer, aggfunc="first")
                if not {"sense", "antisense"} <= set(p.columns):
                    continue
                genes = p.index.to_numpy()
                vs, va = p["sense"].to_numpy(float), p["antisense"].to_numpy(float)
                if np.isfinite(vs).sum() < 15 or np.isfinite(va).sum() < 15:
                    continue
                obs = np.nanmean(vs) - np.nanmean(va)
                k = rng.integers(0, len(genes), size=(a.n_boot, len(genes)))
                with np.errstate(all="ignore"):
                    bs = np.nanmean(vs[k], axis=1) - np.nanmean(va[k], axis=1)
                p2 = 2 * min((1 + np.sum(bs <= 0)) / (1 + len(bs)), (1 + np.sum(bs >= 0)) / (1 + len(bs)))
                o = {"table": label, "subfamily": sf, "layer": layer, "genes_sense": int(np.isfinite(vs).sum()),
                     "genes_antisense": int(np.isfinite(va).sum()), "sense_minus_antisense": obs,
                     "ci_lo": np.nanpercentile(bs, 2.5), "ci_hi": np.nanpercentile(bs, 97.5), "p_boot": min(1.0, p2)}
                rate = rates.get(sf)
                if layer.endswith("poly_noncpg") and rate:
                    for kk in ("sense_minus_antisense", "ci_lo", "ci_hi"):
                        o[f"{kk}_pct"] = 100 * o[kk] / rate
                orows.append(o)
    B = pd.DataFrame(rows)
    O = pd.DataFrame(orows)
    B.to_csv(f"{a.out}_bounds.tsv", sep="\t", index=False)
    O.to_csv(f"{a.out}_orientation.tsv", sep="\t", index=False)

    frows = []
    dfr = [float(x) for x in a.domain_fracs.split(",")]
    srr = [float(x) for x in a.site_reductions.split(",")]
    sel = B[B["layer"].eq("excess_poly_noncpg") & B["orientation"].eq("pooled")] if "bound_pct" in B else B.iloc[:0]
    for r in sel.itertuples(index=False):
        for dd in dfr:
            for rr in srr:
                frows.append({"table": r.table, "subfamily": r.subfamily, "domain_fraction": dd, "site_reduction": rr,
                              "bound_pct": r.bound_pct, "mde80_pct": abs(r.mde80_pct),
                              "max_functional_fraction_95": min(1.0, r.bound_pct / 100 / (dd * rr)),
                              "detectable_fraction_80": min(1.0, abs(r.mde80_pct) / 100 / (dd * rr))})
    F = pd.DataFrame(frows)
    F.to_csv(f"{a.out}_functional.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    show = [c for c in ["table", "subfamily", "orientation", "layer", "genes", "mean", "ci_lo", "ci_hi", "bound_pct",
                        "w_bound_pct", "mde80_pct"] if c in B.columns]
    log(B[B["layer"].str.startswith("excess_")][show].round(4).to_string(index=False))
    if len(O):
        log(O.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
