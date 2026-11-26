#!/usr/bin/env python3
"""Step 6. Constraint profiles along the TE consensus: foreground vs background compartment.

Layers
  population  Mantel-Haenszel log OR of 'site is polymorphic in gnomAD' (fg vs bg), stratified by
              trinucleotide (TE orientation) x CpG status. Negative = fewer variants in fg = constraint.
  singletons  difference in the fraction of singletons among polymorphic sites (fg - bg), within strata;
              positive = rarer alleles in fg = constraint (MAPS-like, less sensitive to mutation rate).
  phylop_*    difference of mean phyloP (fg - bg), stratified by --depth-col quintiles if given;
              positive = more conserved in fg.
For each layer: per-bin and whole-element statistics, permutation p-values (unit labels shuffled,
--perm-scheme global|within_gene), cluster-bootstrap 95% CIs (genes resampled; intergenic units
are their own clusters), BH q-values across all bins of all comparisons, peaks (>= --min-run
consecutive bins with q < --q and the same sign), sense-vs-antisense interaction, and plots.

Centered profiles (columns *_centered, peaks_centered.tsv): each window's effect MINUS the median effect
of all windows of that element, tested with the same permutations/bootstrap. This asks "which segment differs from the rest of
the element" and is robust to uniform offsets (regional mutation rate, exon-wide effects).
--paired-genes keeps only genes that contribute both foreground and background units (within-gene
comparison; needs a gene-based background such as intron).
--gene-table adds per-gene columns (TSV with gene_id without version, e.g. from 11_gene_properties.py)
before --query is applied, so the query can use them (e.g. "stab_t == 3").

Typical runs
  main:          --fg exon --bg intron
  pre-mRNA:      --fg intron --bg intergenic
  null control:  --fg pseudoexon --bg intron   (annotation from step 7)
  sensitivity:   --query "dist_splice >= 50"  /  --query "exon_class == 'mono'"  /  "max_tpm >= 5"
"""
import argparse
import os

import numpy as np
import pandas as pd
from scipy.stats import norm

from ridl_stats import (Agg, bootstrap_weights, call_peaks, load_annot, mh_logor, perm_pvalue,
                        permuted_labels, select_comparison, strat_mean_diff)
from ridl_utils import bh_fdr, log


def pop_stats(agg, label, w):
    S = agg.sums(label, w)
    a, c = S["gn_poly"]
    n1, n0 = S["_n"]
    nb, ns = agg.n_bins, agg.n_strata
    lor_bin, lor_all = mh_logor(a, n1, c, n0, nb, ns)
    s1, s0 = S["gn_singleton"]
    # singleton fraction among polymorphic sites, fg - bg, compared WITHIN strata (context x CpG x ...),
    # because e.g. CpG variants recur and are less often singletons
    sd_bin, sd_all = strat_mean_diff(s1, a, s0, c, nb, ns)
    return {"population": (lor_bin, lor_all), "singletons": (sd_bin, sd_all)}


def cont_stats(agg, label, w, col, name):
    S = agg.sums(label, w)
    s1, s0 = S[col]
    n1, n0 = S["_n"]
    return {name: strat_mean_diff(s1, n1, s0, n0, agg.n_bins, agg.n_strata)}


def run_layerset(agg, fn, n_perm, n_boot, scheme, rng):
    ones = np.ones(len(agg.units))
    obs = fn(agg, agg.label, ones)
    perm = {k: ([], []) for k in obs}
    for _ in range(n_perm):
        r = fn(agg, permuted_labels(agg, rng, scheme), ones)
        for k in r:
            perm[k][0].append(r[k][0])
            perm[k][1].append(r[k][1])
    boot = {k: ([], []) for k in obs}
    for _ in range(n_boot):
        r = fn(agg, agg.label, bootstrap_weights(agg, rng))
        for k in r:
            boot[k][0].append(r[k][0])
            boot[k][1].append(r[k][1])
    res = {}
    for k in obs:
        pb = perm_pvalue(obs[k][0], np.array(perm[k][0])) if n_perm else np.full(agg.n_bins, np.nan)
        pa = perm_pvalue(obs[k][1], np.array(perm[k][1])) if n_perm else np.nan
        if n_boot:
            bb = np.array(boot[k][0])
            lo, hi = np.nanpercentile(bb, 2.5, axis=0), np.nanpercentile(bb, 97.5, axis=0)
            ba = np.array(boot[k][1], dtype=float)
            alo, ahi = np.nanpercentile(ba, 2.5), np.nanpercentile(ba, 97.5)
            with np.errstate(divide="ignore", invalid="ignore"):
                pz = 2 * norm.sf(np.abs(obs[k][0]) / np.nanstd(bb, axis=0))
                pz_all = float(2 * norm.sf(abs(obs[k][1]) / np.nanstd(ba)))
        else:
            lo = hi = pz = np.full(agg.n_bins, np.nan)
            alo = ahi = pz_all = np.nan
        # centered = window effect minus the MEDIAN window effect (robust to a few strongly deviating windows)
        with np.errstate(all="ignore"):
            cen = obs[k][0] - np.nanmedian(obs[k][0])
        if n_perm:
            pb_ = np.array(perm[k][0])
            with np.errstate(all="ignore"):
                pc = perm_pvalue(cen, pb_ - np.nanmedian(pb_, axis=1)[:, None])
        else:
            pc = np.full(agg.n_bins, np.nan)
        if n_boot:
            bb_ = np.array(boot[k][0])
            with np.errstate(all="ignore"):
                bc = bb_ - np.nanmedian(bb_, axis=1)[:, None]
            clo, chi = np.nanpercentile(bc, 2.5, axis=0), np.nanpercentile(bc, 97.5, axis=0)
            with np.errstate(divide="ignore", invalid="ignore"):
                pzc = 2 * norm.sf(np.abs(cen) / np.nanstd(bc, axis=0))
        else:
            clo = chi = pzc = np.full(agg.n_bins, np.nan)
        res[k] = dict(bin=obs[k][0], p=pb, pz=pz, lo=lo, hi=hi, all=obs[k][1], p_all=float(pa),
                      pz_all=pz_all, all_lo=alo, all_hi=ahi, boot=np.array(boot[k][0]) if n_boot else None,
                      cen=cen, p_cen=pc, pz_cen=pzc, cen_lo=clo, cen_hi=chi)
    return res


def strata_codes(d, cols):
    if not cols:
        return None
    key = d[cols[0]].astype(str)
    for c in cols[1:]:
        key = key + "|" + d[c].astype(str)
    return pd.factorize(key)[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annot", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--fg", default="exon")
    ap.add_argument("--bg", default="intron")
    ap.add_argument("--bin-size", type=int, default=10)
    ap.add_argument("--query", default=None)
    ap.add_argument("--subfamilies", default=None, help="comma list (default: all present)")
    ap.add_argument("--pop-strata", default="tri,cpg", help="columns defining MH strata ('' = none)")
    ap.add_argument("--depth-col", default=None, help="stratify phyloP by quintiles of this column")
    ap.add_argument("--cont-strata", default="", help="categorical columns stratifying phyloP layers, e.g. subfamily_orig")
    ap.add_argument("--n-perm", type=int, default=1000)
    ap.add_argument("--n-boot", type=int, default=200)
    ap.add_argument("--perm-scheme", choices=["global", "within_gene"], default="global")
    ap.add_argument("--min-units", type=int, default=20, help="min units per bin in fg and bg")
    ap.add_argument("--p-source", choices=["perm", "boot_z"], default="perm",
                    help="p-values used for BH q: permutation (floor 1/(n_perm+1)) or bootstrap z")
    ap.add_argument("--q", type=float, default=0.1)
    ap.add_argument("--min-run", type=int, default=2)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--no-plots", action="store_true")
    ap.add_argument("--gene-table", default=None, help="TSV with gene_id (no version) + per-gene columns")
    ap.add_argument("--paired-genes", action="store_true",
                    help="keep only genes contributing both fg and bg units (within-gene comparison)")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    rng = np.random.default_rng(a.seed)

    import pyarrow.parquet as pq
    first = a.annot if not os.path.isdir(a.annot) else sorted(
        os.path.join(a.annot, f) for f in os.listdir(a.annot) if f.endswith(".parquet"))[0]
    schema = pq.read_schema(first).names
    need = ["unit_id", "subfamily", "compartment", "orientation", "gene_id", "cons_pos", "gn_callable",
            "gn_poly", "gn_singleton"] + [c for c in schema if c.startswith("phylop_")]
    need += [c for c in a.pop_strata.split(",") if c]
    need += [c for c in a.cont_strata.split(",") if c]
    if a.depth_col:
        need.append(a.depth_col)
    if a.gene_table:
        need += [c for c in schema if a.query and c in a.query]
        df = load_annot(a.annot, columns=need, query=None,
                        subfamilies=a.subfamilies.split(",") if a.subfamilies else None)
        gt = pd.read_csv(a.gene_table, sep="\t")
        gt["gene_id"] = gt["gene_id"].astype(str).str.split(".").str[0]
        gt = gt.drop_duplicates("gene_id").rename(columns={"gene_id": "_gid"})
        df["_gid"] = df["gene_id"].astype(str).str.split(".").str[0]
        df = df.merge(gt, on="_gid", how="left").drop(columns="_gid")
        log(f"  gene table: {gt.shape[1] - 1} columns; bases with a matched gene: "
            f"{df[gt.columns[1]].notna().mean():.1%}")
        if a.query:
            df = df.query(a.query, engine="python")
    else:
        df = load_annot(a.annot, columns=need, query=a.query,
                        subfamilies=a.subfamilies.split(",") if a.subfamilies else None)
    df = df[df["compartment"].isin([a.fg, a.bg])].copy()
    df["bin"] = (df["cons_pos"] - 1) // a.bin_size
    phylop_cols = [c for c in df.columns if c.startswith("phylop_")]
    has_pop = "gn_poly" in df.columns
    pop_strata = [c for c in a.pop_strata.split(",") if c and c in df.columns]
    cont_strata = [c for c in a.cont_strata.split(",") if c and c in df.columns]
    for c in set(a.pop_strata.split(",") + a.cont_strata.split(",")) - set(df.columns) - {""}:
        log(f"  WARNING: stratum column '{c}' not in the annotation -- ignored")
    sfs = a.subfamilies.split(",") if a.subfamilies else sorted(df["subfamily"].unique())
    log(f"layers: {'population, singletons, ' if has_pop else ''}{', '.join(phylop_cols)}")

    rows, whole, boots = [], [], {}
    for sf in sfs:
        dsf = df[df["subfamily"] == sf]
        if dsf.empty:
            continue
        nb = int(dsf["bin"].max()) + 1
        oris = sorted(dsf.loc[dsf["compartment"] == a.fg, "orientation"].unique())
        for ori in oris:
            d = select_comparison(dsf, sf, ori, a.fg, a.bg)
            if a.paired_genes:
                g1 = set(d.loc[d["_lab"] == 1, "gene_id"]) - {"."}
                g0 = set(d.loc[d["_lab"] == 0, "gene_id"]) - {"."}
                if not g0:
                    log("  WARNING: --paired-genes needs a gene-based background (e.g. intron); ignored")
                else:
                    d = d[d["gene_id"].isin(g1 & g0)]
                    log(f"  {sf} {ori}: paired genes {len(g1 & g0)} (of {len(g1)} fg genes)")
            if d["_lab"].nunique() < 2:
                continue
            log(f"  {sf} {ori}: fg units {d.loc[d['_lab'] == 1, 'unit_id'].nunique()}, "
                f"bg units {d.loc[d['_lab'] == 0, 'unit_id'].nunique()}")
            results = {}
            unit_counts = None
            if has_pop:
                dp = d[d["gn_callable"]].copy()
                dp["_s"] = strata_codes(dp, pop_strata) if pop_strata else 0
                for c in ("gn_poly", "gn_singleton"):
                    dp[c] = dp[c].astype(float)
                agg = Agg(dp, ["gn_poly", "gn_singleton"], strata_col="_s", n_bins=nb)
                if a.perm_scheme == "within_gene":
                    mixed = pd.Series(agg.label).groupby(agg.unit_cluster).transform("nunique") > 1
                    if mixed.mean() < 0.5:
                        log(f"    WARNING: only {mixed.mean():.0%} of units sit in genes with both fg and bg "
                            "units; within_gene permutations have little power here (use global)")
                unit_counts = agg.units_per_bin(agg.label)
                results.update(run_layerset(agg, pop_stats, a.n_perm, a.n_boot, a.perm_scheme, rng))
            for col in phylop_cols:
                dc = d[np.isfinite(d[col])].copy()
                if a.depth_col and a.depth_col in dc.columns:
                    nodepth = ~np.isfinite(dc[a.depth_col].astype(float))
                    if nodepth.any():
                        log(f"    {col}: {nodepth.mean():.1%} bases without {a.depth_col} dropped")
                        dc = dc[~nodepth].copy()
                    dc["_q"] = pd.qcut(dc[a.depth_col].rank(method="first"), 5, labels=False)
                    scols = ["_q"] + cont_strata
                else:
                    scols = cont_strata
                if scols:
                    dc["_s"] = strata_codes(dc, scols)
                    sc = "_s"
                else:
                    sc = None
                aggc = Agg(dc, [col], strata_col=sc, n_bins=nb)
                if unit_counts is None:
                    unit_counts = aggc.units_per_bin(aggc.label)
                fn = (lambda ag, lab, w, _c=col: cont_stats(ag, lab, w, _c, _c))
                results.update(run_layerset(aggc, fn, a.n_perm, a.n_boot, a.perm_scheme, rng))
            u1, u0 = unit_counts
            enough = (u1 >= a.min_units) & (u0 >= a.min_units)
            for layer, r in results.items():
                boots[(sf, ori, layer)] = r["boot"]
                whole.append({"layer": layer, "subfamily": sf, "orientation": ori, "effect": r["all"],
                              "ci_lo": r["all_lo"], "ci_hi": r["all_hi"], "p_perm": r["p_all"],
                              "p_boot_z": r["pz_all"],
                              "units_fg": int((d["_lab"] == 1).groupby(d["unit_id"]).any().sum()),
                              "units_bg": int((d["_lab"] == 0).groupby(d["unit_id"]).any().sum())})
                for b in range(nb):
                    ok = enough[b]
                    rows.append({"layer": layer, "subfamily": sf, "orientation": ori, "bin": b,
                                 "cons_from": b * a.bin_size + 1, "cons_to": (b + 1) * a.bin_size,
                                 "units_fg": int(u1[b]), "units_bg": int(u0[b]),
                                 "effect": r["bin"][b] if ok else np.nan,
                                 "ci_lo": r["lo"][b] if ok else np.nan, "ci_hi": r["hi"][b] if ok else np.nan,
                                 "p_perm": r["p"][b] if ok else np.nan,
                                 "p_boot_z": r["pz"][b] if ok else np.nan,
                                 "effect_centered": r["cen"][b] if ok else np.nan,
                                 "cen_lo": r["cen_lo"][b] if ok else np.nan,
                                 "cen_hi": r["cen_hi"][b] if ok else np.nan,
                                 "p_perm_centered": r["p_cen"][b] if ok else np.nan,
                                 "p_boot_z_centered": r["pz_cen"][b] if ok else np.nan})

    prof = pd.DataFrame(rows)
    whole = pd.DataFrame(whole)
    if prof.empty:
        raise SystemExit("nothing to test -- check --fg/--bg and the annotation")
    pcol = "p_perm" if a.p_source == "perm" else "p_boot_z"
    prof["q"] = np.nan
    prof["q_centered"] = np.nan
    for layer, g in prof.groupby("layer"):
        prof.loc[g.index, "q"] = bh_fdr(g[pcol].to_numpy())
        prof.loc[g.index, "q_centered"] = bh_fdr(g[pcol + "_centered"].to_numpy())
    whole["q"] = bh_fdr(whole[pcol].to_numpy())
    prof.to_csv(os.path.join(a.out_dir, "profiles.tsv"), sep="\t", index=False)
    whole.to_csv(os.path.join(a.out_dir, "whole_element.tsv"), sep="\t", index=False)

    peaks = []
    for (layer, sf, ori), g in prof.groupby(["layer", "subfamily", "orientation"]):
        for pk in call_peaks(g, "effect", "q", a.q, a.min_run):
            peaks.append({"layer": layer, "subfamily": sf, "orientation": ori, **pk})
    peaks = pd.DataFrame(peaks)
    peaks.to_csv(os.path.join(a.out_dir, "peaks.tsv"), sep="\t", index=False)
    cpeaks = []
    for (layer, sf, ori), g in prof.groupby(["layer", "subfamily", "orientation"]):
        for pk in call_peaks(g, "effect_centered", "q_centered", a.q, a.min_run):
            cpeaks.append({"layer": layer, "subfamily": sf, "orientation": ori, **pk})
    cpeaks = pd.DataFrame(cpeaks)
    cpeaks.to_csv(os.path.join(a.out_dir, "peaks_centered.tsv"), sep="\t", index=False)

    # sense vs antisense interaction (bootstrap of the difference; independent resamples of the two
    # comparisons -- conservative when they share background units)
    inter = []
    for (layer, sf), g in prof.groupby(["layer", "subfamily"]):
        if set(g["orientation"]) >= {"sense", "antisense"}:
            bs, ba = boots.get((sf, "sense", layer)), boots.get((sf, "antisense", layer))
            gs = g[g["orientation"] == "sense"].set_index("bin")["effect"]
            ga = g[g["orientation"] == "antisense"].set_index("bin")["effect"]
            nbins = min(len(gs), len(ga))
            for b in range(nbins):
                diff = gs.iloc[b] - ga.iloc[b]
                p = np.nan
                lo = hi = np.nan
                if bs is not None and ba is not None and np.isfinite(diff):
                    k = min(len(bs), len(ba))
                    dd = bs[:k, b] - ba[:k, b]
                    dd = dd[np.isfinite(dd)]
                    if len(dd):
                        lo, hi = np.percentile(dd, [2.5, 97.5])
                        p = min(1.0, 2 * min((1 + np.sum(dd <= 0)) / (1 + len(dd)),
                                             (1 + np.sum(dd >= 0)) / (1 + len(dd))))
                inter.append({"layer": layer, "subfamily": sf, "bin": b, "cons_from": b * a.bin_size + 1,
                              "cons_to": (b + 1) * a.bin_size, "sense_minus_antisense": diff,
                              "ci_lo": lo, "ci_hi": hi, "p_boot": p})
    inter = pd.DataFrame(inter)
    if not inter.empty:
        inter["q"] = np.nan
        for layer, g in inter.groupby("layer"):
            inter.loc[g.index, "q"] = bh_fdr(g["p_boot"].to_numpy())
    inter.to_csv(os.path.join(a.out_dir, "strand_interaction.tsv"), sep="\t", index=False)

    log("whole-element results:\n" + whole.to_string(index=False))
    log(f"peaks: {len(peaks)}")
    if len(peaks):
        log(peaks.to_string(index=False))
    log(f"centered peaks (segments differing from the rest of the element): {len(cpeaks)}")
    if len(cpeaks):
        log(cpeaks.to_string(index=False))
    if not a.no_plots:
        plot_profiles(prof, a, whole)


def plot_profiles(prof, a, whole=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ylab = {"population": "MH log OR polymorphic\n(fg vs bg; <0 = constraint)",
            "singletons": "singleton fraction\nfg - bg (>0 = constraint)"}
    for sf, g in prof.groupby("subfamily"):
        layers = sorted(g["layer"].unique(), key=lambda x: (x != "population", x != "singletons", x))
        oris = sorted(g["orientation"].unique())
        fig, axes = plt.subplots(len(layers) + 1, len(oris), figsize=(6 * len(oris), 2.4 * (len(layers) + 1)),
                                 squeeze=False, sharex=True)
        for j, ori in enumerate(oris):
            for i, layer in enumerate(layers):
                ax = axes[i, j]
                d = g[(g["layer"] == layer) & (g["orientation"] == ori)].sort_values("bin")
                x = (d["cons_from"] + d["cons_to"]) / 2
                ax.axhline(0, color="grey", lw=0.8)
                ax.fill_between(x, d["ci_lo"], d["ci_hi"], color="tab:blue", alpha=0.2, lw=0)
                ax.plot(x, d["effect"], color="tab:blue", lw=1.2)
                sig = d["q"] < a.q
                ax.scatter(x[sig], d.loc[sig, "effect"], color="tab:red", s=14, zorder=3)
                if "effect_centered" in d.columns:
                    med = (d["effect"] - d["effect_centered"]).dropna()
                    if len(med):
                        ax.axhline(med.iloc[0], color="tab:blue", ls=":", lw=1)
                if "q_centered" in d.columns:
                    cs = d["q_centered"] < a.q
                    ax.scatter(x[cs], d.loc[cs, "effect"], facecolors="none", edgecolors="black", s=46, zorder=4)
                ax.set_ylabel(ylab.get(layer, f"{layer}\nfg - bg (>0 = constraint)"), fontsize=8)
                if i == 0:
                    ax.set_title(f"{sf}  {a.fg} vs {a.bg}  ({ori})", fontsize=10)
            ax = axes[-1, j]
            d = g[(g["layer"] == layers[0]) & (g["orientation"] == ori)].sort_values("bin")
            x = (d["cons_from"] + d["cons_to"]) / 2
            ax.plot(x, d["units_fg"], label=f"{a.fg} units")
            ax.plot(x, d["units_bg"], label=f"{a.bg} units")
            ax.axhline(a.min_units, color="grey", ls=":", lw=0.8)
            ax.set_yscale("log")
            ax.set_xlabel("consensus position")
            ax.legend(fontsize=7)
        fig.suptitle(f"red dot: window differs from 0 (q < {a.q}); black ring: window differs from the "
                     f"element's median window (centered, q < {a.q}); dotted line: median window", fontsize=8, y=0.995)
        fig.tight_layout(rect=(0, 0, 1, 0.985))
        fig.savefig(os.path.join(a.out_dir, f"profile_{sf}.png"), dpi=130)
        plt.close(fig)


if __name__ == "__main__":
    main()
