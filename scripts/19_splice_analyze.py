#!/usr/bin/env python3
"""Splice-site constraint: are the GT/AG dinucleotides of lncRNA splice sites under selection, and those
that lie inside transposable elements in particular?

Input: positions.tsv.gz from 18_splice_prepare.py. Every position is annotated with gnomAD (polymorphic,
singleton, callable from the all-sites allele-number table), phyloP tracks and sequence context. Then:
  flank comparisons   lnc_nonTE_site vs lnc_nonTE_flank and pc_site vs pc_flank, stratified by
                      trinucleotide x CpG; clusters for the bootstrap = splice sites
  matched comparisons lnc_TE_site vs lnc_TE_ctrl (same subfamily, consensus position and base in unused
                      copies), stratified by that match; clusters = strata. Also by TE family, site type,
                      orientation and for the antisense-Alu acceptor hotspots.
Effects: Mantel-Haenszel log odds ratio of a site being polymorphic (negative = fewer variants at splice
sites = constraint), stratified difference in phyloP (positive = more conserved), stratified difference in
the share of singletons among variants (positive = rarer variants). 95% intervals by cluster bootstrap.

Outputs: <out>_positions.tsv.gz (annotated), <out>_stats.tsv
"""
import argparse

import numpy as np
import pandas as pd
import pysam
import pyBigWig

from ridl_utils import log, merge_intervals

COMP = str.maketrans("ACGTN", "TGCAN")


def gnomad_lookup(pos, vcfs):
    ac = pd.Series(0, index=pos.index, dtype=np.int64)
    files = {}
    for p in vcfs:
        try:
            cs = pysam.TabixFile(p).contigs
        except (OSError, ValueError):
            cs = list(pysam.VariantFile(p).header.contigs)
        for c in cs:
            files.setdefault(c, []).append(p)
    for c, g in pos.groupby("chrom"):
        found = {}
        iv = merge_intervals(pd.DataFrame({"chrom": c, "start": g["pos"] - 200, "end": g["pos"] + 200}))
        for p in files.get(c, []):
            vf = pysam.VariantFile(p)
            for r in iv.itertuples(index=False):
                for rec in vf.fetch(c, max(0, int(r.start)), int(r.end)):
                    a_ = rec.info.get("AC")
                    n = int(sum(a_)) if isinstance(a_, tuple) else int(a_ or 0)
                    found[rec.pos - 1] = max(found.get(rec.pos - 1, 0), n)
        ac.loc[g.index] = g["pos"].map(found).fillna(0).astype(np.int64).to_numpy()
    return ac


def an_lookup(pos, an_tsv):
    tb = pysam.TabixFile(an_tsv)
    an = pd.Series(np.nan, index=pos.index)
    for c, g in pos.groupby("chrom"):
        found = {}
        iv = merge_intervals(pd.DataFrame({"chrom": c, "start": g["pos"] - 5, "end": g["pos"] + 5}))
        for r in iv.itertuples(index=False):
            try:
                for line in tb.fetch(c, max(0, int(r.start)), int(r.end)):
                    f = line.split("\t")
                    found[int(f[1]) - 1] = float(f[2])
            except ValueError:
                pass
        an.loc[g.index] = g["pos"].map(found).to_numpy()
    return an


def phylop_lookup(pos, path):
    bw = pyBigWig.open(path)
    v = pd.Series(np.nan, index=pos.index)
    for c, g in pos.groupby("chrom"):
        vals = []
        for p in g["pos"].to_numpy():
            x = bw.values(c, int(p), int(p) + 1)[0]
            vals.append(np.nan if x is None else x)
        v.loc[g.index] = vals
    return v


def compare(df, lab, stratum, cluster, metrics, n_boot, rng):
    """df rows: lab (1 site / 0 control). Returns dict of estimates and bootstrap CIs."""
    d = df.copy()
    d["_s"] = pd.factorize(d[stratum])[0]
    d["_c"] = pd.factorize(d[cluster])[0]
    ns, nc = d["_s"].max() + 1, d["_c"].max() + 1
    s_idx, c_idx = d["_s"].to_numpy(), d["_c"].to_numpy()
    L = d[lab].to_numpy().astype(bool)
    call = d["callable"].to_numpy().astype(bool)
    poly = d["poly"].to_numpy().astype(bool) & call
    sing = d["singleton"].to_numpy().astype(bool) & call

    def est(w):
        out = {}
        bc = lambda m: np.bincount(s_idx, weights=w * m, minlength=ns)
        A, B = bc(L & call & poly), bc(L & call & ~poly)
        C, D = bc(~L & call & poly), bc(~L & call & ~poly)
        n = A + B + C + D
        with np.errstate(divide="ignore", invalid="ignore"):
            R, S = np.nansum(np.where(n > 0, A * D / n, 0)), np.nansum(np.where(n > 0, B * C / n, 0))
            out["mh_log_or"] = np.log(R / S) if R > 0 and S > 0 else np.nan
        # singleton share among variants, stratified difference
        s1, p1, s0, p0 = bc(L & sing), bc(L & poly), bc(~L & sing), bc(~L & poly)
        h = np.where((p1 > 0) & (p0 > 0), p1 * p0 / np.maximum(p1 + p0, 1e-9), 0)
        with np.errstate(divide="ignore", invalid="ignore"):
            out["singleton_diff"] = (np.nansum(h * (np.where(p1 > 0, s1 / p1, 0) - np.where(p0 > 0, s0 / p0, 0)))
                                     / h.sum()) if h.sum() > 0 else np.nan
        for m in metrics:
            v = d[m].to_numpy()
            ok = np.isfinite(v)
            x1 = np.bincount(s_idx, weights=w * np.where(L & ok, v, 0), minlength=ns)
            x0 = np.bincount(s_idx, weights=w * np.where(~L & ok, v, 0), minlength=ns)
            n1, n0 = bc(L & ok), bc(~L & ok)
            hh = np.where((n1 > 0) & (n0 > 0), n1 * n0 / np.maximum(n1 + n0, 1e-9), 0)
            with np.errstate(divide="ignore", invalid="ignore"):
                out[f"{m}_diff"] = (np.nansum(hh * (np.where(n1 > 0, x1 / n1, 0) - np.where(n0 > 0, x0 / n0, 0)))
                                    / hh.sum()) if hh.sum() > 0 else np.nan
        return out

    obs = est(np.ones(len(d)))
    boots = []
    for _ in range(n_boot):
        cnt = np.bincount(rng.integers(0, nc, nc), minlength=nc).astype(float)
        boots.append(est(cnt[c_idx]))
    res = {"clusters": nc, "bases_site": int(L.sum()), "bases_ctrl": int((~L).sum()),
           "poly_site": poly[L].sum() / max(call[L].sum(), 1), "poly_ctrl": poly[~L].sum() / max(call[~L].sum(), 1)}
    for m in metrics:
        res[f"{m}_site"] = np.nanmean(d.loc[L, m])
        res[f"{m}_ctrl"] = np.nanmean(d.loc[~L, m])
    for k, v in obs.items():
        b = np.array([x[k] for x in boots], dtype=float)
        res[k] = v
        res[f"{k}_lo"], res[f"{k}_hi"] = (np.nanpercentile(b, [2.5, 97.5]) if np.isfinite(b).any() else (np.nan, np.nan))
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--positions", required=True)
    ap.add_argument("--gnomad-vcf", action="append", default=[])
    ap.add_argument("--an-tsv", required=True)
    ap.add_argument("--an-min-frac", type=float, default=0.8)
    ap.add_argument("--phylop", action="append", default=[], help="name=path.bw")
    ap.add_argument("--genome-fasta", required=True)
    ap.add_argument("--n-boot", type=int, default=200)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    pos = pd.read_csv(a.positions, sep="\t", low_memory=False)
    log(f"positions: {len(pos)}")

    fa = pysam.FastaFile(a.genome_fasta)
    tri, cpg = [], []
    for c, p in zip(pos["chrom"], pos["pos"]):
        s = fa.fetch(c, max(0, int(p) - 1), int(p) + 2).upper()
        tri.append(s)
        cpg.append(s[1:3] == "CG" or s[0:2] == "CG")
    pos["tri"], pos["cpg"] = tri, cpg
    pos["ac"] = gnomad_lookup(pos, a.gnomad_vcf)
    pos["an"] = an_lookup(pos, a.an_tsv)
    an_max = pos.groupby("chrom")["an"].transform("max")
    pos["callable"] = pos["an"] >= a.an_min_frac * an_max
    pos["poly"] = pos["ac"] > 0
    pos["singleton"] = pos["ac"] == 1
    metrics = []
    for item in a.phylop:
        n, p = item.split("=", 1)
        pos[f"phylop_{n}"] = phylop_lookup(pos, p)
        metrics.append(f"phylop_{n}")
    pos["ctx"] = pos["tri"] + "_" + pos["cpg"].astype(str)
    pos.to_csv(f"{a.out}_positions.tsv.gz", sep="\t", index=False)
    log("callable / polymorphic by set:\n" + pos.groupby("set").agg(
        n=("pos", "size"), callable=("callable", "mean"), poly=("poly", "mean")).round(4).to_string())

    rows = []
    for name, s_set, c_set in (("lncRNA non-TE splice sites vs intronic flank", "lnc_nonTE_site", "lnc_nonTE_flank"),
                               ("protein-coding splice sites vs intronic flank", "pc_site", "pc_flank")):
        d = pos[pos["set"].isin([s_set, c_set])].copy()
        if d.empty:
            continue
        d["lab"] = (d["set"] == s_set).astype(int)
        for typ in ("all", "donor", "acceptor"):
            dd = d if typ == "all" else d[d["site_type"] == typ]
            r = compare(dd, "lab", "ctx", "site_id", metrics, a.n_boot, rng)
            rows.append({"comparison": name, "subset": typ, **r})

    te = pos[pos["set"] == "lnc_TE_site"].copy()
    ctrl = pos[pos["set"] == "lnc_TE_ctrl"].copy()
    hot = ((te["rep_family"] == "Alu") & (te["te_orientation"] == "antisense") & (te["site_type"] == "acceptor")
           & (te["cons_pos"].between(266, 288) | te["cons_pos"].between(106, 124)))
    subsets = {"all": te.index, "Alu": te.index[te["rep_family"] == "Alu"], "MIR": te.index[te["rep_family"] == "MIR"],
               "L2": te.index[te["rep_family"] == "L2"], "donor": te.index[te["site_type"] == "donor"],
               "acceptor": te.index[te["site_type"] == "acceptor"],
               "TE antisense to transcript": te.index[te["te_orientation"] == "antisense"],
               "TE sense to transcript": te.index[te["te_orientation"] == "sense"],
               "antisense-Alu acceptor hotspots": te.index[hot]}
    for sub, idx in subsets.items():
        s = te.loc[idx]
        if s.empty:
            continue
        c = ctrl[ctrl["stratum"].isin(set(s["stratum"]))]
        d = pd.concat([s.assign(lab=1), c.assign(lab=0)], ignore_index=True)
        if c.empty:
            continue
        r = compare(d, "lab", "stratum", "stratum", metrics, a.n_boot, rng)
        rows.append({"comparison": "lncRNA TE-derived splice sites vs same positions in unused copies",
                     "subset": sub, **r})
    st = pd.DataFrame(rows)
    st.to_csv(f"{a.out}_stats.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    show = ["comparison", "subset", "clusters", "bases_site", "poly_site", "poly_ctrl", "mh_log_or", "mh_log_or_lo",
            "mh_log_or_hi"] + [c for c in st.columns if c.endswith("_diff")]
    log(st[show].round(4).to_string(index=False))


if __name__ == "__main__":
    main()
