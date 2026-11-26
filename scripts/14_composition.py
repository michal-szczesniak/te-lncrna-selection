#!/usr/bin/env python3
"""Composition check: are exon - intron differences explained by sequence composition (GC, CpG)?

Works on a step-4 annotation (typically 'nonTE' windows of equal length). For every unit (window or copy
part) it computes GC content, mean phyloP of each track, and the polymorphism rate at callable sites.
Exon and intron units are then compared
  (1) unmatched,
  (2) at the base level for non-CpG and for CpG positions separately,
  (3) GC-matched: units are binned into GC quintiles and intron means are re-weighted to the exon GC
      distribution (difference = sum over bins of exon share x (exon mean - intron mean)).
95% intervals come from resampling genes. With --paired-genes only genes with both exon and intron units
are used. If the GC-matched or non-CpG differences shrink towards zero, composition explains the gap.

Outputs: <out>_summary.tsv, <out>_bins.tsv
"""
import argparse

import numpy as np
import pandas as pd

from ridl_stats import load_annot
from ridl_utils import log


def gc_matched(u, col, bins):
    d = 0.0
    nex = u.loc[u["compartment"] == "exon", col].notna().sum()
    if nex == 0:
        return np.nan
    for b in range(bins):
        ub = u[u["gcbin"] == b]
        e = ub.loc[ub["compartment"] == "exon", col].dropna()
        i = ub.loc[ub["compartment"] == "intron", col].dropna()
        if len(e) == 0:
            continue
        if len(i) == 0:
            return np.nan
        d += len(e) / nex * (e.mean() - i.mean())
    return d


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annot", required=True)
    ap.add_argument("--subfamilies", default="nonTE")
    ap.add_argument("--query", default=None)
    ap.add_argument("--paired-genes", action="store_true")
    ap.add_argument("--bins", type=int, default=5)
    ap.add_argument("--n-boot", type=int, default=200)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)

    import glob, os
    import pyarrow.parquet as pq
    first = sorted(glob.glob(os.path.join(a.annot, "*.parquet")))[0] if os.path.isdir(a.annot) else a.annot
    schema = pq.read_schema(first).names
    pcols = [c for c in schema if c.startswith("phylop_")]
    cols = ["unit_id", "gene_id", "compartment", "subfamily", "ref", "cpg"] + pcols + \
        [c for c in ("gn_poly", "gn_callable") if c in schema]
    df = load_annot(a.annot, columns=cols, query=a.query, subfamilies=a.subfamilies.split(","))
    df = df[df["compartment"].isin(["exon", "intron"])].copy()
    df["gene_id"] = df["gene_id"].astype(str)
    if a.paired_genes:
        g = df.groupby("gene_id")["compartment"].nunique()
        df = df[df["gene_id"].isin(g[g == 2].index)]
    df["is_gc"] = df["ref"].astype(str).isin(["G", "C"])
    has_pop = "gn_poly" in df.columns
    if has_pop:
        df["poly_c"] = np.where(df["gn_callable"], df["gn_poly"].astype(float), np.nan)
        df["poly_c_noncpg"] = np.where(df["gn_callable"] & ~df["cpg"], df["gn_poly"].astype(float), np.nan)
    metrics = list(pcols) + (["poly_c", "poly_c_noncpg"] if has_pop else [])
    for c in pcols:
        df[f"{c}_noncpg"] = df[c].where(~df["cpg"])
    metrics += [f"{c}_noncpg" for c in pcols]

    keys = ["unit_id", "gene_id", "compartment"]
    u = df.groupby(keys, observed=True).agg(gc=("is_gc", "mean"), cpg=("cpg", "mean"),
                                           **{m: (m, "mean") for m in metrics}).reset_index()
    u["gcbin"] = pd.qcut(u["gc"].rank(method="first"), a.bins, labels=False)
    log(f"units: exon {int((u.compartment == 'exon').sum())}, intron {int((u.compartment == 'intron').sum())}; "
        f"median GC exon {u.loc[u.compartment == 'exon', 'gc'].median():.3f}, "
        f"intron {u.loc[u.compartment == 'intron', 'gc'].median():.3f}; "
        f"CpG share exon {u.loc[u.compartment == 'exon', 'cpg'].mean():.4f}, "
        f"intron {u.loc[u.compartment == 'intron', 'cpg'].mean():.4f}")

    genes = u["gene_id"].unique()
    gidx = {g: np.flatnonzero(u["gene_id"].to_numpy() == g) for g in genes}

    def stats(uu):
        out = {}
        for m in metrics:
            e = uu.loc[uu["compartment"] == "exon", m].dropna()
            i = uu.loc[uu["compartment"] == "intron", m].dropna()
            out[(m, "unmatched")] = e.mean() - i.mean()
            out[(m, "gc_matched")] = gc_matched(uu, m, a.bins)
        return out

    obs = stats(u)
    boots = []
    for _ in range(a.n_boot):
        pick = rng.choice(genes, len(genes), replace=True)
        boots.append(stats(u.iloc[np.concatenate([gidx[g] for g in pick])]))
    rows = []
    for (m, kind), v in obs.items():
        bv = np.array([b[(m, kind)] for b in boots], dtype=float)
        lo, hi = np.nanpercentile(bv, [2.5, 97.5]) if np.isfinite(bv).any() else (np.nan, np.nan)
        rows.append({"metric": m, "comparison": kind, "exon_minus_intron": v, "ci_lo": lo, "ci_hi": hi})
    s = pd.DataFrame(rows)
    s.to_csv(f"{a.out}_summary.tsv", sep="\t", index=False)

    brows = []
    for b in range(a.bins):
        ub = u[u["gcbin"] == b]
        r = {"gc_bin": b, "gc_range": f"{ub['gc'].min():.2f}-{ub['gc'].max():.2f}",
             "units_exon": int((ub.compartment == "exon").sum()), "units_intron": int((ub.compartment == "intron").sum())}
        for m in metrics:
            r[f"{m}_exon"] = ub.loc[ub.compartment == "exon", m].mean()
            r[f"{m}_intron"] = ub.loc[ub.compartment == "intron", m].mean()
        brows.append(r)
    pd.DataFrame(brows).to_csv(f"{a.out}_bins.tsv", sep="\t", index=False)
    pd.set_option("display.width", 200)
    log(s.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
