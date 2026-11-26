#!/usr/bin/env python3
"""Gene-level dose test: does the WITHIN-GENE exon - intron difference grow with RNA properties?

For each gene with enough exonic and intronic bases of a TE group (or of 'nonTE' windows), the script
computes exon - intron differences of mean phyloP (each phyloP column) and of the polymorphism rate at
callable non-CpG sites. It then correlates these per-gene differences with gene properties (Spearman),
with bootstrap 95% CIs over genes. Continuous correlations use every gene and are more powerful than
tertiles.

--baseline <per-gene table from a nonTE run>: also tests the EXCESS of the TE difference over the
difference of ordinary sequence in the same gene (TE diff - nonTE diff). An RNA-level TE effect predicts a
positive correlation of the phyloP excess with stability/maturity; an exon-wide effect does not.

Outputs: <out>_genes.tsv (per-gene differences), <out>_trend.tsv (correlations)
"""
import argparse

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from ridl_stats import load_annot
from ridl_utils import log


def per_gene(df, pcols, min_bases):
    rows = []
    for (sf, ori, gid), g in df.groupby(["subfamily", "orientation", "gene_id"], observed=True):
        ex, it = g[g["compartment"] == "exon"], g[g["compartment"] == "intron"]
        if len(ex) < min_bases or len(it) < min_bases:
            continue
        r = {"subfamily": sf, "orientation": ori, "gene_id": gid, "n_exon": len(ex), "n_intron": len(it)}
        for c in pcols:
            r[f"diff_{c}"] = ex[c].mean() - it[c].mean()
        if "gn_poly" in g.columns:
            e = ex[ex["gn_callable"] & ~ex["cpg"]]
            i = it[it["gn_callable"] & ~it["cpg"]]
            r["diff_poly_noncpg"] = (e["gn_poly"].mean() - i["gn_poly"].mean()) if len(e) and len(i) else np.nan
        rows.append(r)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annot", required=True)
    ap.add_argument("--subfamilies", required=True, help="comma list; each analysed separately")
    ap.add_argument("--gene-table", required=True, help="from 11_gene_properties.py")
    ap.add_argument("--props", default="stab,nuc,splice,synth,tau")
    ap.add_argument("--query", default=None)
    ap.add_argument("--min-bases", type=int, default=30)
    ap.add_argument("--n-boot", type=int, default=1000)
    ap.add_argument("--baseline", default=None, help="<out>_genes.tsv of a nonTE run")
    ap.add_argument("--pool-orientations", action="store_true", help="also analyse sense+antisense together")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)

    import glob, os
    import pyarrow.parquet as pq
    first = sorted(glob.glob(os.path.join(a.annot, "*.parquet")))[0] if os.path.isdir(a.annot) else a.annot
    schema = pq.read_schema(first).names
    pcols = [c for c in schema if c.startswith("phylop_")]
    cols = ["subfamily", "compartment", "orientation", "gene_id", "chrom", "max_tpm"] + pcols + \
        [c for c in ("gn_poly", "gn_callable", "cpg") if c in schema]
    df = load_annot(a.annot, columns=cols, query=a.query, subfamilies=a.subfamilies.split(","))
    df = df[df["compartment"].isin(["exon", "intron"])]
    for c in ("subfamily", "orientation", "gene_id"):
        df[c] = df[c].astype(str)
    if a.pool_orientations:
        df = pd.concat([df, df.assign(orientation="pooled")], ignore_index=True)
    g = per_gene(df, pcols, a.min_bases)
    if g.empty:
        raise SystemExit("no gene has enough exonic and intronic bases")
    gtpm = df.groupby("gene_id")["max_tpm"].first()
    g["gtex"] = np.log2(g["gene_id"].map(gtpm).astype(float) + 0.1)

    gt = pd.read_csv(a.gene_table, sep="\t")
    gt["gene_id"] = gt["gene_id"].astype(str).str.split(".").str[0]
    g["_gid"] = g["gene_id"].str.split(".").str[0]
    props = [p for p in a.props.split(",") if p in gt.columns]
    g = g.merge(gt[["gene_id"] + props].rename(columns={"gene_id": "_gid"}), on="_gid", how="left")
    props = props + ["gtex"]
    layers = [c for c in g.columns if c.startswith("diff_")]

    if a.baseline:
        b = pd.read_csv(a.baseline, sep="\t")
        b = b[b["orientation"].isin(["sense", "pooled"])].drop_duplicates("gene_id")
        b = b.set_index("gene_id")
        for c in layers:
            if c in b.columns:
                g[c.replace("diff_", "excess_")] = g[c] - g["gene_id"].map(b[c])
        layers += [c.replace("diff_", "excess_") for c in layers if c in b.columns]
    g.drop(columns="_gid").to_csv(f"{a.out}_genes.tsv", sep="\t", index=False)

    rows = []
    for (sf, ori), s in g.groupby(["subfamily", "orientation"]):
        for layer in layers:
            for p in props:
                d = s[[layer, p]].dropna()
                if len(d) < 15:
                    continue
                rho, pv = spearmanr(d[layer], d[p])
                x, y = d[layer].to_numpy(), d[p].to_numpy()
                boots = []
                for _ in range(a.n_boot):
                    k = rng.integers(0, len(d), len(d))
                    r_ = spearmanr(x[k], y[k])[0]
                    if np.isfinite(r_):
                        boots.append(r_)
                lo, hi = (np.percentile(boots, [2.5, 97.5]) if boots else (np.nan, np.nan))
                rows.append({"subfamily": sf, "orientation": ori, "layer": layer, "property": p, "genes": len(d),
                             "median_diff": float(np.median(x)), "spearman_rho": rho, "ci_lo": lo, "ci_hi": hi,
                             "p": pv})
    t = pd.DataFrame(rows)
    t.to_csv(f"{a.out}_trend.tsv", sep="\t", index=False)
    pd.set_option("display.width", 220)
    log(t.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
