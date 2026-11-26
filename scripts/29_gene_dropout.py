#!/usr/bin/env python3
"""Revision step 29. Where does the exonic deficit of variants across all genes come from, and which genes
drop out of the same-gene comparison?

For each group and orientation, exonic (fg) and intronic (bg) units are split by whether their gene contributes
both kinds of units ('paired') or only one ('unpaired'). Whole-element Mantel-Haenszel log ORs of a site
being variable (callable sites; strata --pop-strata) with 95% intervals from --n-boot resamples of genes:
  exon all      vs intron all       the all-genes comparison of the paper (Fig. 2b, open symbols)
  exon paired   vs intron paired    the same-gene comparison (Fig. 2b, filled symbols)
  exon all      vs intron paired    is the shift on the intronic side?
  exon paired   vs intron all       ... or on the exonic side?
  exon unpaired vs exon paired      do exons of genes without intronic copies carry fewer variants?
  intron unpaired vs intron paired  do introns of genes without exonic copies carry more variants?
Gene characteristics (paired vs unpaired genes with exonic units): share of single-exon genes, median
maximum GTEx TPM, gene span, units per gene and any numeric columns of --gene-table.

Outputs: <out>_decomposition.tsv, <out>_genes.tsv
"""
import argparse

import numpy as np
import pandas as pd

from ridl_stats import Agg, bootstrap_weights, load_annot, mh_logor, select_comparison
from ridl_utils import log


def codes(d, cols):
    key = d[cols[0]].astype(str)
    for c in cols[1:]:
        key = key + "|" + d[c].astype(str)
    return pd.factorize(key)[0]


def whole(d, strata, n_boot, rng):
    d = d.copy()
    d["bin"] = 0
    d["_s"] = codes(d, strata)
    d["gn_poly"] = d["gn_poly"].astype(float)
    ag = Agg(d, ["gn_poly"], strata_col="_s", n_bins=1)
    ones = np.ones(len(ag.units))

    def est(w):
        S = ag.sums(ag.label, w)
        a, c = S["gn_poly"]
        n1, n0 = S["_n"]
        return mh_logor(a, n1, c, n0, 1, ag.n_strata)[1]
    obs = est(ones)
    b = np.array([est(bootstrap_weights(ag, rng)) for _ in range(n_boot)], dtype=float)
    return obs, np.nanpercentile(b, 2.5), np.nanpercentile(b, 97.5), int((ag.label == 1).sum()), int((ag.label == 0).sum())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annot", required=True)
    ap.add_argument("--subfamilies", required=True, help="comma list; each analysed separately")
    ap.add_argument("--query", default=None)
    ap.add_argument("--pop-strata", default="tri,cpg,subfamily_orig")
    ap.add_argument("--lnc-genes", default=None, help="step-1 lnc_genes.tsv.gz (gene span)")
    ap.add_argument("--gene-table", default=None, help="per-gene properties (gene_id without version)")
    ap.add_argument("--n-boot", type=int, default=200)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    strata = [c for c in a.pop_strata.split(",") if c]
    cols = ["unit_id", "subfamily", "compartment", "orientation", "gene_id", "exon_class", "max_tpm", "gn_callable",
            "gn_poly"] + strata
    span = None
    if a.lnc_genes:
        g = pd.read_csv(a.lnc_genes, sep="\t")
        span = (g["end"] - g["start"]).groupby(g["gene_id"]).first()
    props = None
    if a.gene_table:
        props = pd.read_csv(a.gene_table, sep="\t")
        props["gene_id"] = props["gene_id"].astype(str).str.split(".").str[0]
        props = props.drop_duplicates("gene_id").set_index("gene_id")
        props = props.select_dtypes("number")
    rows, grows = [], []
    for grp in a.subfamilies.split(","):
        df = load_annot(a.annot, columns=cols, query=a.query, subfamilies=[grp])
        df = df[df["compartment"].isin(["exon", "intron"]) & df["gn_callable"].astype(bool)]
        for ori in sorted(df.loc[df["compartment"] == "exon", "orientation"].unique()):
            d = select_comparison(df, grp, ori, "exon", "intron")
            g1 = set(d.loc[d["_lab"] == 1, "gene_id"]) - {"."}
            g0 = set(d.loc[d["_lab"] == 0, "gene_id"]) - {"."}
            paired = g1 & g0
            d["_paired"] = d["gene_id"].isin(paired)
            ex, it = d[d["_lab"] == 1], d[d["_lab"] == 0]
            designs = {
                "exon all vs intron all": (ex, it),
                "exon paired vs intron paired": (ex[ex["_paired"]], it[it["_paired"]]),
                "exon all vs intron paired": (ex, it[it["_paired"]]),
                "exon paired vs intron all": (ex[ex["_paired"]], it),
                "exon unpaired vs exon paired": (ex[~ex["_paired"]], ex[ex["_paired"]].assign(_lab=0)),
                "intron unpaired vs intron paired": (it[~it["_paired"]].assign(_lab=1), it[it["_paired"]]),
            }
            for name, (f, b) in designs.items():
                if f.empty or b.empty:
                    continue
                x = pd.concat([f.assign(_lab=1), b.assign(_lab=0)], ignore_index=True)
                # units shared by both sides of an exon-vs-exon design are disambiguated by the label in Agg
                e, lo, hi, nf, nb_ = whole(x, strata, a.n_boot, rng)
                rows.append({"group": grp, "orientation": ori, "design": name, "units_fg": nf, "units_bg": nb_,
                             "genes_fg": f["gene_id"].nunique(), "genes_bg": b["gene_id"].nunique(),
                             "mh_log_or": e, "ci_lo": lo, "ci_hi": hi})
            log(f"  {grp} {ori}: genes with exonic units {len(g1)}, paired {len(paired)}")
            # characteristics of genes with exonic units, paired vs unpaired
            ge = ex.groupby("gene_id").agg(units=("unit_id", "nunique"), mono=("exon_class", lambda s: (s == "mono").mean() > 0.5),
                                           max_tpm=("max_tpm", "first"))
            ge["paired"] = ge.index.isin(paired)
            if span is not None:
                ge["span"] = ge.index.map(span)
            if props is not None:
                gid = ge.index.astype(str).str.split(".").str[0]
                for c in props.columns:
                    ge[c] = props[c].reindex(gid).to_numpy()
            for pv, s in ge.groupby("paired"):
                r = {"group": grp, "orientation": ori, "paired": bool(pv), "genes": len(s),
                     "share_single_exon": s["mono"].mean(), "median_units": s["units"].median(),
                     "median_max_tpm": s["max_tpm"].median()}
                for c in [c for c in s.columns if c not in ("units", "mono", "max_tpm", "paired")]:
                    r[f"median_{c}"] = s[c].median()
                grows.append(r)
    R = pd.DataFrame(rows)
    G = pd.DataFrame(grows)
    R.to_csv(f"{a.out}_decomposition.tsv", sep="\t", index=False)
    G.to_csv(f"{a.out}_genes.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    log(R.round(4).to_string(index=False))
    log(G.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
