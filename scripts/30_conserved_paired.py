#!/usr/bin/env python3
"""Revision step 30. Conserved stretches of single TE copies, exon vs intron, with the paper's own controls:
the same genes and matched GC content (the analysis of Fig. 4a pooled all genes and did not match GC).

Input: <out>_units.tsv.gz of 15_conserved_fragments.py (one row per unit with its best --window-nt stretch of
each phyloP track) and the genome FASTA (GC of each unit). For each family (--level) a unit is 'conserved' if
its best stretch exceeds --threshold. Exonic and intronic units are compared by a Mantel-Haenszel odds ratio
stratified by unit-length quartile x GC tertile, (a) across all genes and (b) only within genes that contain
both exonic and intronic units of that family. 95% intervals from --n-boot resamples of genes. Families
without any conserved exonic unit are reported but not tested. Use --track prim for Alu.

Output: <out>_paired.tsv
"""
import argparse

import numpy as np
import pandas as pd
import pysam

from ridl_utils import log


def mh_or(conserved, exon, strata, w=None):
    w = np.ones(len(conserved)) if w is None else w
    s = pd.factorize(strata)[0]
    ns = s.max() + 1 if len(s) else 0
    bc = lambda m: np.bincount(s, weights=w * m, minlength=ns)
    A, B, C, D = bc(conserved & exon), bc(~conserved & exon), bc(conserved & ~exon), bc(~conserved & ~exon)
    n = A + B + C + D
    with np.errstate(divide="ignore", invalid="ignore"):
        R = np.nansum(np.where(n > 0, A * D / n, 0))
        S = np.nansum(np.where(n > 0, B * C / n, 0))
    return R / S if S > 0 else np.nan


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--units", required=True)
    ap.add_argument("--genome-fasta", required=True)
    ap.add_argument("--track", default="p447")
    ap.add_argument("--threshold", type=float, default=1.5)
    ap.add_argument("--level", default="subfamily_orig", choices=["subfamily_orig", "subfamily"])
    ap.add_argument("--families", default=None, help="regex on family names (default: all)")
    ap.add_argument("--exclude-chroms", default="chrX")
    ap.add_argument("--n-boot", type=int, default=500)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    u = pd.read_csv(a.units, sep="\t", low_memory=False)
    u = u[u["compartment"].isin(["exon", "intron"]) & ~u["chrom"].isin(a.exclude_chroms.split(","))].copy()
    if a.families:
        u = u[u[a.level].astype(str).str.contains(a.families, regex=True)]
    col = f"best_phylop_{a.track}"
    u = u.dropna(subset=[col])
    fa = pysam.FastaFile(a.genome_fasta)
    gc = []
    for c, s, e in zip(u["chrom"], u["start"], u["end"]):
        q = fa.fetch(c, int(s), int(e)).upper()
        gc.append((q.count("G") + q.count("C")) / max(len(q), 1))
    u["gc"] = gc
    rows = []
    for fam, g in u.groupby(a.level):
        g = g.copy()
        if g["compartment"].nunique() < 2 or len(g) < 40:
            continue
        g["lq"] = pd.qcut(g["n_bases"].rank(method="first"), 4, labels=False)
        g["gq"] = pd.qcut(g["gc"].rank(method="first"), 3, labels=False)
        g["_st"] = g["lq"].astype(str) + "|" + g["gq"].astype(str)
        both = set(g.loc[g["compartment"] == "exon", "gene_id"]) & set(g.loc[g["compartment"] == "intron", "gene_id"])
        for design, s in (("all genes", g), ("same genes", g[g["gene_id"].isin(both)])):
            if s["compartment"].nunique() < 2:
                continue
            cons = (s[col] > a.threshold).to_numpy()
            ex = (s["compartment"] == "exon").to_numpy()
            st = s["_st"].to_numpy()
            r = {"family": fam, "track": a.track, "threshold": a.threshold, "design": design,
                 "units_exon": int(ex.sum()), "units_intron": int((~ex).sum()),
                 "share_exon": cons[ex].mean() if ex.any() else np.nan,
                 "share_intron": cons[~ex].mean() if (~ex).any() else np.nan,
                 "conserved_exon": int((cons & ex).sum()), "conserved_intron": int((cons & ~ex).sum()),
                 "genes": s["gene_id"].nunique()}
            if r["conserved_exon"] == 0 and r["conserved_intron"] == 0:
                r.update(mh_or=np.nan, ci_lo=np.nan, ci_hi=np.nan, note="no conserved units")
            elif r["conserved_exon"] == 0:
                r.update(mh_or=0.0, ci_lo=np.nan, ci_hi=np.nan, note="no conserved exonic unit (not testable)")
            elif r["conserved_intron"] == 0:
                r.update(mh_or=np.inf, ci_lo=np.nan, ci_hi=np.nan, note="no conserved intronic unit (not testable)")
            else:
                obs = mh_or(cons, ex, st)
                genes, inv = np.unique(s["gene_id"].astype(str).to_numpy(), return_inverse=True)
                bs = []
                for _ in range(a.n_boot):
                    cnt = np.bincount(rng.integers(0, len(genes), len(genes)), minlength=len(genes)).astype(float)
                    bs.append(mh_or(cons, ex, st, cnt[inv]))
                bs = np.array(bs, dtype=float)
                r.update(mh_or=obs, ci_lo=np.nanpercentile(bs, 2.5), ci_hi=np.nanpercentile(bs, 97.5), note="")
            rows.append(r)
    out = pd.DataFrame(rows)
    out.to_csv(f"{a.out}_paired.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    log(out.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
