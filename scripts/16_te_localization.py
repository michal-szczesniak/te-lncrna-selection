#!/usr/bin/env python3
"""Replication: do exonic TE fragments go with nuclear enrichment of lncRNAs (Carlevaro-Fita et al. 2019)?

Per gene (reference transcript from step 1: full exons, introns between them) the script measures how many
bases of each TE category lie in exons and in introns (RepeatMasker .out), and relates exonic and intronic
TE density to the nucleus/cytoplasm ratio (--gene-table column 'nuc', e.g. from 11_gene_properties.py).
An RNA-level effect predicts an association for EXONIC density only; intronic density is the control,
because introns are not part of the mature RNA.

Model per category (OLS, robust HC1 standard errors):
  nuc ~ z(exonic density) + z(intronic density) + log10(exon length) + log10(intron length + 1)
        + exonic GC + expression covariate (--expr-col)
Also a dose table (mean nuc by exonic bp of the category) and, with --units (from 15_conserved_fragments),
a model with the presence of a CONSERVED exonic fragment of MIR/L2 (best 20-nt phyloP > --cons-thr).

Outputs: <out>_genes.tsv.gz, <out>_models.tsv, <out>_dose.tsv
"""
import argparse

import numpy as np
import pandas as pd
import pysam
from scipy.stats import norm

from ridl_utils import log, open_text

CATS = {
    "L2b_MIRb_MIRc": lambda n, c, f: n in ("L2b", "MIRb", "MIRc"),
    "MIR_all": lambda n, c, f: f == "MIR",
    "L2_all": lambda n, c, f: f == "L2",
    "Alu": lambda n, c, f: f == "Alu",
    "L1": lambda n, c, f: f == "L1",
    "LTR": lambda n, c, f: c == "LTR",
    "DNA": lambda n, c, f: c == "DNA",
}


def ols_hc1(y, X, names):
    n, k = X.shape
    XtX_inv = np.linalg.pinv(X.T @ X)
    b = XtX_inv @ X.T @ y
    e = y - X @ b
    meat = (X * e[:, None] ** 2).T @ X
    V = XtX_inv @ meat @ XtX_inv * n / max(n - k, 1)
    se = np.sqrt(np.diag(V))
    p = 2 * norm.sf(np.abs(b / se))
    return pd.DataFrame({"term": names, "coef": b, "se": se, "p": p})


def overlap_bp(iv, reps):
    """iv: DataFrame gene_id, chrom, start, end, strand; reps: per chrom (starts, ends, cats, strands).
    Returns dict (gene_id, category, orientation) -> bp."""
    out = {}
    for r in iv.itertuples(index=False):
        if r.chrom not in reps:
            continue
        rs, re_, rc, rst = reps[r.chrom]
        # repeats are sorted by start; their ends are not monotonic -> start 50 kb upstream (longer than any
        # repeat record in these categories) and scan forward
        i = int(np.searchsorted(rs, r.start - 50000))
        while i < len(rs) and rs[i] < r.end:
            if re_[i] > r.start:
                bp = min(re_[i], r.end) - max(rs[i], r.start)
                ori = "sense" if rst[i] == r.strand else "antisense"
                for c in rc[i]:
                    for key in ((r.gene_id, c, ori), (r.gene_id, c, "all")):
                        out[key] = out.get(key, 0) + bp
            i += 1
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref-tx", required=True, help="step-1 ref_tx.tsv.gz")
    ap.add_argument("--rmsk-out", required=True)
    ap.add_argument("--genome-fasta", required=True)
    ap.add_argument("--gene-table", required=True, help="TSV with gene_id (no version) and 'nuc'")
    ap.add_argument("--expr-col", default="synth")
    ap.add_argument("--nuc-col", default="nuc", help="column of --gene-table with log2 nucleus/cytoplasm")
    ap.add_argument("--min-exons", type=int, default=1, help="e.g. 2 = multi-exonic genes only")
    ap.add_argument("--units", default=None, help="<out>_units.tsv.gz from 15_conserved_fragments.py")
    ap.add_argument("--cons-thr", type=float, default=1.5)
    ap.add_argument("--chroms", default=",".join(f"chr{i}" for i in range(1, 23)))
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    chroms = set(a.chroms.split(","))

    tx = pd.read_csv(a.ref_tx, sep="\t", dtype={"chrom": str})
    tx = tx[tx["chrom"].isin(chroms)].sort_values(["gene_id", "tx_offset"])
    ex = tx[["gene_id", "chrom", "start", "end", "strand"]].copy()
    introns = []
    for gid, g in tx.groupby("gene_id"):
        s = g.sort_values("start")
        st, en = s["start"].to_numpy(), s["end"].to_numpy()
        for k in range(len(s) - 1):
            if st[k + 1] > en[k]:
                introns.append((gid, s["chrom"].iloc[0], en[k], st[k + 1], s["strand"].iloc[0]))
    intr = pd.DataFrame(introns, columns=["gene_id", "chrom", "start", "end", "strand"])
    log(f"genes {ex['gene_id'].nunique()}, exons {len(ex)}, introns {len(intr)}")

    rows = []
    with open_text(a.rmsk_out) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 15 or not f[0].isdigit() or f[4] not in chroms:
                continue
            cls_fam = f[10].split("/")
            cls, fam = cls_fam[0], (cls_fam[1] if len(cls_fam) > 1 else "")
            cats = tuple(c for c, fn in CATS.items() if fn(f[9], cls, fam))
            if cats:
                rows.append((f[4], int(f[5]) - 1, int(f[6]), "+" if f[8] == "+" else "-", cats))
    rp = pd.DataFrame(rows, columns=["chrom", "start", "end", "strand", "cats"]).sort_values(["chrom", "start"])
    reps = {c: (g["start"].to_numpy(), g["end"].to_numpy(), list(g["cats"]), g["strand"].to_numpy())
            for c, g in rp.groupby("chrom")}
    log(f"repeat records in categories: {len(rp)}")

    ov_ex, ov_in = overlap_bp(ex, reps), overlap_bp(intr, reps)
    fa = pysam.FastaFile(a.genome_fasta)
    gene = ex.groupby("gene_id").agg(chrom=("chrom", "first"), strand=("strand", "first"),
                                     n_exons=("start", "size")).reset_index()
    ex["len"] = ex["end"] - ex["start"]
    gene["exon_len"] = gene["gene_id"].map(ex.groupby("gene_id")["len"].sum())
    intr["len"] = intr["end"] - intr["start"]
    gene["intron_len"] = gene["gene_id"].map(intr.groupby("gene_id")["len"].sum()).fillna(0)
    gc = {}
    for r in ex.itertuples(index=False):
        sq = fa.fetch(r.chrom, int(r.start), int(r.end)).upper()
        g0 = gc.get(r.gene_id, [0, 0])
        gc[r.gene_id] = [g0[0] + sq.count("G") + sq.count("C"), g0[1] + len(sq)]
    gene["exon_gc"] = gene["gene_id"].map({k: v[0] / max(v[1], 1) for k, v in gc.items()})
    for c in CATS:
        for ori in ("all", "sense", "antisense"):
            gene[f"ex_{c}_{ori}"] = [ov_ex.get((g, c, ori), 0) for g in gene["gene_id"]]
            gene[f"in_{c}_{ori}"] = [ov_in.get((g, c, ori), 0) for g in gene["gene_id"]]

    if a.units:
        u = pd.read_csv(a.units, sep="\t", usecols=["subfamily", "compartment", "gene_id", "best_phylop_p447"])
        cons = u[(u["compartment"] == "exon") & u["subfamily"].isin(["MIR*", "L2abc", "MIR3", "L2d*"])
                 & (u["best_phylop_p447"] > a.cons_thr)]
        gene["ex_conserved_MIR_L2"] = gene["gene_id"].isin(set(cons["gene_id"])).astype(int)

    gt = pd.read_csv(a.gene_table, sep="\t")
    gt["gene_id"] = gt["gene_id"].astype(str).str.split(".").str[0]
    gene["_gid"] = gene["gene_id"].str.split(".").str[0]
    gt = gt.rename(columns={a.nuc_col: "nuc"})
    keep = ["nuc", a.expr_col] if a.expr_col in gt.columns else ["nuc"]
    gene = gene.merge(gt[["gene_id"] + keep].rename(columns={"gene_id": "_gid"}).drop_duplicates("_gid"),
                      on="_gid", how="left").drop(columns="_gid")
    gene.to_csv(f"{a.out}_genes.tsv.gz", sep="\t", index=False, float_format="%.5g")

    d = gene.dropna(subset=["nuc", "exon_gc"] + ([a.expr_col] if a.expr_col in gene.columns else [])).copy()
    d = d[d["n_exons"] >= a.min_exons]
    log(f"genes with nucleus/cytoplasm values and covariates: {len(d)}")
    d["log_exlen"] = np.log10(d["exon_len"])
    d["log_inlen"] = np.log10(d["intron_len"] + 1)
    base = ["log_exlen", "log_inlen", "exon_gc"] + ([a.expr_col] if a.expr_col in d.columns else [])
    z = lambda x: (x - x.mean()) / (x.std() if x.std() > 0 else 1)
    res = []
    for c in CATS:
        for ori in ("all", "sense", "antisense"):
            dex = d[f"ex_{c}_{ori}"] / d["exon_len"]
            din = d[f"in_{c}_{ori}"] / d["intron_len"].replace(0, np.nan)
            din = din.fillna(0)
            X = np.column_stack([np.ones(len(d)), z(dex), z(din)] + [z(d[b]) for b in base])
            m = ols_hc1(d["nuc"].to_numpy(float), X, ["intercept", "exonic_density", "intronic_density"] + base)
            for term in ("exonic_density", "intronic_density"):
                r = m[m.term == term].iloc[0]
                res.append({"category": c, "orientation": ori, "term": term, "coef_per_sd": r.coef, "se": r.se,
                            "p": r.p, "genes": len(d), "genes_with_exonic": int((d[f"ex_{c}_{ori}"] > 0).sum())})
    if "ex_conserved_MIR_L2" in d.columns:
        X = np.column_stack([np.ones(len(d)), d["ex_conserved_MIR_L2"]] + [z(d[b]) for b in base])
        m = ols_hc1(d["nuc"].to_numpy(float), X, ["intercept", "conserved_exonic_MIR_L2"] + base)
        r = m[m.term == "conserved_exonic_MIR_L2"].iloc[0]
        res.append({"category": f"conserved MIR/L2 (best 20-nt phyloP > {a.cons_thr})", "orientation": "all",
                    "term": "exonic_presence", "coef_per_sd": r.coef, "se": r.se, "p": r.p, "genes": len(d),
                    "genes_with_exonic": int(d["ex_conserved_MIR_L2"].sum())})
    res = pd.DataFrame(res)
    res.to_csv(f"{a.out}_models.tsv", sep="\t", index=False)

    bins = [-1, 0, 100, 300, 10 ** 9]
    labels = ["0", "1-100", "101-300", ">300"]
    dose = []
    for c in ("L2b_MIRb_MIRc", "Alu"):
        for part in ("ex", "in"):
            b = pd.cut(d[f"{part}_{c}_all"], bins, labels=labels)
            t = d.groupby(b, observed=True)["nuc"].agg(["size", "mean", "median"]).reset_index()
            t.insert(0, "where", "exon" if part == "ex" else "intron")
            t.insert(0, "category", c)
            dose.append(t.rename(columns={f"{part}_{c}_all": "bp"}))
    dose = pd.concat(dose, ignore_index=True)
    dose.columns = ["category", "where", "bp", "genes", "mean_nuc", "median_nuc"]
    dose.to_csv(f"{a.out}_dose.tsv", sep="\t", index=False)
    pd.set_option("display.width", 200)
    log(res[res.orientation == "all"].round(4).to_string(index=False))
    log(dose.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
