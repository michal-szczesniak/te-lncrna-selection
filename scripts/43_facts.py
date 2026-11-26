#!/usr/bin/env python3
"""Revision step 43 (patch 6). Counts quoted in the manuscript that no result table holds, recomputed from the
pipeline outputs, so that every number in the text has a script behind it.

  copies        Alu, MIR and L2 copies placed on consensus coordinates (step-2 te_copies.tsv.gz, autosomes, all
                compartments; units of a copy that spans an exon boundary count once) and their subfamilies
  mapped        the same after mapping and grouping (step-3 copies_grouped.tsv.gz, if given)
  genes         kept lncRNA genes (step 1) and genes with exonic sequence in the compartments
  exonic        exonic sequence left after filtering, in bp and as a share of the exons of the kept genes
                (step-20 attrition table, if given)

Output: <out>.tsv (quantity, value, source).
"""
import argparse
import os

import pandas as pd

from ridl_utils import log

AUTO = [f"chr{i}" for i in range(1, 23)]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--copies", required=True, help="step-2 te_copies.tsv.gz of the main set")
    ap.add_argument("--copies-grouped", default=None, help="step-3 copies_grouped.tsv.gz")
    ap.add_argument("--step1", required=True, help="step-1 directory (lnc_genes.tsv.gz, compartments.tsv.gz)")
    ap.add_argument("--attrition", default=None, help="step-20 attrition.tsv")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rows = []

    def add(q, v, src):
        rows.append({"quantity": q, "value": v, "source": src})

    c = pd.read_csv(a.copies, sep="\t", low_memory=False)
    c = c[c["chrom"].isin(AUTO)]
    fam = c["subfamily"].astype(str).str.extract(r"^(Alu|MIR|L2)", expand=False)
    c = c[fam.notna()]
    add("Alu/MIR/L2 copies, autosomes, all compartments", int(c["copy_id"].nunique()), os.path.basename(a.copies))
    add("Alu/MIR/L2 subfamilies", int(c["subfamily"].nunique()), os.path.basename(a.copies))
    for comp, d in c.groupby("compartment"):
        add(f"copies with a {comp} unit", int(d["copy_id"].nunique()), os.path.basename(a.copies))
    if a.copies_grouped and os.path.exists(a.copies_grouped):
        g = pd.read_csv(a.copies_grouped, sep="\t", low_memory=False)
        g = g[g["chrom"].isin(AUTO)]
        add("copies after mapping and grouping, autosomes", int(g["copy_id"].nunique()),
            os.path.basename(a.copies_grouped))
        if "subfamily_orig" in g:
            add("subfamilies after grouping", int(g["subfamily_orig"].nunique()), os.path.basename(a.copies_grouped))
    lnc = pd.read_csv(os.path.join(a.step1, "lnc_genes.tsv.gz"), sep="\t", low_memory=False)
    lnc_a = lnc[lnc["chrom"].isin(AUTO)]
    add("autosomal lncRNA genes in the annotation", len(lnc_a), "lnc_genes.tsv.gz")
    add("kept lncRNA genes (autosomes)", int((lnc_a["drop_reason"].fillna("") == "").sum()), "lnc_genes.tsv.gz")
    for r, n in lnc_a["drop_reason"].fillna("").value_counts().items():
        if r:
            add(f"genes removed: {r}", int(n), "lnc_genes.tsv.gz")
    comp = pd.read_csv(os.path.join(a.step1, "compartments.tsv.gz"), sep="\t", low_memory=False)
    ex = comp[(comp["compartment"] == "exon") & comp["chrom"].isin(AUTO)]
    add("genes with exonic sequence after filtering (autosomes)", int(ex["gene_id"].nunique()), "compartments.tsv.gz")
    add("exonic sequence after filtering, bp (autosomes)", int((ex["end"] - ex["start"]).sum()), "compartments.tsv.gz")
    if a.attrition and os.path.exists(a.attrition):
        t = pd.read_csv(a.attrition, sep="\t")
        t = t[~t.iloc[:, 0].astype(str).str.startswith("ALTERNATIVE")]  # rows describing an alternative filter
        log("attrition table:\n" + t.to_string(index=False))
        for col in t.columns:
            if "share" in col.lower() or "left" in col.lower():
                add(f"attrition, last step, {col}", t[col].iloc[-1], os.path.basename(a.attrition))
    out = pd.DataFrame(rows)
    out.to_csv(f"{a.out}.tsv", sep="\t", index=False)
    log("\n" + out.to_string(index=False))


if __name__ == "__main__":
    main()
