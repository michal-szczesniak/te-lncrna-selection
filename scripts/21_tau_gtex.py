#!/usr/bin/env python3
"""Tissue specificity (tau) from GTEx median TPM, replacing the tau column of the gene-property table.

tau = sum_i (1 - x_i / max(x)) / (n - 1), with x = log2(median TPM + 1) over tissues (Yanai et al. 2005);
cultured cells (tissue names starting with "Cells") are excluded by default. Genes with zero expression in all
tissues get no value. The previous column is kept as tau_export, and tertiles (tau_t = 1, 2, 3) are recomputed
over genes with a value, as in 11_gene_properties.py.

Output: the gene-property table with the new tau and tau_t; prints the agreement with the previous values.
"""
import argparse
import re

import numpy as np
import pandas as pd

from ridl_utils import log


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gct", required=True, help="GTEx gene median TPM (.gct.gz)")
    ap.add_argument("--gene-table", required=True, help="gene_props.tsv from 11_gene_properties.py")
    ap.add_argument("--exclude-regex", default=r"^Cells")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    g = pd.read_csv(a.gct, sep="\t", skiprows=2)
    tissues = [c for c in g.columns[2:] if not re.search(a.exclude_regex, c)]
    log(f"{len(g):,} genes, {len(tissues)} tissues used ({len(g.columns) - 2 - len(tissues)} excluded)")
    x = np.log2(g[tissues].astype(float) + 1)
    mx = x.max(axis=1)
    tau = (1 - x.div(mx.where(mx > 0), axis=0)).sum(axis=1, min_count=1) / (len(tissues) - 1)
    tau = pd.Series(tau.values, index=g["Name"].str.replace(r"\.\d+$", "", regex=True))
    tau = tau[~tau.index.duplicated()]

    gp = pd.read_csv(a.gene_table, sep="\t")
    if "tau" in gp.columns:
        gp = gp.rename(columns={"tau": "tau_export"})
    gp["tau"] = gp["gene_id"].map(tau)
    ok = gp["tau"].notna()
    gp["tau_t"] = np.nan
    gp.loc[ok, "tau_t"] = pd.qcut(gp.loc[ok, "tau"].rank(method="first"), 3, labels=[1, 2, 3]).astype(float)
    log(f"tau assigned to {ok.sum():,} of {len(gp):,} genes")
    if "tau_export" in gp.columns:
        both = gp[["tau", "tau_export"]].dropna()
        if len(both) > 2:
            rho = both["tau"].corr(both["tau_export"], method="spearman")
            log(f"agreement with previous tau: Spearman rho = {rho:.3f} over {len(both):,} genes; "
                f"median {both['tau'].median():.3f} vs {both['tau_export'].median():.3f}")
    gp.to_csv(a.out, sep="\t", index=False)
    log(f"wrote {a.out}")


if __name__ == "__main__":
    main()
