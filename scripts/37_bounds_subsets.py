#!/usr/bin/env python3
"""Revision step 37 (patch 5). Limits on TE-specific constraint in subsets of lncRNAs more likely to be functional.

If RNA-level TE function were concentrated in the minority of functional lncRNAs, averaging over all genes would
dilute it. Without an independent list of functional lncRNAs, expression is used as the proxy: highly expressed
and broadly expressed lncRNAs are more often conserved and functional. The per-gene table of step 12
(te_vs_nonte_genes.tsv: exon-intron differences of TE copies and their excess over repeat-free sequence of the same
gene) is split into gene subsets and 23_bounds.py is rerun on each:
  expr_top25 / expr_top10   genes in the top quarter / tenth of GTEx v10 maximum expression (column 'gtex',
                            log2 max TPM), quantiles over the genes of the table
  broad_tau25               genes in the lowest quarter of tissue specificity (column 'tau'), if present
  all                       all genes (as in the paper)
With --well-supported (patch 6), lncRNAs with independent support are added:
  gtex_any                  genes with a GTEx v10 expression estimate (genes absent from GTEx removed)
  lincrna_v29               genes already classified as lincRNAs in GENCODE release 29 (--gene-class,
                            column project_class == 'lincRNA')
  lincrna_v29_gtex          both
Thresholds are computed once over the distinct genes of the table, so the same genes enter every TE group.

Outputs: <out>_<subset>_genes.tsv (filtered tables), the 23_bounds.py outputs <out>_bounds.tsv,
<out>_orientation.tsv, <out>_functional.tsv (table column = subset), and <out>_subsets.tsv (thresholds, genes).
"""
import argparse
import os
import subprocess
import sys

import numpy as np
import pandas as pd

from ridl_utils import log


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes", required=True, help="12_gene_trend.py --baseline <out>_genes.tsv")
    ap.add_argument("--comp-dir", required=True, help="as --comp-dir of 23_bounds.py (intronic variant rates)")
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--well-supported", action="store_true", help="add the gtex_any / lincrna_v29 subsets")
    ap.add_argument("--gene-class", default=None, help="gene_class.tsv (gene_id, project_class) for lincrna_v29")
    ap.add_argument("--skip-expression-subsets", action="store_true",
                    help="only 'all' and the --well-supported subsets")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    g = pd.read_csv(a.genes, sep="\t", low_memory=False)
    per_gene = g.drop_duplicates("gene_id").set_index("gene_id")
    subsets, info = {"all": g}, [{"subset": "all", "genes": int(g.gene_id.nunique()), "rule": "all genes"}]
    if "gtex" in per_gene and not a.skip_expression_subsets:
        x = per_gene["gtex"].astype(float)
        for name, q in (("expr_top25", 0.75), ("expr_top10", 0.90)):
            thr = float(np.nanquantile(x, q))
            keep = set(x.index[x >= thr])
            subsets[name] = g[g["gene_id"].isin(keep)]
            info.append({"subset": name, "genes": len(keep), "rule": f"log2(max TPM + 0.1) >= {thr:.3f}"})
    else:
        log("WARNING: no 'gtex' column; expression subsets skipped")
    if "tau" in per_gene and not a.skip_expression_subsets:
        t = per_gene["tau"].astype(float)
        thr = float(np.nanquantile(t, 0.25))
        keep = set(t.index[t <= thr])
        subsets["broad_tau25"] = g[g["gene_id"].isin(keep)]
        info.append({"subset": "broad_tau25", "genes": len(keep), "rule": f"tau <= {thr:.3f}"})
    else:
        log("no 'tau' column; breadth subset skipped")
    if a.well_supported:
        strip = lambda v: pd.Series(v).astype(str).str.replace(r"\.\d+(_PAR_Y)?$", "", regex=True)
        gid = strip(g["gene_id"]).to_numpy()
        have = set()
        if "gtex" in per_gene:
            have = set(strip(per_gene.index[per_gene["gtex"].astype(float).notna()]))
            subsets["gtex_any"] = g[np.isin(gid, list(have))]
            info.append({"subset": "gtex_any", "genes": len(have & set(gid)), "rule": "GTEx v10 estimate present"})
        else:
            log("WARNING: no 'gtex' column; gtex_any skipped")
        if a.gene_class:
            gc = pd.read_csv(a.gene_class, sep="\t", low_memory=False)
            linc = set(strip(gc.loc[gc["project_class"].astype(str) == "lincRNA", "gene_id"]))
            subsets["lincrna_v29"] = g[np.isin(gid, list(linc))]
            info.append({"subset": "lincrna_v29", "genes": len(linc & set(gid)), "rule": "GENCODE 29 lincRNA"})
            if have:
                both = linc & have
                subsets["lincrna_v29_gtex"] = g[np.isin(gid, list(both))]
                info.append({"subset": "lincrna_v29_gtex", "genes": len(both & set(gid)),
                             "rule": "GENCODE 29 lincRNA with GTEx estimate"})
        subsets = {k: v for k, v in subsets.items() if v["gene_id"].nunique() >= 15 or k == "all"}
    args = []
    for name, d in subsets.items():
        f = f"{a.out}_{name}_genes.tsv"
        d.to_csv(f, sep="\t", index=False)
        args += ["--genes", f"{name}={f}"]
    pd.DataFrame(info).to_csv(f"{a.out}_subsets.tsv", sep="\t", index=False)
    log(pd.DataFrame(info).to_string(index=False))
    here = os.path.dirname(os.path.abspath(__file__))
    cmd = [sys.executable, os.path.join(here, "23_bounds.py")] + args + \
        ["--comp-dir", a.comp_dir, "--n-boot", str(a.n_boot), "--out", a.out]
    log(" ".join(cmd))
    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
