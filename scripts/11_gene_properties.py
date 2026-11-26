#!/usr/bin/env python3
"""Per-gene RNA properties for dose analyses, from the RBP-project export (GENCODE v29, control samples).

Each property is the mean of the HepG2 and K562 values (or the one available):
  stab    stability_6h_log2      relative stability (BruChase 6 h / Bru)
  synth   synthesis_log2rpkm     synthesis (Bru-seq)
  nuc     nuc_cyto_log2          log2 nucleus / cytoplasm (poly(A) fractionation)
  splice  splicing_eff           -log2(intron / exon signal density); only where splicing_intron_reads >= 10
  expr    log2(tpm_polyA + 0.1)  expression in control samples
  tau     tau_gtex               tissue specificity
Tertiles (<prop>_t = 1, 2, 3) are computed over lncRNA genes (project classes) with a value.
Output: TSV with gene_id (no version) and these columns; pass it to step 6 as --gene-table.
"""
import argparse

import numpy as np
import pandas as pd

from ridl_utils import log

PROPS = {"stab": "stability_6h_log2", "synth": "synthesis_log2rpkm", "nuc": "nuc_cyto_log2",
         "splice": "splicing_eff", "expr": "tpm_polyA"}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--export", required=True, help="gene_properties_controls_v29.tsv.gz")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-intron-reads", type=float, default=10)
    a = ap.parse_args()
    d = pd.read_csv(a.export, sep="\t", low_memory=False)
    out = pd.DataFrame({"gene_id": d["gene_id"].astype(str).str.split(".").str[0]})
    lnc = d["is_lncRNA_project"].astype(str).str.lower().eq("true").to_numpy()
    for name, base in PROPS.items():
        vals = []
        for cell in ("HepG2", "K562"):
            c = f"{base}_{cell}"
            if c not in d.columns:
                continue
            v = d[c].astype(float)
            if name == "expr":
                v = np.log2(v + 0.1)
            if name == "splice":
                r = d.get(f"splicing_intron_reads_{cell}")
                if r is not None:
                    v = v.where(r.astype(float) >= a.min_intron_reads)
            vals.append(v)
        out[name] = pd.concat(vals, axis=1).mean(axis=1, skipna=True) if vals else np.nan
    out["tau"] = d["tau_gtex"].astype(float) if "tau_gtex" in d.columns else np.nan
    for name in list(PROPS) + ["tau"]:
        ref = out.loc[lnc & out[name].notna(), name]
        if len(ref) < 30:
            out[f"{name}_t"] = np.nan
            continue
        q1, q2 = ref.quantile([1 / 3, 2 / 3]).to_numpy()
        out[f"{name}_t"] = np.select([out[name] <= q1, out[name] <= q2, out[name] > q2], [1, 2, 3], np.nan)
        out.loc[out[name].isna(), f"{name}_t"] = np.nan
        log(f"  {name}: lncRNA with value {len(ref)}, tertile cut-offs {q1:.3f} / {q2:.3f}")
    out = out.drop_duplicates("gene_id")
    out.to_csv(a.out, sep="\t", index=False)
    log(f"wrote {len(out)} genes -> {a.out}")


if __name__ == "__main__":
    main()
