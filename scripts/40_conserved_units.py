#!/usr/bin/env python3
"""Revision step 40 (patch 5b). Step-15 unit table (best-stretch phyloP per unit) for annotations that lack the
'subfamily_orig' column, such as the repeat-free windows of step 9.

Uses unit_table() of 15_conserved_fragments.py unchanged; a missing 'subfamily_orig' is filled with 'subfamily'.
Output: <out>_units.tsv.gz, the same layout as step 15 (units shorter than --min-bases dropped, chrX excluded).
"""
import argparse
import glob
import importlib.util
import os

import pandas as pd
import pyarrow.parquet as pq

from ridl_utils import log


def load_15():
    here = os.path.dirname(os.path.abspath(__file__))
    spec = importlib.util.spec_from_file_location("cons15", os.path.join(here, "15_conserved_fragments.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annot", required=True)
    ap.add_argument("--exclude-chroms", default="chrX")
    ap.add_argument("--window", type=int, default=20)
    ap.add_argument("--min-bases", type=int, default=30)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    m15 = load_15()
    excl = set(c for c in a.exclude_chroms.split(",") if c)
    files = sorted(glob.glob(os.path.join(a.annot, "*.parquet")))
    if not files:
        raise SystemExit(f"no parquet files in {a.annot}")
    schema = pq.read_schema(files[0]).names
    pcols = [c for c in schema if c.startswith("phylop_")]
    want = ["unit_id", "subfamily", "subfamily_orig", "compartment", "orientation", "gene_id", "chrom", "pos"] + pcols
    cols = [c for c in want if c in schema]
    missing = set(want) - set(cols)
    if missing - {"subfamily_orig"}:
        raise SystemExit(f"annotation lacks columns {sorted(missing - {'subfamily_orig'})}")
    parts = []
    for f in files:
        df = pd.read_parquet(f, columns=cols)
        if df.empty or df["chrom"].iloc[0] in excl:
            continue
        if "subfamily_orig" not in df:
            df["subfamily_orig"] = df["subfamily"]
        parts.append(m15.unit_table(df, pcols, a.window))
        log(f"  {os.path.basename(f)}: {len(parts[-1])} units")
    u = pd.concat(parts, ignore_index=True)
    u = u[u["n_bases"] >= a.min_bases]
    u.to_csv(f"{a.out}_units.tsv.gz", sep="\t", index=False, float_format="%.4g")
    log(f"{len(u)} units: " + ", ".join(f"{k} {v}" for k, v in u["compartment"].value_counts().items()))


if __name__ == "__main__":
    main()
