#!/usr/bin/env python3
"""Compare the result tables of two runs (e.g. the original work/ and the frozen rerun).

Every .tsv file (also .tsv.gz) under the new work directory is matched by relative path in the old one and
compared column by column: numeric columns by the largest absolute and relative difference, other columns
by exact equality. Files that are missing in either run are listed. Large per-base tables are skipped.

Usage: compare_runs.py OLD_WORK NEW_WORK [--tol 1e-9] [--out report.tsv]
Tables found only in the old run (e.g. outputs of superseded runs) are counted but not listed.
Exit status 0 if every shared table agrees within --tol, 1 otherwise.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

SKIP_DIRS = {"synthetic", "synthetic_shift"}
MAX_BYTES = 200_000_000


def tables(root):
    out = {}
    for d, dirs, files in os.walk(root):
        dirs[:] = [x for x in dirs if x not in SKIP_DIRS]
        for f in files:
            if f.endswith((".tsv", ".tsv.gz")):
                p = os.path.join(d, f)
                if os.path.getsize(p) <= MAX_BYTES:
                    out[os.path.relpath(p, root)] = p
    return out


def compare(a, b, tol):
    try:
        x = pd.read_csv(a, sep="\t", low_memory=False)
        y = pd.read_csv(b, sep="\t", low_memory=False)
    except pd.errors.EmptyDataError:
        return ("same", "both empty", 0.0) if os.path.getsize(a) == os.path.getsize(b) else ("DIFF", "one empty", np.nan)
    if list(x.columns) != list(y.columns):
        return "DIFF", f"columns differ: {sorted(set(x.columns) ^ set(y.columns))[:6]}", np.nan
    if len(x) != len(y):
        return "DIFF", f"rows {len(x)} vs {len(y)}", np.nan
    worst, where = 0.0, ""
    for c in x.columns:
        u, v = x[c], y[c]
        if pd.api.types.is_numeric_dtype(u) and pd.api.types.is_numeric_dtype(v):
            u, v = u.to_numpy(float), v.to_numpy(float)
            nan_mismatch = np.isnan(u) != np.isnan(v)
            if nan_mismatch.any():
                return "DIFF", f"{c}: missing values differ in {int(nan_mismatch.sum())} rows", np.nan
            ok = ~np.isnan(u)
            if ok.any():
                d = np.abs(u[ok] - v[ok]) / np.maximum(1.0, np.abs(u[ok]))
                if d.max() > worst:
                    worst, where = float(d.max()), c
        elif not u.astype(str).equals(v.astype(str)):
            return "DIFF", f"{c}: values differ in {int((u.astype(str) != v.astype(str)).sum())} rows", np.nan
    return ("same" if worst <= tol else "DIFF"), (f"largest difference in {where}" if where else ""), worst


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("old")
    ap.add_argument("new")
    ap.add_argument("--tol", type=float, default=1e-9, help="largest allowed relative difference (min scale 1)")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    old, new = tables(a.old), tables(a.new)
    rows = []
    for rel in sorted(set(old) | set(new)):
        if rel not in old:
            rows.append((rel, "only new", "", np.nan))
        elif rel not in new:
            rows.append((rel, "only old", "", np.nan))
        else:
            rows.append((rel, *compare(old[rel], new[rel], a.tol)))
    rep = pd.DataFrame(rows, columns=["table", "status", "note", "max_rel_diff"])
    pd.set_option("display.width", 200)
    pd.set_option("display.max_colwidth", 90)
    print(rep.status.value_counts().to_string())
    bad = rep[rep.status.isin(["DIFF", "only new"])]
    if len(bad):
        print(bad.to_string(index=False))
    if a.out:
        rep.to_csv(a.out, sep="\t", index=False)
    sys.exit(0 if (rep.status == "DIFF").sum() == 0 else 1)


if __name__ == "__main__":
    main()
