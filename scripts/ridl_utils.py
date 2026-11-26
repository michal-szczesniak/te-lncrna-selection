"""Shared helpers for the RIDL 2.0 pilot scripts.

Conventions used everywhere:
  * genomic intervals are 0-based, half-open (BED-like): [start, end)
  * consensus positions are 1-based (RepeatMasker convention)
  * strand is '+' or '-'
"""
import gzip
import re
import sys

import numpy as np
import pandas as pd

DEFAULT_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def log(msg):
    print(msg, file=sys.stderr, flush=True)


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith((".gz", ".bgz")) else open(path)


def revcomp(seq):
    return "".join(COMP.get(b, "N") for b in reversed(seq.upper()))


def parse_chroms(arg):
    if arg is None:
        return DEFAULT_CHROMS
    return [c.strip() for c in arg.split(",") if c.strip()]


# ----------------------------------------------------------------------------- GTF
def read_gtf(path, chroms):
    """Return (genes, transcripts, exons) DataFrames with 0-based starts."""
    cols = ["chrom", "src", "feature", "start", "end", "score", "strand", "frame", "attr"]
    df = pd.read_csv(path, sep="\t", comment="#", header=None, names=cols,
                     usecols=["chrom", "feature", "start", "end", "strand", "attr"],
                     dtype={"chrom": str, "feature": str, "attr": str})
    df = df[df["chrom"].isin(chroms) & df["feature"].isin(["gene", "transcript", "exon"])].copy()
    df["start"] = df["start"].astype(np.int64) - 1
    df["end"] = df["end"].astype(np.int64)

    def grab(key):
        return df["attr"].str.extract(rf'{key} "([^"]+)"', expand=False)

    df["gene_id"] = grab("gene_id")
    gt = grab("gene_type")
    if gt.isna().all():
        gt = grab("gene_biotype")  # Ensembl naming
    df["gene_type"] = gt
    df["gene_name"] = grab("gene_name")
    df["transcript_id"] = grab("transcript_id")
    df["canonical"] = df["attr"].str.contains('tag "Ensembl_canonical"', regex=False)
    df = df.drop(columns="attr")
    genes = df[df["feature"] == "gene"].drop(columns=["feature", "transcript_id", "canonical"])
    tx = df[df["feature"] == "transcript"].drop(columns="feature")
    ex = df[df["feature"] == "exon"].drop(columns=["feature", "canonical"])
    return genes.reset_index(drop=True), tx.reset_index(drop=True), ex.reset_index(drop=True)


# ----------------------------------------------------------------------- intervals
def merge_intervals(df, by=("chrom",)):
    """Merge overlapping/adjacent intervals within groups. Returns chrom,start,end (+ group cols)."""
    out = []
    by = list(by)
    for key, g in df.sort_values(by + ["start"]).groupby(by, sort=False):
        s = g["start"].to_numpy()
        e = g["end"].to_numpy()
        cummax_e = np.maximum.accumulate(e)
        new = np.ones(len(s), dtype=bool)
        new[1:] = s[1:] > cummax_e[:-1]
        grp = np.cumsum(new) - 1
        ms = np.minimum.reduceat(s, np.flatnonzero(new))
        me = np.maximum.reduceat(e, np.flatnonzero(new))
        part = pd.DataFrame({"start": ms, "end": me})
        key = key if isinstance(key, tuple) else (key,)
        for col, val in zip(by, key):
            part[col] = val
        out.append(part)
    if not out:
        return pd.DataFrame(columns=by + ["start", "end"])
    res = pd.concat(out, ignore_index=True)
    return res[by + ["start", "end"]]


def index_by_chrom(merged):
    """merged: non-overlapping intervals -> {chrom: (starts, ends)} sorted."""
    idx = {}
    for chrom, g in merged.groupby("chrom"):
        g = g.sort_values("start")
        idx[chrom] = (g["start"].to_numpy(), g["end"].to_numpy())
    return idx


def subtract(a, b_index):
    """a: DataFrame with chrom,start,end (+any other cols, may overlap each other).
    b_index: dict from index_by_chrom of a *merged* set. Returns pieces of a not covered by b."""
    rows = []
    cols = [c for c in a.columns if c not in ("start", "end")]
    for rec in a.itertuples(index=False):
        d = rec._asdict()
        s0, e0 = d["start"], d["end"]
        if d["chrom"] not in b_index:
            rows.append(d)
            continue
        bs, be = b_index[d["chrom"]]
        i = np.searchsorted(be, s0, side="right")  # first b with end > s0
        cur = s0
        while i < len(bs) and bs[i] < e0:
            if bs[i] > cur:
                piece = {k: d[k] for k in cols}
                piece.update(start=cur, end=int(bs[i]))
                rows.append(piece)
            cur = max(cur, int(be[i]))
            i += 1
        if cur < e0:
            piece = {k: d[k] for k in cols}
            piece.update(start=cur, end=e0)
            rows.append(piece)
    if not rows:
        return a.iloc[0:0].copy()
    out = pd.DataFrame(rows)
    out = out[out["end"] > out["start"]]
    return out[list(a.columns)].reset_index(drop=True)


def overlaps_any(df, b_index):
    """Boolean array: does each interval of df overlap any interval in b_index."""
    res = np.zeros(len(df), dtype=bool)
    for chrom, g in df.groupby("chrom"):
        if chrom not in b_index:
            continue
        bs, be = b_index[chrom]
        i = np.searchsorted(be, g["start"].to_numpy(), side="right")
        ok = i < len(bs)
        hit = np.zeros(len(g), dtype=bool)
        hit[ok] = bs[i[ok]] < g["end"].to_numpy()[ok]
        res[df.index.get_indexer(g.index)] = hit
    return res


def points_to_intervals(chrom, pos, pad):
    return pd.DataFrame({"chrom": chrom, "start": np.maximum(np.asarray(pos) - pad, 0),
                         "end": np.asarray(pos) + pad})


def read_bed3(paths):
    parts = []
    for p in paths or []:
        d = pd.read_csv(p, sep="\t", header=None, comment="#", usecols=[0, 1, 2],
                        names=["chrom", "start", "end"], dtype={"chrom": str})
        parts.append(d)
    if not parts:
        return pd.DataFrame(columns=["chrom", "start", "end"])
    return pd.concat(parts, ignore_index=True)


# ------------------------------------------------------------------- tables / IO
def write_table(df, path):
    if str(path).endswith(".parquet"):
        df.to_parquet(path, index=False)
    else:
        df.to_csv(path, sep="\t", index=False)


def read_table(path):
    if str(path).endswith(".parquet"):
        return pd.read_parquet(path)
    return pd.read_csv(path, sep="\t")


def bh_fdr(p):
    p = np.asarray(p, dtype=float)
    q = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    if ok.sum() == 0:
        return q
    pv = p[ok]
    order = np.argsort(pv)
    ranked = pv[order] * len(pv) / (np.arange(len(pv)) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(pv)
    out[order] = np.minimum(ranked, 1.0)
    q[ok] = out
    return q
