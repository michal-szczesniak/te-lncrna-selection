#!/usr/bin/env python3
"""Step 2. Select TE copies of chosen subfamilies and classify them.

A copy (one RepeatMasker record) is assigned to a compartment only if it lies ENTIRELY
inside a single compartment interval from step 1. Orientation = TE strand relative to
the host lncRNA strand (sense/antisense); intergenic copies get '.'.

Supported annotation formats:
  --rmsk-format ucsc   UCSC rmsk.txt(.gz): bin swScore milliDiv milliDel milliIns genoName
                       genoStart genoEnd strand repName repClass repFamily repStart repEnd repLeft id
                       NOTE: for '-' strand UCSC stores repStart = -(bases left), repLeft = consensus
                       start. The UCSC 'id' column is not an insertion ID, so each record is its own unit.
  --rmsk-format rmout  RepeatMasker .out(.gz) (e.g. hg38.fa.out from repeatmasker.org). The ID column
                       groups fragments of one insertion into one unit (used for resampling).
Use rmout together with the matching .align file if you plan exact mapping in step 3.

--assign whole   (default) a copy enters only if it lies entirely inside one compartment interval.
--assign perbase a copy enters if >= --min-overlap of its bases fall into compartments; copy-level fields
                 describe the dominant compartment, and step 4 (--compartments) assigns every base its
                 own compartment/orientation, so the exonic part of a copy straddling an exon edge counts.

Outputs: te_copies.tsv.gz, gnomad_regions.bed (merged copy regions for subsetting VCFs),
         copy_counts.tsv (counts per subfamily x compartment x orientation)
"""
import argparse
import os
import re

import numpy as np
import pandas as pd

from ridl_utils import (index_by_chrom, log, merge_intervals, open_text, overlaps_any,
                        parse_chroms, read_bed3)

UCSC_COLS = ["bin", "swScore", "milliDiv", "milliDel", "milliIns", "chrom", "start", "end", "strand",
             "repName", "repClass", "repFamily", "repStart", "repEnd", "repLeft", "id"]


def read_rmsk_ucsc(path, chroms, want):
    rows = []
    for chunk in pd.read_csv(path, sep="\t", header=None, names=UCSC_COLS, chunksize=500_000,
                             dtype={"chrom": str, "repName": str}):
        chunk = chunk[chunk["chrom"].isin(chroms) & chunk["repName"].map(want).astype(bool)]
        rows.append(chunk)
    df = pd.concat(rows, ignore_index=True)
    plus = df["strand"] == "+"
    # sanity check of the strand convention
    bad = (plus & (df["repStart"] <= 0)) | (~plus & (df["repLeft"] <= 0))
    if bad.mean() > 0.01:
        log(f"  WARNING: {bad.mean():.1%} records violate the expected UCSC repStart/repLeft "
            "convention -- check the file format!")
    out = pd.DataFrame({
        "chrom": df["chrom"], "start": df["start"].astype(np.int64), "end": df["end"].astype(np.int64),
        "te_strand": df["strand"], "subfamily": df["repName"], "rep_class": df["repClass"],
        "rep_family": df["repFamily"],
        "cons_start": np.where(plus, df["repStart"], df["repLeft"]).astype(np.int64),
        "cons_end": df["repEnd"].astype(np.int64),
        "cons_left": np.where(plus, -df["repLeft"], -df["repStart"]).astype(np.int64),
        "div_pct": df["milliDiv"] / 10.0,
    })
    out["unit_id"] = ["u%d" % i for i in range(len(out))]
    return out


_paren = re.compile(r"^\((\d+)\)$")


def _num(x):
    m = _paren.match(x)
    return int(m.group(1)) if m else int(x)


def read_rmsk_out(path, chroms, want):
    recs = []
    with open_text(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 15 or not f[0].isdigit():
                continue
            chrom, name = f[4], f[9]
            if chrom not in chroms or not want(name):
                continue
            cls_fam = f[10].split("/")
            plus = f[8] == "+"
            if plus:
                cs, ce, cl = _num(f[11]), _num(f[12]), _num(f[13])
            else:  # 'C': (left) end begin
                cl, ce, cs = _num(f[11]), _num(f[12]), _num(f[13])
            recs.append((chrom, int(f[5]) - 1, int(f[6]), "+" if plus else "-", name, cls_fam[0],
                         cls_fam[-1], cs, ce, cl, float(f[1]), f[14]))
    out = pd.DataFrame(recs, columns=["chrom", "start", "end", "te_strand", "subfamily", "rep_class",
                                      "rep_family", "cons_start", "cons_end", "cons_left", "div_pct",
                                      "rm_id"])
    out["unit_id"] = out["chrom"] + ":" + out["rm_id"].astype(str)
    return out.drop(columns="rm_id")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rmsk", required=True)
    ap.add_argument("--rmsk-format", choices=["ucsc", "rmout"], default="ucsc")
    ap.add_argument("--subfamilies", default=None, help="comma list, e.g. AluSx,AluSx1,AluJb,MIRb,L2b")
    ap.add_argument("--subfamily-regex", default=None, help="alternative: regex on names, e.g. '^(AluS|MIR)'")
    ap.add_argument("--assign", choices=["whole", "perbase"], default="whole")
    ap.add_argument("--min-overlap", type=int, default=10, help="perbase: min bases inside compartments")
    ap.add_argument("--max-intronic", type=int, default=0,
                    help="random subsample of (dominantly) intronic copies per subfamily (0 = keep all)")
    ap.add_argument("--compartments", required=True, help="compartments.tsv.gz from step 1")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--chroms", default=None)
    ap.add_argument("--min-len", type=int, default=30)
    ap.add_argument("--exclude-bed", action="append", default=[],
                    help="drop copies overlapping these regions (e.g. ENCODE cCREs, sensitivity runs)")
    ap.add_argument("--region-flank", type=int, default=2, help="flank for gnomad_regions.bed")
    ap.add_argument("--max-intergenic", type=int, default=5000,
                    help="random subsample of intergenic copies per subfamily (0 = keep all)")
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    chroms = set(parse_chroms(a.chroms))
    if not a.subfamilies and not a.subfamily_regex:
        ap.error("give --subfamilies or --subfamily-regex")
    subfams = set(s.strip() for s in (a.subfamilies or "").split(",") if s.strip())
    rx = re.compile(a.subfamily_regex) if a.subfamily_regex else None

    def want(name):
        return name in subfams or (rx is not None and rx.search(name) is not None)

    log("reading repeat annotation ...")
    reader = read_rmsk_ucsc if a.rmsk_format == "ucsc" else read_rmsk_out
    te = reader(a.rmsk, chroms, want)
    log(f"  {len(te)} records of {te['subfamily'].nunique()} subfamilies")
    missing = subfams - set(te["subfamily"])
    if missing:
        log(f"  WARNING: no records for: {sorted(missing)}")
    te = te[(te["end"] - te["start"]) >= a.min_len].reset_index(drop=True)
    te["copy_id"] = ["c%07d" % i for i in range(len(te))]

    if a.exclude_bed:
        excl = index_by_chrom(merge_intervals(read_bed3(a.exclude_bed)))
        hit = overlaps_any(te[["chrom", "start", "end"]], excl)
        log(f"  excluding {hit.sum()} copies overlapping --exclude-bed")
        te = te[~hit].reset_index(drop=True)

    comp = pd.read_csv(a.compartments, sep="\t", dtype={"chrom": str})
    te["compartment"] = "none"
    for col in ("gene_id", "gene_name", "strand", "exon_class"):
        te[col] = "."
    te["max_tpm"] = np.nan
    te["overlap_bp"] = 0
    for chrom, g in te.groupby("chrom"):
        c = comp[comp["chrom"] == chrom].sort_values("start")
        if c.empty:
            continue
        cs, ce = c["start"].to_numpy(), c["end"].to_numpy()
        gs, ge = g["start"].to_numpy(), g["end"].to_numpy()
        if a.assign == "whole":
            i = np.searchsorted(cs, gs, side="right") - 1
            ok = (i >= 0)
            ok[ok] = ce[i[ok]] >= ge[ok]
            # reject if the next interval also starts before the copy end (overlapping compartments)
            nxt = np.clip(i + 1, 0, len(cs) - 1)
            ok &= ~((i + 1 < len(cs)) & (cs[nxt] < ge))
            best = i
            ovl = np.where(ok, ge - gs, 0)
        else:
            i0 = np.searchsorted(ce, gs, side="right")  # first interval ending after copy start
            best = np.full(len(g), -1)
            ovl = np.zeros(len(g), dtype=np.int64)
            for k in range(len(g)):
                i, top = i0[k], 0
                while i < len(cs) and cs[i] < ge[k]:
                    ov = min(ce[i], ge[k]) - max(cs[i], gs[k])
                    ovl[k] += ov
                    if ov > top:
                        top, best[k] = ov, i
                    i += 1
            ok = (best >= 0) & (ovl >= a.min_overlap)
        rows = g.index[ok]
        ci = c.iloc[best[ok]]
        te.loc[rows, "compartment"] = ci["compartment"].to_numpy()
        te.loc[rows, "overlap_bp"] = ovl[ok]
        for col in ("gene_id", "gene_name", "strand", "exon_class", "max_tpm"):
            te.loc[rows, col] = ci[col].to_numpy()
    te = te.rename(columns={"strand": "gene_strand"})
    te["orientation"] = np.where(te["compartment"].isin(["exon", "intron"]),
                                 np.where(te["te_strand"] == te["gene_strand"], "sense", "antisense"), ".")
    n_none = (te["compartment"] == "none").sum()
    te = te[te["compartment"] != "none"].reset_index(drop=True)
    log(f"  classified copies: {len(te)}  (not inside a single compartment: {n_none})")
    rng = np.random.default_rng(a.seed)
    for comp_name, cap in (("intergenic", a.max_intergenic), ("intron", a.max_intronic)):
        if cap <= 0:
            continue
        keep = np.ones(len(te), dtype=bool)
        for sf, g in te[te["compartment"] == comp_name].groupby("subfamily"):
            if len(g) > cap:
                drop = rng.choice(g.index.to_numpy(), size=len(g) - cap, replace=False)
                keep[drop] = False
        log(f"  {comp_name} subsampling: dropped {(~keep).sum()} copies")
        te = te[keep].reset_index(drop=True)

    cols = ["copy_id", "unit_id", "chrom", "start", "end", "te_strand", "subfamily", "rep_class",
            "rep_family", "cons_start", "cons_end", "cons_left", "div_pct", "compartment", "orientation",
            "gene_id", "gene_name", "gene_strand", "exon_class", "max_tpm", "overlap_bp"]
    te[cols].to_csv(os.path.join(a.out_dir, "te_copies.tsv.gz"), sep="\t", index=False)

    reg = te[["chrom", "start", "end"]].copy()
    reg["start"] = (reg["start"] - a.region_flank).clip(lower=0)
    reg["end"] = reg["end"] + a.region_flank
    merge_intervals(reg).sort_values(["chrom", "start"]).to_csv(
        os.path.join(a.out_dir, "gnomad_regions.bed"), sep="\t", index=False, header=False)

    counts = (te.assign(bp=te["end"] - te["start"])
                .groupby(["subfamily", "compartment", "orientation"])
                .agg(copies=("copy_id", "size"), units=("unit_id", "nunique"),
                     genes=("gene_id", "nunique"), bp=("bp", "sum"))
                .reset_index())
    counts.to_csv(os.path.join(a.out_dir, "copy_counts.tsv"), sep="\t", index=False)
    log("copy counts:\n" + counts.to_string(index=False))


if __name__ == "__main__":
    main()
