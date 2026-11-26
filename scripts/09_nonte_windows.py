#!/usr/bin/env python3
"""Control. Windows of ordinary (non-repeat) sequence in the same compartments.

Tiles lncRNA exonic and intronic compartments (and intergenic regions) into windows of --window nt that
overlap no RepeatMasker annotation of any class, and writes them in the step-2 copies format as
subfamily 'nonTE' (orientation 'sense' for exon/intron windows, consensus positions 1..window).
Then run: step 3 --mode linear, gnomAD subsetting with the written gnomad_regions.bed, step 4 with
--compartments, and step 6 with --subfamilies nonTE --pop-strata tri,cpg --bin-size <window>.
The exon-vs-intron (and intron-vs-intergenic) difference of ordinary lncRNA sequence is the baseline
against which the whole-element TE results must be judged.

Outputs: te_copies.tsv.gz (windows), gnomad_regions.bed, window_counts.tsv
"""
import argparse
import os

import numpy as np
import pandas as pd

from ridl_utils import (index_by_chrom, log, merge_intervals, open_text, parse_chroms, subtract)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--compartments", required=True)
    ap.add_argument("--rmsk-out", required=True, help="RepeatMasker .out(.gz): ALL records are masked")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--window", type=int, default=100)
    ap.add_argument("--max-per-comp", type=int, default=20000, help="random windows kept per compartment")
    ap.add_argument("--compartment-set", default="exon,intron,intergenic")
    ap.add_argument("--chroms", default=None)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    chroms = set(parse_chroms(a.chroms))
    keep = a.compartment_set.split(",")
    rng = np.random.default_rng(a.seed)

    comp = pd.read_csv(a.compartments, sep="\t", dtype={"chrom": str})
    comp = comp[comp["compartment"].isin(keep) & comp["chrom"].isin(chroms)].reset_index(drop=True)

    log("reading all repeat annotations ...")
    rows = []
    with open_text(a.rmsk_out) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 15 or not f[0].isdigit() or f[4] not in chroms:
                continue
            rows.append((f[4], int(f[5]) - 1, int(f[6])))
    rep = pd.DataFrame(rows, columns=["chrom", "start", "end"])
    log(f"  {len(rep)} repeat records")
    pieces = subtract(comp, index_by_chrom(merge_intervals(rep)))
    pieces = pieces[(pieces["end"] - pieces["start"]) >= a.window].reset_index(drop=True)

    W = a.window
    n = ((pieces["end"] - pieces["start"]) // W).to_numpy()
    rep_idx = np.repeat(np.arange(len(pieces)), n)
    k = np.arange(n.sum()) - np.repeat(np.cumsum(n) - n, n)
    win = pieces.iloc[rep_idx].reset_index(drop=True)
    win["start"] = win["start"].to_numpy() + k * W
    win["end"] = win["start"] + W
    out = []
    for c, g in win.groupby("compartment"):
        if len(g) > a.max_per_comp:
            g = g.iloc[np.sort(rng.choice(len(g), a.max_per_comp, replace=False))]
        out.append(g)
        log(f"  {c}: {len(win[win['compartment'] == c])} windows available, {len(g)} kept")
    win = pd.concat(out, ignore_index=True).sort_values(["chrom", "start"]).reset_index(drop=True)

    genic = win["compartment"].isin(["exon", "intron"]).to_numpy()
    res = pd.DataFrame({
        "copy_id": ["w%07d" % i for i in range(len(win))],
        "chrom": win["chrom"], "start": win["start"], "end": win["end"],
        "te_strand": np.where(genic, win["strand"], "+"),
        "subfamily": "nonTE", "rep_class": "none", "rep_family": "none",
        "cons_start": 1, "cons_end": W, "cons_left": 0, "div_pct": 0.0,
        "compartment": win["compartment"], "orientation": np.where(genic, "sense", "."),
        "gene_id": win["gene_id"], "gene_name": win["gene_name"], "gene_strand": win["strand"],
        "exon_class": win["exon_class"], "max_tpm": win["max_tpm"], "overlap_bp": W})
    res.insert(1, "unit_id", res["copy_id"])
    res.to_csv(os.path.join(a.out_dir, "te_copies.tsv.gz"), sep="\t", index=False)
    merge_intervals(res[["chrom", "start", "end"]]).sort_values(["chrom", "start"]).to_csv(
        os.path.join(a.out_dir, "gnomad_regions.bed"), sep="\t", index=False, header=False)
    cnt = res.groupby("compartment").agg(windows=("copy_id", "size"), genes=("gene_id", "nunique"))
    cnt.to_csv(os.path.join(a.out_dir, "window_counts.tsv"), sep="\t")
    log("windows written:\n" + cnt.to_string())


if __name__ == "__main__":
    main()
