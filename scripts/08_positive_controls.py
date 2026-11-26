#!/usr/bin/env python3
"""Step 8. Positive-control case studies.

--controls is a BED (chrom, start, end, name) with genomic coordinates of TE fragments with known
function (e.g. the MIR in DRAIC exon 5, the TE in CARMEN, ANRIL RIDLs, SIRLOIN-containing Alus).
For each control the script reports which annotated TE bases (step 4) it covers, grouped by copy,
compartment and orientation: the consensus positions / bins, the profile values of those bins for
every layer, and whether any covered bin lies in a peak from step 6.

Note: a control whose subfamily was not selected in step 2 will not be found -- add it there.
"""
import argparse

import numpy as np
import pandas as pd

from ridl_stats import load_annot
from ridl_utils import log


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--controls", required=True, help="BED: chrom start end name")
    ap.add_argument("--annot", required=True, help="annotation directory from step 4")
    ap.add_argument("--profiles", required=True, help="profiles.tsv from step 6")
    ap.add_argument("--peaks", required=True, help="peaks.tsv from step 6")
    ap.add_argument("--bin-size", type=int, default=10)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    ctl = pd.read_csv(a.controls, sep="\t", header=None, usecols=[0, 1, 2, 3],
                      names=["chrom", "start", "end", "name"], dtype={"chrom": str})
    ann = load_annot(a.annot, columns=["copy_id", "chrom", "pos", "cons_pos", "subfamily", "compartment",
                                       "orientation", "gene_id"])
    prof = pd.read_csv(a.profiles, sep="\t")
    try:
        peaks = pd.read_csv(a.peaks, sep="\t")
    except pd.errors.EmptyDataError:
        peaks = pd.DataFrame(columns=["layer", "subfamily", "orientation", "bin_from", "bin_to"])
    if peaks.empty:
        peaks = pd.DataFrame(columns=["layer", "subfamily", "orientation", "bin_from", "bin_to",
                                      "cons_from", "cons_to"])

    out = []
    for c in ctl.itertuples(index=False):
        hit = ann[(ann["chrom"] == c.chrom) & (ann["pos"] >= c.start) & (ann["pos"] < c.end)]
        if hit.empty:
            out.append({"control": c.name, "status": "no annotated TE bases (subfamily not selected?)"})
            continue
        for (cid, sf, comp, ori, gid), b in hit.groupby(["copy_id", "subfamily", "compartment", "orientation",
                                                          "gene_id"]):
            bins = sorted(set(((b["cons_pos"] - 1) // a.bin_size).tolist()))
            rec = {"control": c.name, "status": "ok", "copy_id": cid, "subfamily": sf, "compartment": comp,
                   "orientation": ori, "gene_id": gid, "n_bases": len(b),
                   "cons_from": int(b["cons_pos"].min()), "cons_to": int(b["cons_pos"].max()),
                   "bins": ",".join(map(str, bins))}
            p = prof[(prof["subfamily"] == sf) & (prof["orientation"] == ori) & prof["bin"].isin(bins)]
            for layer, g in p.groupby("layer"):
                rec[f"{layer}_mean_effect"] = float(np.nanmean(g["effect"])) if g["effect"].notna().any() else np.nan
                rec[f"{layer}_min_q"] = float(np.nanmin(g["q"])) if g["q"].notna().any() else np.nan
            pk = peaks[(peaks["subfamily"] == sf) & (peaks["orientation"] == ori)]
            in_peak = [f"{r.layer}:{r.cons_from}-{r.cons_to}" for r in pk.itertuples(index=False)
                       if any(r.bin_from <= x <= r.bin_to for x in bins)]
            rec["in_peaks"] = ";".join(in_peak) if in_peak else ""
            out.append(rec)
    res = pd.DataFrame(out)
    res.to_csv(a.out, sep="\t", index=False)
    log(res.to_string(index=False))


if __name__ == "__main__":
    main()
