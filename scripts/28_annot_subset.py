#!/usr/bin/env python3
"""Revision step 28. Write a filtered copy of a step-4 annotation, with unit-level classes, for sensitivity runs
of 06_profiles.py / 12_gene_trend.py / 14_composition.py.

Added columns (all units of the selected groups are kept; filter afterwards with --query in step 6):
  unit_class      exonic units: 'boundary' if the copy's genomic span [start-pad, end+pad) contains a splice
                  site of its host gene (the TE then defines an exon edge), else 'embedded' (the exon edges lie
                  outside the TE: the classical RIDL situation). Intronic/intergenic units: 'background'.
  min_dist_unit   minimum distance (nt) of any base of the unit to a splice site of the host gene
  alu_structure   Alu groups only (AluS*, AluJ*, AluY*): consensus span of the whole copy (all compartments)
                  covers the left arm, the right arm, both ('dimer') or neither ('partial'); arms are given
                  on the group's reference consensus (--left-arm / --right-arm, default AluSx coordinates)
  phylop_consmatch  1 if the base equals the (grouped) consensus base, 0 otherwise; NaN at CpG sites and where
                  the consensus base is unknown. Named phylop_* so that step 6 and step 12 treat it as a
                  continuous layer: exon - intron > 0 means fewer substitutions since insertion in exonic copies
                  (constraint integrated over the age of the insertion; informative for young Alu elements).

Typical use (one group per call):
--relabel LABEL sets 'subfamily' of all selected bases to LABEL (subfamily_orig is kept), so that several groups
can be analysed as one class (e.g. all LTR groups) by steps 6, 12 and 14 with --subfamilies LABEL; use it with
--pop-strata ...,subfamily_orig so that the mix of subfamilies is still controlled.

  28_annot_subset.py --annot work/annot --subfamilies 'AluS*' --copies work/step3/copies_grouped.tsv.gz \
      --splice-sites work/step1_all/splice_sites.tsv.gz --query "chrom != 'chrX'" --out-dir work/annot_rev/AluS
  06_profiles.py --annot work/annot_rev/AluS --query "unit_class != 'boundary'" ...
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from ridl_utils import log

ALU_GROUPS = ("AluS*", "AluJ*", "AluY*")


def parse_range(s):
    a, b = s.split("-")
    return int(a), int(b)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annot", required=True)
    ap.add_argument("--subfamilies", required=True, help="comma list of groups (annotation 'subfamily' values) or ALL")
    ap.add_argument("--copies", required=True, help="step-3b copies_grouped.tsv.gz (copy coordinates)")
    ap.add_argument("--splice-sites", required=True, help="step-1 splice_sites.tsv.gz (chrom, pos, gene_id)")
    ap.add_argument("--query", default=None)
    ap.add_argument("--pad", type=int, default=10, help="splice sites within this distance of the copy count")
    ap.add_argument("--left-arm", default="20-110", help="consensus range that must be covered (left Alu arm)")
    ap.add_argument("--right-arm", default="150-270", help="consensus range that must be covered (right Alu arm)")
    ap.add_argument("--relabel", default=None, help="set subfamily of all selected bases to this label")
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    groups = [g for g in a.subfamilies.split(",") if g]
    la, ra = parse_range(a.left_arm), parse_range(a.right_arm)

    cp = pd.read_csv(a.copies, sep="\t", usecols=["copy_id", "chrom", "start", "end"], dtype={"chrom": str})
    ss = pd.read_csv(a.splice_sites, sep="\t", dtype={"chrom": str})
    ss_by_gene = {g: np.sort(d["pos"].to_numpy()) for g, d in ss.groupby("gene_id")}

    files = sorted(glob.glob(os.path.join(a.annot, "*.parquet"))) if os.path.isdir(a.annot) else [a.annot]
    if groups == ["ALL"]:
        groups = sorted(set().union(*[set(pq.read_table(f, columns=["subfamily"]).column(0).to_pylist())
                                      for f in files]))
        log(f"all groups: {len(groups)}")
    total = {}
    for f in files:
        d = pd.read_parquet(f, filters=[("subfamily", "in", groups)])
        if d.empty:
            continue
        if a.query:
            d = d.query(a.query, engine="python")
            if d.empty:
                continue
        # ---------------------------------------------------------------- unit class (exon edge inside the TE?)
        u = d.groupby(["unit_id", "copy_id", "compartment", "gene_id"], observed=True).agg(
            min_dist_unit=("dist_splice", "min")).reset_index()
        u = u.merge(cp, on="copy_id", how="left")
        cls = np.full(len(u), "background", dtype=object)
        ex = (u["compartment"] == "exon").to_numpy()
        for i in np.flatnonzero(ex):
            sp = ss_by_gene.get(u.at[i, "gene_id"])
            s, e = u.at[i, "start"] - a.pad, u.at[i, "end"] + a.pad
            hit = False
            if sp is not None and len(sp) and np.isfinite(s):
                k = np.searchsorted(sp, s)
                hit = k < len(sp) and sp[k] < e
            cls[i] = "boundary" if hit else "embedded"
        u["unit_class"] = cls
        d = d.merge(u[["unit_id", "copy_id", "compartment", "unit_class", "min_dist_unit"]],
                    on=["unit_id", "copy_id", "compartment"], how="left")
        # ---------------------------------------------------------------- Alu structure (whole copy)
        span = d.groupby("copy_id", observed=True)["cons_pos"].agg(["min", "max"])
        left = (span["min"] <= la[0]) & (span["max"] >= la[1])
        right = (span["min"] <= ra[0]) & (span["max"] >= ra[1])
        struct = pd.Series(np.select([left & right, left, right], ["dimer", "left", "right"], "partial"),
                           index=span.index)
        d["alu_structure"] = np.where(d["subfamily"].astype(str).isin(ALU_GROUPS),
                                      d["copy_id"].map(struct).to_numpy(), ".")
        # ---------------------------------------------------------------- consensus-match layer
        if "match" in d.columns and "cons_base" in d.columns:
            ok = (~d["cpg"].astype(bool)) & d["cons_base"].astype(str).isin(list("ACGT"))
            d["phylop_consmatch"] = np.where(ok, d["match"].astype(float), np.nan).astype(np.float32)
        else:
            log("  WARNING: no consensus columns (run step 4 with --consensus-fasta); phylop_consmatch not added")
        for c in ("unit_class", "alu_structure"):
            d[c] = d[c].astype(str)
        if a.relabel:
            d["subfamily"] = a.relabel
        d.to_parquet(os.path.join(a.out_dir, os.path.basename(f)), index=False)
        for k, v in d.drop_duplicates(["unit_id", "compartment"]).groupby(["compartment", "unit_class"]).size().items():
            total[k] = total.get(k, 0) + int(v)
        log(f"  {os.path.basename(f)}: {len(d)} bases")
    log("units by compartment and class: " + ", ".join(f"{k[0]}/{k[1]} {v}" for k, v in sorted(total.items())))
    log("done")


if __name__ == "__main__":
    main()
