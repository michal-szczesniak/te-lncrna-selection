#!/usr/bin/env python3
"""How much lncRNA exonic sequence each step-1 filter removes (Supplementary Table S7).

Starts from the union of exons of the kept lncRNA genes (step-1 lnc_genes.tsv.gz, empty drop_reason) and
applies the step-1 exon filters in the same way as 01_build_compartments.py:
  splice   +-exon_splice_pad around every internal exon boundary of any annotated transcript
  tss      +-exon_tss_pad around every transcription start of any annotated transcript
  tes      +-tes_pad around every transcript end
  mask     segmental duplications (--mask-bed)
  other    spans of all other annotated genes (any type, either strand)
For each filter it reports the sequence removed when applied alone, and the sequence and genes left when
the filters are applied one after another. Also summarises gene-level selection from lnc_genes.tsv.gz.

Output: <out>.tsv (exon attrition) and <out>_genes.tsv (gene selection)
"""
import argparse

import pandas as pd

from ridl_utils import index_by_chrom, log, merge_intervals, points_to_intervals, read_bed3, read_gtf, subtract


def splice_and_ends(tx, ex):
    ex = ex.sort_values(["transcript_id", "start"])
    n_ex = ex.groupby("transcript_id")["start"].transform("size")
    first = ~ex["transcript_id"].duplicated(keep="first")
    last = ~ex["transcript_id"].duplicated(keep="last")
    multi = n_ex > 1
    s = ex.loc[multi & ~first, ["chrom", "start"]].rename(columns={"start": "pos"})
    e = ex.loc[multi & ~last, ["chrom", "end"]].rename(columns={"end": "pos"})
    splice = pd.concat([s, e], ignore_index=True).drop_duplicates()
    tss = tx["start"].where(tx["strand"] == "+", tx["end"])
    tes = tx["end"].where(tx["strand"] == "+", tx["start"])
    return splice, pd.DataFrame({"chrom": tx["chrom"], "tss": tss, "tes": tes})


def bp(df):
    return int((df["end"] - df["start"]).sum()) if len(df) else 0


def remove_other_genes(comp, all_spans):
    pieces = []
    for chrom, g in comp.groupby("chrom"):
        others = all_spans[all_spans["chrom"] == chrom]
        for gid, gg in g.groupby("gene_id"):
            o = others[others["gene_id"] != gid]
            lo, hi = gg["start"].min(), gg["end"].max()
            o = o[(o["end"] > lo) & (o["start"] < hi)]
            if len(o):
                gg = subtract(gg, index_by_chrom(merge_intervals(o[["chrom", "start", "end"]])))
            pieces.append(gg)
    return pd.concat(pieces, ignore_index=True) if pieces else comp.iloc[0:0]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--lnc-genes", required=True, help="step-1 lnc_genes.tsv.gz")
    ap.add_argument("--mask-bed", action="append", default=[])
    ap.add_argument("--exon-splice-pad", type=int, default=10)
    ap.add_argument("--exon-tss-pad", type=int, default=500)
    ap.add_argument("--tes-pad", type=int, default=100)
    ap.add_argument("--chroms", default=",".join(f"chr{i}" for i in range(1, 23)))
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    chroms = a.chroms.split(",")
    genes, tx, ex = read_gtf(a.gtf, chroms)
    lg = pd.read_csv(a.lnc_genes, sep="\t")
    lg = lg[lg["chrom"].isin(chroms)] if "chrom" in lg.columns else lg
    kept = set(lg.loc[lg["drop_reason"].fillna("") == "", "gene_id"])
    gsum = lg["drop_reason"].fillna("kept").replace("", "kept").value_counts().rename("genes").reset_index()
    gsum.columns = ["status", "genes"]
    gsum.to_csv(f"{a.out}_genes.tsv", sep="\t", index=False)
    log("gene selection (chromosomes analysed):\n" + gsum.to_string(index=False))

    kex = ex[ex["gene_id"].isin(kept)]
    start = merge_intervals(kex[["chrom", "start", "end", "gene_id"]], by=("gene_id", "chrom"))
    splice, ends = splice_and_ends(tx, ex)
    masks = read_bed3(a.mask_bed)
    masks = masks[masks["chrom"].isin(chroms)]

    def idx(d):
        return index_by_chrom(merge_intervals(d[["chrom", "start", "end"]])) if len(d) else {}
    filters = [
        ("splice sites", idx(points_to_intervals(splice["chrom"], splice["pos"], a.exon_splice_pad))),
        ("transcription starts", idx(points_to_intervals(ends["chrom"], ends["tss"], a.exon_tss_pad))),
        ("transcript ends", idx(points_to_intervals(ends["chrom"], ends["tes"], a.tes_pad))),
        ("segmental duplications", idx(masks)),
    ]
    all_spans = genes[["chrom", "start", "end", "gene_id"]]
    total = bp(start)
    rows = [{"step": "exons of kept lncRNA genes", "removed_alone_bp": 0, "left_bp": total,
             "left_share": 1.0, "genes_with_exons": start["gene_id"].nunique()}]
    cur = start
    for name, ix in filters:
        alone = total - bp(subtract(start, ix))
        cur = subtract(cur, ix)
        rows.append({"step": f"minus {name}", "removed_alone_bp": alone, "left_bp": bp(cur),
                     "left_share": bp(cur) / total, "genes_with_exons": cur["gene_id"].nunique()})
        log(f"{name}: alone -{alone:,} bp; left {bp(cur):,} bp")
    alone = total - bp(remove_other_genes(start, all_spans))
    cur = remove_other_genes(cur, all_spans)
    rows.append({"step": "minus other genes' spans", "removed_alone_bp": alone, "left_bp": bp(cur),
                 "left_share": bp(cur) / total, "genes_with_exons": cur["gene_id"].nunique()})
    other_ex = ex[["chrom", "start", "end", "gene_id"]]
    before = subtract(subtract(subtract(subtract(start, filters[0][1]), filters[1][1]), filters[2][1]), filters[3][1])
    alone = total - bp(remove_other_genes(start, other_ex))
    alt = remove_other_genes(before, other_ex)
    rows.append({"step": "ALTERNATIVE: minus other genes' exons only", "removed_alone_bp": alone,
                 "left_bp": bp(alt), "left_share": bp(alt) / total, "genes_with_exons": alt["gene_id"].nunique()})
    t = pd.DataFrame(rows)
    t.to_csv(f"{a.out}.tsv", sep="\t", index=False)
    pd.set_option("display.width", 200)
    log("exonic sequence attrition:\n" + t.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
