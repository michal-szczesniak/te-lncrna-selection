#!/usr/bin/env python3
"""Step 4. Annotate every mapped base with sequence context, conservation and population data.

Per-base columns written (one parquet part per chromosome in --out-dir):
  copy/unit/compartment/orientation/gene info, pos, cons_pos,
  ref (genome + strand), te_base (TE orientation), tri (trinucleotide in TE orientation),
  cpg (base is part of a CpG in the current genome), cons_base/match (if --consensus-fasta),
  dist_splice (distance to nearest splice site of the host gene; exon/intron only),
  tx_relpos (relative 5'->3' position in the host reference transcript; exonic bases only),
  phylop_<name> for every --phylop name=path, depth (if --depth-tsv),
  gn_callable, gn_poly, gn_ac, gn_an, gn_singleton (gnomAD SNVs, FILTER=PASS).

Callability: a base is callable unless (i) it overlaps --mask-bed, (ii) mappability < threshold
(--mappability-bw), or (iii) AN at the site is below --an-min-frac * AN_max. AN for non-variant
sites comes from --an-tsv (tabix-indexed chrom,pos,AN; columns set by --an-cols). Without
--an-tsv only variant sites can be checked, which is weaker -- see README.

With --compartments (use with step 2 --assign perbase) every base gets its own compartment, gene,
exon_class, max_tpm and orientation from the step-1 table; bases outside all compartments are dropped.

For speed, first subset gnomAD to the regions from step 2, e.g.
  bcftools view -R gnomad_regions.bed -f PASS -v snps -Oz -o gnomad.subset.vcf.gz <gnomad_chr.vcf.bgz>
"""
import argparse
import os
from collections import defaultdict

import numpy as np
import pandas as pd
import pysam

from ridl_utils import index_by_chrom, log, merge_intervals, read_bed3, revcomp

COPY_COLS = ["copy_id", "unit_id", "subfamily", "compartment", "orientation", "gene_id", "exon_class",
             "max_tpm", "te_strand", "chrom", "start", "end", "div_pct"]


def load_consensus(path):
    seqs = {}
    if not path:
        return seqs
    with pysam.FastxFile(path) as fh:
        for rec in fh:
            seqs[rec.name.split("#")[0]] = rec.sequence.upper()
    return seqs


class Gnomad:
    def __init__(self, vcfs, an_tsv, an_cols, an_max, an_min_frac):
        self.files = {}
        for p in vcfs:
            # map each chromosome to the file whose INDEX contains it: per-chromosome files
            # (e.g. gnomAD) list all chromosomes in their headers
            try:
                contigs = pysam.TabixFile(p).contigs
            except (OSError, ValueError):
                contigs = list(pysam.VariantFile(p).header.contigs)
            for c in contigs:
                if c in self.files and self.files[c] != p:
                    log(f"  WARNING: {c} present in {self.files[c]} and {p}; using the first")
                self.files.setdefault(c, p)
        self.handles = {}
        self.an_tsv = pysam.TabixFile(an_tsv) if an_tsv else None
        self.an_cols = an_cols
        self.an_max = an_max
        self.an_min_frac = an_min_frac
        self.seen_an = []

    def _vf(self, chrom):
        p = self.files.get(chrom)
        if p is None:
            return None
        if p not in self.handles:
            self.handles[p] = pysam.VariantFile(p)
        return self.handles[p]

    def region(self, chrom, start, end):
        """Return dicts pos->(ac, an) for PASS SNVs and pos->an from the AN table (0-based pos)."""
        var = {}
        vf = self._vf(chrom)
        if vf is not None:
            for rec in vf.fetch(chrom, start, end):
                filt = list(rec.filter.keys())
                if filt and filt != ["PASS"]:
                    continue
                if len(rec.ref) != 1 or any(alt is None or len(alt) != 1 for alt in rec.alts or []):
                    continue
                ac = rec.info.get("AC", 0)
                ac = int(sum(ac)) if isinstance(ac, (tuple, list)) else int(ac or 0)
                an = rec.info.get("AN", None)
                p0 = rec.pos - 1
                prev = var.get(p0, (0, an))
                var[p0] = (prev[0] + ac, an if an is not None else prev[1])
                if an is not None:
                    self.seen_an.append(an)
        an_map = {}
        if self.an_tsv is not None:
            cpos, can = self.an_cols
            try:
                for row in self.an_tsv.fetch(chrom, start, end, parser=pysam.asTuple()):
                    an_map[int(row[cpos]) - 1] = int(row[can])
            except ValueError:
                pass
        return var, an_map


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bases", required=True, help="parquet from step 3")
    ap.add_argument("--copies", required=True, help="te_copies.tsv.gz from step 2")
    ap.add_argument("--genome-fasta", required=True, help="indexed FASTA (samtools faidx)")
    ap.add_argument("--splice-sites", required=True)
    ap.add_argument("--ref-tx", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--consensus-fasta", default=None)
    ap.add_argument("--phylop", action="append", default=[], help="name=path.bw (repeatable)")
    ap.add_argument("--depth-tsv", default=None, help="chrom,pos(0-based),depth[,...] from helper_maf_depth.py")
    ap.add_argument("--gnomad-vcf", action="append", default=[], help="bgzipped+indexed VCF (repeatable)")
    ap.add_argument("--an-tsv", default=None)
    ap.add_argument("--an-cols", default="1,2", help="0-based column indices of pos,AN in --an-tsv")
    ap.add_argument("--an-max", type=float, default=None, help="max AN (default: max observed)")
    ap.add_argument("--an-min-frac", type=float, default=0.8)
    ap.add_argument("--mask-bed", action="append", default=[])
    ap.add_argument("--mappability-bw", default=None)
    ap.add_argument("--min-mappability", type=float, default=1.0)
    ap.add_argument("--compartments", default=None, help="per-base compartment assignment (step-1 table)")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)

    bases = pd.read_parquet(a.bases)
    cp = pd.read_csv(a.copies, sep="\t", dtype={"chrom": str})
    cp = cp[COPY_COLS + [c for c in ("subfamily_orig",) if c in cp.columns]]
    comp_by_chrom = {}
    if a.compartments:
        comp = pd.read_csv(a.compartments, sep="\t", dtype={"chrom": str})
        comp_by_chrom = {c: g.sort_values("start").reset_index(drop=True) for c, g in comp.groupby("chrom")}
    cons = load_consensus(a.consensus_fasta)
    fa = pysam.FastaFile(a.genome_fasta)

    import pyBigWig
    bws = {}
    for item in a.phylop:
        name, path = item.split("=", 1)
        bws[name] = pyBigWig.open(path)
    mbw = pyBigWig.open(a.mappability_bw) if a.mappability_bw else None
    masks = read_bed3(a.mask_bed)
    mask_idx = index_by_chrom(merge_intervals(masks)) if len(masks) else {}
    gn = Gnomad(a.gnomad_vcf, a.an_tsv, tuple(int(x) for x in a.an_cols.split(",")), a.an_max,
                a.an_min_frac) if a.gnomad_vcf else None
    depth = None
    if a.depth_tsv:
        depth = pd.read_csv(a.depth_tsv, sep="\t", dtype={"chrom": str})

    ss = pd.read_csv(a.splice_sites, sep="\t", dtype={"chrom": str})
    ss_by_gene = {g: np.sort(d["pos"].to_numpy()) for g, d in ss.groupby("gene_id")}
    rtx = pd.read_csv(a.ref_tx, sep="\t", dtype={"chrom": str})
    rtx_by_gene = {g: d.sort_values("start") for g, d in rtx.groupby("gene_id")}

    for chrom in cp["chrom"].unique():
        cpc = cp[cp["chrom"] == chrom]
        B = bases[bases["copy_id"].isin(set(cpc["copy_id"]))].merge(cpc, on="copy_id", how="inner")
        if B.empty:
            continue
        B = B.sort_values("pos").reset_index(drop=True)
        if a.compartments:
            c = comp_by_chrom.get(chrom)
            if c is None:
                continue
            p = B["pos"].to_numpy()
            k = np.searchsorted(c["start"].to_numpy(), p, side="right") - 1
            inside = (k >= 0) & (p < c["end"].to_numpy()[np.clip(k, 0, len(c) - 1)])
            B = B[inside].reset_index(drop=True)
            ci = c.iloc[k[inside]].reset_index(drop=True)
            for col in ("compartment", "gene_id", "exon_class", "max_tpm"):
                B[col] = ci[col].to_numpy()
            gstr = ci["strand"].to_numpy()
            B["orientation"] = np.where(np.isin(B["compartment"].to_numpy(), ["exon", "intron"]),
                                        np.where(B["te_strand"].to_numpy() == gstr, "sense", "antisense"), ".")
            if B.empty:
                continue
        log(f"  {chrom}: {len(B)} bases")
        n = len(B)
        ref = np.empty(n, dtype=object)
        tri = np.empty(n, dtype=object)
        cpg = np.zeros(n, dtype=bool)
        ph = {k: np.full(n, np.nan, dtype=np.float32) for k in bws}
        mapq = np.ones(n, dtype=bool)
        gn_ac = np.zeros(n, dtype=np.int64)
        gn_an = np.full(n, np.nan)
        an_site = np.full(n, np.nan)

        for cid, idx in B.groupby("copy_id").indices.items():
            r0 = B.iloc[idx[0]]
            s, e = int(r0["start"]), int(r0["end"])
            seq = fa.fetch(chrom, max(s - 1, 0), e + 1).upper()
            off = 1 if s > 0 else 0
            pos = B["pos"].to_numpy()[idx]
            rel = pos - s + off
            ref[idx] = [seq[i] for i in rel]
            tri_plus = [seq[i - 1:i + 2] if 0 < i < len(seq) - 1 else "NNN" for i in rel]
            if r0["te_strand"] == "-":
                tri[idx] = [revcomp(t) for t in tri_plus]
            else:
                tri[idx] = tri_plus
            cpg[idx] = [(seq[i] == "C" and i + 1 < len(seq) and seq[i + 1] == "G") or
                        (seq[i] == "G" and i > 0 and seq[i - 1] == "C") for i in rel]
            for k, bw in ph.items():
                try:
                    vals = np.array(bws[k].values(chrom, s, e, numpy=True), dtype=np.float32)
                    bw[idx] = vals[pos - s]
                except RuntimeError:
                    pass  # chromosome/interval absent from the bigWig -> NaN
            if mbw is not None:
                mv = np.array(mbw.values(chrom, s, e, numpy=True), dtype=float)
                mapq[idx] = np.nan_to_num(mv[pos - s], nan=0.0) >= a.min_mappability
            if gn is not None:
                var, an_map = gn.region(chrom, s, e)
                for j, p in zip(idx, pos):
                    if p in var:
                        gn_ac[j] = var[p][0]
                        if var[p][1] is not None:
                            gn_an[j] = var[p][1]
                    if p in an_map:
                        an_site[j] = an_map[p]

        B["ref"] = ref
        B["te_base"] = np.where(B["te_strand"] == "-", [revcomp(x) for x in ref], ref)
        B["tri"] = tri
        B["cpg"] = cpg
        if cons:
            cb = [cons.get(sf, "")[cp_ - 1] if 0 < cp_ <= len(cons.get(sf, "")) else "N"
                  for sf, cp_ in zip(B["subfamily"], B["cons_pos"])]
            B["cons_base"] = cb
            B["match"] = B["cons_base"] == B["te_base"]

        # distance to nearest splice site of the host gene
        dist = np.full(n, np.nan)
        relpos = np.full(n, np.nan)
        for gid, idx in B.groupby("gene_id").indices.items():
            if gid == ".":
                continue
            pos = B["pos"].to_numpy()[idx]
            sp = ss_by_gene.get(gid)
            if sp is not None and len(sp):
                j = np.searchsorted(sp, pos)
                left = np.abs(pos - sp[np.clip(j - 1, 0, len(sp) - 1)])
                right = np.abs(sp[np.clip(j, 0, len(sp) - 1)] - pos)
                dist[idx] = np.minimum(left, right)
            rt = rtx_by_gene.get(gid)
            if rt is not None:
                ex_idx = idx[(B["compartment"].to_numpy()[idx] == "exon")]
                if len(ex_idx):
                    p = B["pos"].to_numpy()[ex_idx]
                    st, en = rt["start"].to_numpy(), rt["end"].to_numpy()
                    k = np.searchsorted(st, p, side="right") - 1
                    ok = (k >= 0) & (p < en[np.clip(k, 0, len(en) - 1)])
                    kk = np.clip(k, 0, len(en) - 1)
                    plus = rt["strand"].iloc[0] == "+"
                    within = np.where(plus, p - st[kk], en[kk] - 1 - p)
                    txc = rt["tx_offset"].to_numpy()[kk] + within
                    relpos[ex_idx[ok]] = txc[ok] / rt["tx_len"].iloc[0]
        B["dist_splice"] = dist
        B["tx_relpos"] = relpos
        for k, v in ph.items():
            B[f"phylop_{k}"] = v
        if depth is not None:
            d = depth[depth["chrom"] == chrom].drop(columns="chrom")
            B = B.merge(d, on="pos", how="left")

        # callability
        callable_ = mapq.copy()
        if chrom in mask_idx:
            ms, me = mask_idx[chrom]
            j = np.searchsorted(me, B["pos"].to_numpy(), side="right")
            ok = j < len(ms)
            inmask = np.zeros(n, dtype=bool)
            inmask[ok] = ms[j[ok]] <= B["pos"].to_numpy()[ok]
            callable_ &= ~inmask
        if gn is not None:
            an_max = a.an_max or (max(gn.seen_an) if gn.seen_an else np.nan)
            an_use = np.where(np.isnan(an_site), gn_an, an_site)
            low = an_use < a.an_min_frac * an_max  # NaN -> False (unknown treated as callable)
            callable_ &= ~low
            if gn.an_tsv is not None:
                callable_ &= ~np.isnan(an_site)  # with an AN table, sites absent from it are uncallable
            B["gn_ac"] = gn_ac
            B["gn_an"] = an_use
            B["gn_poly"] = gn_ac > 0
            B["gn_singleton"] = gn_ac == 1
        B["gn_callable"] = callable_
        msg = f"  {chrom}: callable {callable_.mean():.1%}"
        if gn is not None:
            if chrom not in gn.files:
                log(f"  WARNING: no gnomAD file contains {chrom}")
            msg += f", polymorphic among callable {B.loc[callable_, 'gn_poly'].mean():.1%}"
        log(msg)
        B.to_parquet(os.path.join(a.out_dir, f"part-{chrom}.parquet"), index=False)
    log("done")


if __name__ == "__main__":
    main()
