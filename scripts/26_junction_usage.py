#!/usr/bin/env python3
"""Revision step 26. How much is each annotated splice junction used? GTEx v10 junction read counts for the
introns of the splice-site test (or any list of introns).

Input: GTEx junction table (GCT: '#1.2', dimensions, header 'Name Description <samples>') streamed once;
optional sample attributes (SAMPID, SMTSD) for per-tissue means. Introns come from the site ids of
19_splice_analyze.py positions ('chrom:start-end<strand><d|a>', 0-based half-open intron) or from --introns
(chrom, start, end, strand).
Junction names are parsed leniently: chr1_100_200, chr1:100-200, chr1:100:200, with or without the 'chr'
prefix and with an optional strand ('_+', ':+', ':1', ':2'); if the Name column does not parse, the
Description column is tried. The coordinate convention (first/last intron base, exon ends, 0- or 1-based) is
checked against the data: all offsets of -1, 0, +1 on both ends are tried and the one with most matches is
used. If nothing matches, the script stops with example names instead of writing an empty table.

Per intron: total reads, mean reads per sample, samples with >= 1 and >= 5 reads, highest tissue mean and
number of tissues with mean >= 1 read (with attributes), and relative use: the share of reads of all junctions
that share this intron's donor (donor_share) or acceptor (acceptor_share), i.e. how often the splice site is
used by this junction rather than by alternatives (PSI-like, independent of expression level).

Output: <out>.tsv.gz
"""
import argparse
import gzip
import re
import sys

import numpy as np
import pandas as pd

from ridl_utils import log

JRE = re.compile(r"^(?:chr)?(?P<c>[0-9]{1,2}|X|Y|M|MT)[_:](?P<s>\d+)[-_:](?P<e>\d+)(?:[_:](?P<st>[-+012]))?$")


def parse_junction(name, desc):
    """(chrom, a, b, strand or None) from a junction name, trying Name, then Description; None if neither."""
    for x in (name, desc):
        m = JRE.match(x.strip())
        if m:
            c = m.group("c")
            c = "chrM" if c == "MT" else f"chr{c}"
            st = {"1": "+", "2": "-", "+": "+", "-": "-"}.get(m.group("st"))
            return c, int(m.group("s")), int(m.group("e")), st
    return None


def introns_from_positions(path):
    p = pd.read_csv(path, sep="\t", usecols=["set", "site_id"], low_memory=False)
    p = p[p["set"].isin(["lnc_TE_site", "lnc_nonTE_site", "pc_site", "lnc_nonTE_flank", "pc_flank"])]
    sid = p["site_id"].drop_duplicates().str[:-1]  # drop d/a suffix
    m = sid.str.extract(r"^(?P<chrom>[^:]+):(?P<start>\d+)-(?P<end>\d+)(?P<strand>[+-])$")
    m = m.dropna().drop_duplicates()
    m["start"] = m["start"].astype(np.int64)
    m["end"] = m["end"].astype(np.int64)
    return m.reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gct", required=True)
    ap.add_argument("--attributes", default=None, help="GTEx SampleAttributesDS.txt (SAMPID, SMTSD)")
    ap.add_argument("--positions", default=None, help="19_splice_analyze.py <out>_positions.tsv.gz")
    ap.add_argument("--introns", default=None, help="TSV with chrom, start, end, strand (0-based half-open)")
    ap.add_argument("--max-rows", type=int, default=0, help="stop after this many rows (testing)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    if a.positions:
        it = introns_from_positions(a.positions)
    elif a.introns:
        it = pd.read_csv(a.introns, sep="\t", dtype={"chrom": str})
    else:
        raise SystemExit("give --positions or --introns")
    log(f"introns to look up: {len(it)}")
    it["chrom"] = it["chrom"].astype(str)
    sgn = lambda d: f"{d:+d}" if d else ""
    conv = {f"s{sgn(ds)},e{sgn(de)}": (ds, de) for ds in (1, 0, -1) for de in (0, 1, -1)}
    keys = {name: {} for name in conv}
    for name, (ds, de) in conv.items():
        for i, (c, s, e) in enumerate(zip(it["chrom"], it["start"], it["end"])):
            keys[name].setdefault((c, int(s) + ds, int(e) + de), []).append(i)
    strand_of = it["strand"].to_numpy()

    with gzip.open(a.gct, "rt") as fh:
        first = fh.readline()
        if not first.startswith("#"):
            raise SystemExit("not a GCT file")
        fh.readline()
        header = fh.readline().rstrip("\n").split("\t")
        samples = header[2:]
        log(f"samples: {len(samples)}")
        tissue_of = None
        if a.attributes:
            at = pd.read_csv(a.attributes, sep="\t", usecols=["SAMPID", "SMTSD"], dtype=str)
            tmap = dict(zip(at["SAMPID"], at["SMTSD"]))
            tissue_of = pd.Series([tmap.get(s, "NA") for s in samples])
            tissues, t_idx = np.unique(tissue_of.to_numpy(), return_inverse=True)
            t_n = np.bincount(t_idx, minlength=len(tissues))
            log(f"tissues: {len(tissues)} (samples without tissue: {(tissue_of == 'NA').sum()})")
        # left/right read totals per (chrom, coordinate, strand or None), for the PSI-like shares
        left, right = {}, {}
        vec = {name: {} for name in conv}  # per convention: intron index -> summed read vector, junction names
        # only rows that share a coordinate with one of our introns are parsed (shares need their totals)
        want_l = {(c, int(s) + d) for c, s in zip(it["chrom"], it["start"]) for d in (-1, 0, 1)}
        want_r = {(c, int(e) + d) for c, e in zip(it["chrom"], it["end"]) for d in (-1, 0, 1)}
        n = n_unparsed = n_hit = 0
        examples = []
        for line in fh:
            n += 1
            name, desc, rest = line.rstrip("\n").split("\t", 2)
            if len(examples) < 5:
                examples.append(f"{name} | {desc}")
            p = parse_junction(name, desc)
            if p is None:
                n_unparsed += 1
                continue
            c, s1, e1, st = p
            hit_l, hit_r = (c, s1) in want_l, (c, e1) in want_r
            if not (hit_l or hit_r):
                continue
            n_hit += 1
            v = np.fromstring(rest, dtype=float, sep="\t")
            if len(v) != len(samples):
                v = np.array(rest.split("\t"), dtype=float)
            tot = float(v.sum())
            if hit_l:
                left[(c, s1, st)] = left.get((c, s1, st), 0) + tot
            if hit_r:
                right[(c, e1, st)] = right.get((c, e1, st), 0) + tot
            for cn, kd in keys.items():
                for i in kd.get((c, s1, e1), ()):
                    if st is not None and st != strand_of[i]:
                        continue
                    if i in vec[cn]:
                        vec[cn][i][0] += v
                        vec[cn][i][1].append(name)
                    else:
                        vec[cn][i] = [v.copy(), [name]]
            if n % 200000 == 0:
                log(f"  {n} rows")
            if a.max_rows and n >= a.max_rows:
                break
    log("first junction names (Name | Description): " + "; ".join(examples))
    log(f"rows read: {n}; names not parsed: {n_unparsed}; rows sharing a coordinate with our introns: {n_hit}")
    log("matches by convention: " + ", ".join(f"{k} {len(v)}" for k, v in vec.items()))
    best = max(vec, key=lambda k: len(vec[k]))
    if not vec[best]:
        sys.exit("ERROR: no intron matched any junction; check the name format above and the intron coordinates")
    ds, de = conv[best]
    log(f"using convention {best} (intron start {ds:+d}, end {de:+d} relative to 0-based half-open)")
    found = {}
    for i, (v, names) in vec[best].items():
        tot = float(v.sum())
        r = {"gtex_name": ",".join(names), "total_reads": tot, "mean_reads": tot / len(v),
             "samples_ge1": int((v >= 1).sum()), "samples_ge5": int((v >= 5).sum())}
        if tissue_of is not None:
            tm = np.bincount(t_idx, weights=v, minlength=len(tissues)) / np.maximum(t_n, 1)
            r["max_tissue_mean"] = float(tm.max())
            r["max_tissue"] = tissues[int(tm.argmax())]
            r["tissues_mean_ge1"] = int((tm >= 1).sum())
        found[i] = r
    F = pd.DataFrame.from_dict(found, orient="index")
    out = it.join(F, how="left")
    out["matched"] = out["gtex_name"].notna()
    for c in ("total_reads", "samples_ge1", "samples_ge5"):
        out[c] = out[c].fillna(0)
    out["mean_reads"] = out["mean_reads"].fillna(0.0)
    if "max_tissue_mean" in out:
        out["max_tissue_mean"] = out["max_tissue_mean"].fillna(0.0)
        out["tissues_mean_ge1"] = out["tissues_mean_ge1"].fillna(0)
    # totals over unstranded junctions and junctions on the intron's strand
    ls = np.array([left.get((c, int(s) + ds, None), 0) + left.get((c, int(s) + ds, st), 0)
                   for c, s, st in zip(out["chrom"], out["start"], out["strand"])], dtype=float)
    rs = np.array([right.get((c, int(e) + de, None), 0) + right.get((c, int(e) + de, st), 0)
                   for c, e, st in zip(out["chrom"], out["end"], out["strand"])], dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        lsh = np.where(ls > 0, out["total_reads"] / ls, np.nan)
        rsh = np.where(rs > 0, out["total_reads"] / rs, np.nan)
    plus = out["strand"].to_numpy() == "+"
    out["donor_share"] = np.where(plus, lsh, rsh)
    out["acceptor_share"] = np.where(plus, rsh, lsh)
    out["intron_key"] = out["chrom"] + ":" + out["start"].astype(str) + "-" + out["end"].astype(str) + out["strand"]
    out.to_csv(a.out if a.out.endswith(".gz") else f"{a.out}.tsv.gz", sep="\t", index=False)
    log(f"matched introns: {out['matched'].mean():.1%}; median total reads of matched: "
        f"{out.loc[out['matched'], 'total_reads'].median():.0f}")


if __name__ == "__main__":
    main()
