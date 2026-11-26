#!/usr/bin/env python3
"""Step 3. Map each genomic base of each TE copy to a 1-based consensus position.

Modes
  linear  Interpolate genomic positions onto the consensus span [cons_start, cons_end]
          (respecting strand). Fast, no extra files; small indels blur positions by a few nt,
          acceptable for ~10-nt bins in a first pass.
  align   Exact mapping from a RepeatMasker .align(.gz) file (e.g. hg38.fa.align.gz from
          repeatmasker.org). Must come from the SAME RepeatMasker run as the annotation used in
          step 2 (use --rmsk-format rmout with the matching .out file); copies are matched by
          chrom/start/end/subfamily and the match rate is reported. Bases inserted relative to the
          consensus get no consensus position and are dropped.

Output: parquet with copy_id, pos (0-based genomic), cons_pos (1-based).
Use --compare-to <other_output.parquet> to report agreement between two mappings.
In align mode, --consensus-out writes each selected subfamily's consensus reconstructed from the
consensus lines of its alignments (majority base per position; needed by step 3b for grouping).
Alignments are matched to copies by chromosome, subfamily and OVERLAP (RepeatMasker's .out intervals
can be a few-dozen nt longer/shorter than the alignment); only bases inside the copy are kept.
Split LINE models in the .align (e.g. L2b_5end / L2b_orf2 / L2b_3end, reported as 'L2b' in the .out)
are matched to the base name; only the model given by --split-model (default 3end) is used, so that
all positions of a subfamily share one consensus coordinate system.
"""
import argparse
import re

import numpy as np
import pandas as pd

from ridl_utils import log, open_text

HEADER = re.compile(r"^\s*(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+(\S+)\s+(\d+)\s+(\d+)\s+\((\d+)\)\s+(C\s+)?(\S+)\s+")
ALN = re.compile(r"^(C\s+|\s+)?(\S+)\s+(\d+)\s+([A-Za-z\-*?]+)\s+(\d+)\s*$")


def map_linear(cp):
    L = (cp["end"] - cp["start"]).to_numpy()
    cs = cp["cons_start"].to_numpy()
    ce = cp["cons_end"].to_numpy()
    rep = np.repeat(np.arange(len(cp)), L)
    k = np.arange(L.sum()) - np.repeat(np.cumsum(L) - L, L)  # 0..L-1 within copy
    Lg = L[rep]
    Lc = (ce - cs + 1)[rep]
    frac = np.where(Lg > 1, k / np.maximum(Lg - 1, 1), 0.0)
    off = np.rint(frac * (Lc - 1)).astype(np.int64)
    plus = (cp["te_strand"].to_numpy() == "+")[rep]
    cons = np.where(plus, cs[rep] + off, ce[rep] - off)
    pos = cp["start"].to_numpy()[rep] + k
    return pd.DataFrame({"copy_id": cp["copy_id"].to_numpy()[rep], "pos": pos.astype(np.int64),
                         "cons_pos": cons.astype(np.int32)})


def iter_align_blocks(path, wanted=None):
    """Yield (header_match, query_seq, query_start, subj_seq, subj_start) for every alignment
    (only for repeat names in `wanted`, if given; other blocks are skipped cheaply)."""
    hdr, q, s, qs, ss, turn = None, [], [], None, None, 0
    with open_text(path) as fh:
        for line in fh:
            if line[:1].isdigit():
                m = HEADER.match(line)
                if m:
                    if hdr is not None and q:
                        yield hdr, "".join(q), qs, "".join(s), ss
                    hdr, q, s, qs, ss, turn = m, [], [], None, None, 0
                    if wanted is not None and m.group(10).split("#")[0] not in wanted:
                        hdr = None
                    continue
            if hdr is None:
                continue
            if line.startswith(("Matrix", "Kimura", "Transitions", "Gap_init")):
                continue
            am = ALN.match(line.rstrip("\n"))
            if not am:
                continue
            if turn == 0:
                if qs is None:
                    qs = int(am.group(3))
                q.append(am.group(4))
                turn = 1
            else:
                if ss is None:
                    ss = int(am.group(3))
                s.append(am.group(4))
                turn = 0
    if hdr is not None and q:
        yield hdr, "".join(q), qs, "".join(s), ss


BASE_IDX = np.full(256, -1, dtype=np.int64)
for _i, _b in enumerate(b"ACGT"):
    BASE_IDX[_b] = _i
    BASE_IDX[ord(chr(_b).lower())] = _i
COMPL_IDX = np.array([3, 2, 1, 0])


SPLIT = re.compile(r"^(.+)_(5end|orf2|3end)$")


def base_name(raw, split_model):
    """Map an .align repeat name to the .out subfamily name; None if it is a split model we skip."""
    m = SPLIT.match(raw)
    if m:
        return m.group(1) if m.group(2) == split_model else None
    return raw


def map_align(cp, align_path, cons_counts=None, cons_max=3000, split_model="3end"):
    idx = {}
    for (chrom, sf), g in cp.groupby(["chrom", "subfamily"]):
        g = g.sort_values("start")
        idx[(chrom, sf)] = (g["start"].to_numpy(), g["end"].to_numpy(), g["copy_id"].to_numpy())
    wanted = set(cp["subfamily"])
    wanted_raw = set(wanted) | {f"{w}_{split_model}" for w in wanted}
    out_c, out_p, out_r = [], [], []
    seen_copies = set()
    n_cons = {}
    for hdr, qseq, qstart, sseq, sstart in iter_align_blocks(align_path, wanted_raw):
        name = base_name(hdr.group(10).split("#")[0], split_model)
        if name not in wanted:
            continue
        comp = hdr.group(9) is not None
        step = -1 if comp else 1
        chrom, qb, qe = hdr.group(5), int(hdr.group(6)), int(hdr.group(7))
        cid, cstart, cend = None, None, None
        ent = idx.get((chrom, name))
        if ent is not None:
            ss_, ee_, ids_ = ent
            i = np.searchsorted(ee_, qb - 1, side="right")
            best = 0
            while i < len(ss_) and ss_[i] < qe:
                ov = min(ee_[i], qe) - max(ss_[i], qb - 1)
                if ov > best:
                    best, cid, cstart, cend = ov, ids_[i], ss_[i], ee_[i]
                i += 1
        want_cons = cons_counts is not None and n_cons.get(name, 0) < cons_max
        if cid is None and not want_cons:
            continue
        qa = np.frombuffer(qseq.encode(), dtype="S1")
        sa = np.frombuffer(sseq.encode(), dtype="S1")
        if len(qa) != len(sa):
            log(f"  WARNING: unequal alignment lengths ({chrom}:{qb}-{qe} {name}), skipped")
            continue
        qgap = qa == b"-"
        sgap = sa == b"-"
        qpos = qstart + np.cumsum(~qgap) - 1          # 1-based genomic pos at each column
        spos = sstart + step * (np.cumsum(~sgap) - 1)  # consensus pos at each column
        if want_cons:
            n_cons[name] = n_cons.get(name, 0) + 1
            bi = BASE_IDX[np.frombuffer(sseq.encode(), dtype=np.uint8)]
            ok = ~sgap & (bi >= 0) & (spos >= 1)
            if comp:  # consensus line of a complement alignment shows the complemented consensus
                bi = np.where(bi >= 0, COMPL_IDX[np.clip(bi, 0, 3)], -1)
            cc = cons_counts.setdefault(name, np.zeros((1, 4), dtype=np.int64))
            mx = int(spos[ok].max()) if ok.any() else 0
            if mx >= len(cc):
                cc = np.vstack([cc, np.zeros((mx + 1 - len(cc), 4), dtype=np.int64)])
                cons_counts[name] = cc
            np.add.at(cc, (spos[ok], bi[ok]), 1)
        if cid is None:
            continue
        seen_copies.add(cid)
        both = ~qgap & ~sgap & (qpos - 1 >= cstart) & (qpos - 1 < cend)
        out_c.append(np.full(both.sum(), cid, dtype=object))
        out_p.append(qpos[both] - 1)
        out_r.append(spos[both])
    n_seen = len(seen_copies)
    log(f"  matched alignments: {n_seen} / {len(cp)} copies ({n_seen / max(len(cp), 1):.1%})")
    if not out_c:
        return pd.DataFrame(columns=["copy_id", "pos", "cons_pos"])
    res = pd.DataFrame({"copy_id": np.concatenate(out_c), "pos": np.concatenate(out_p).astype(np.int64),
                        "cons_pos": np.concatenate(out_r).astype(np.int32)})
    dup = res.duplicated(["copy_id", "pos"])
    if dup.any():
        log(f"  {int(dup.sum())} bases covered by >1 alignment of the same copy: first kept")
        res = res[~dup]
    rate = cp.assign(m=cp["copy_id"].isin(seen_copies)).groupby("subfamily")["m"].mean()
    low = rate[rate < 0.9]
    if len(low):
        log("  subfamilies with <90% matched copies: " +
            ", ".join(f"{k} {v:.0%}" for k, v in low.sort_values().items()))
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--copies", required=True, help="te_copies.tsv.gz from step 2")
    ap.add_argument("--mode", choices=["linear", "align"], default="linear")
    ap.add_argument("--align-file", default=None)
    ap.add_argument("--out", required=True, help="output .parquet")
    ap.add_argument("--compare-to", default=None, help="another mapping parquet to compare with")
    ap.add_argument("--consensus-out", default=None, help="align mode: write reconstructed consensus FASTA")
    ap.add_argument("--consensus-max-aln", type=int, default=3000,
                    help="alignments per subfamily used for consensus reconstruction")
    ap.add_argument("--split-model", default="3end", help="which split LINE model to use (5end/orf2/3end)")
    a = ap.parse_args()

    cp = pd.read_csv(a.copies, sep="\t", dtype={"chrom": str})
    if a.mode == "linear":
        res = map_linear(cp)
    else:
        if not a.align_file:
            ap.error("--align-file is required for --mode align")
        counts = {} if a.consensus_out else None
        res = map_align(cp, a.align_file, counts, a.consensus_max_aln, a.split_model)
        if a.consensus_out:
            with open(a.consensus_out, "w") as fh:
                for name, cc in sorted(counts.items()):
                    seq = "".join("ACGT"[r.argmax()] if r.sum() else "N" for r in cc[1:])
                    fh.write(f">{name}\n{seq}\n")
            log(f"  wrote reconstructed consensus for {len(counts)} subfamilies -> {a.consensus_out}")
    bad = (res["cons_pos"] < 1).sum()
    if bad:
        log(f"  WARNING: {bad} bases with cons_pos < 1 removed")
        res = res[res["cons_pos"] >= 1]
    res.to_parquet(a.out, index=False)
    log(f"  wrote {len(res)} mapped bases from {res['copy_id'].nunique()} copies -> {a.out}")

    if a.compare_to:
        other = pd.read_parquet(a.compare_to)
        m = res.merge(other, on=["copy_id", "pos"], suffixes=("", "_other"))
        d = (m["cons_pos"] - m["cons_pos_other"]).abs()
        log(f"  comparison on {len(m)} shared bases: identical {np.mean(d == 0):.1%}, "
            f"|diff|<=2: {np.mean(d <= 2):.1%}, |diff|<=5: {np.mean(d <= 5):.1%}, median |diff| {d.median()}")


if __name__ == "__main__":
    main()
