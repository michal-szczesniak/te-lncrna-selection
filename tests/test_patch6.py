#!/usr/bin/env python3
"""Tests of the patch-6 steps (41 canonical exons, 42 starts in LTRs, 43 facts, 37 --well-supported) and of the
figure and table scripts on small synthetic inputs. No network, about a minute.

Usage: python tests/test_patch6.py <scripts dir> <scratch dir> [<manuscript dir with figures/ and tables/>]
Ends with 'ALL PATCH 6 CHECKS PASSED'.
"""
import gzip
import os
import subprocess
import sys

import numpy as np
import pandas as pd

S, TMP = sys.argv[1], sys.argv[2]
M = sys.argv[3] if len(sys.argv) > 3 else None
os.makedirs(TMP, exist_ok=True)
FAIL = []
rng = np.random.default_rng(6)


def check(name, ok, detail=""):
    print(f"[{'ok' if ok else 'FAIL'}] {name} {detail}")
    if not ok:
        FAIL.append(name)


def run(args):
    env = dict(os.environ, PYTHONPATH=S + os.pathsep + os.environ.get("PYTHONPATH", ""), MPLBACKEND="Agg")
    r = subprocess.run([sys.executable] + args, capture_output=True, text=True, env=env)
    if r.returncode != 0:
        print(r.stdout[-3000:], r.stderr[-3000:])
    return r.returncode == 0


def gtf_line(chrom, feat, s, e, strand, gid, tid=None, gtype="lncRNA", canonical=False):
    attr = f'gene_id "{gid}"; gene_type "{gtype}"; gene_name "{gid}";'
    if tid:
        attr += f' transcript_id "{tid}";'
    if canonical:
        attr += ' tag "Ensembl_canonical";'
    return "\t".join([chrom, "T", feat, str(s), str(e), ".", strand, ".", attr])


# ------------------------------------------------------------------ synthetic GTF
# G1 (+): canonical T1 exons 1001-1200, 3001-3200; alternative T2 adds exon 2001-2100 and starts at 1101
# G2 (-): canonical T3 exons 10001-10300; T4 starts at 10251 (another start)
lines = [gtf_line("chr1", "gene", 1001, 3200, "+", "G1"),
         gtf_line("chr1", "transcript", 1001, 3200, "+", "G1", "T1", canonical=True),
         gtf_line("chr1", "exon", 1001, 1200, "+", "G1", "T1"), gtf_line("chr1", "exon", 3001, 3200, "+", "G1", "T1"),
         gtf_line("chr1", "transcript", 1101, 3200, "+", "G1", "T2"),
         gtf_line("chr1", "exon", 1101, 1200, "+", "G1", "T2"), gtf_line("chr1", "exon", 2001, 2100, "+", "G1", "T2"),
         gtf_line("chr1", "exon", 3001, 3200, "+", "G1", "T2"),
         gtf_line("chr1", "gene", 10001, 10300, "-", "G2"),
         gtf_line("chr1", "transcript", 10001, 10300, "-", "G2", "T3", canonical=True),
         gtf_line("chr1", "exon", 10001, 10300, "-", "G2", "T3"),
         gtf_line("chr1", "transcript", 10001, 10251, "-", "G2", "T4"),
         gtf_line("chr1", "exon", 10001, 10251, "-", "G2", "T4")]
with gzip.open(f"{TMP}/test.gtf.gz", "wt") as fh:
    fh.write("##test\n" + "\n".join(lines) + "\n")

# ------------------------------------------------------------------ 41 canonical exons
os.makedirs(f"{TMP}/annot", exist_ok=True)
rows = []
for gid, comp, lo, hi in (("G1", "exon", 1000, 1200), ("G1", "exon", 2000, 2100), ("G1", "intron", 1500, 1600),
                          ("G1", "exon", 3000, 3200), ("G2", "exon", 10000, 10300)):
    for p in range(lo, hi):
        rows.append({"chrom": "chr1", "pos": p, "compartment": comp, "gene_id": gid,
                     "unit_id": f"{gid}_{comp}_{lo}", "copy_id": f"c{lo}"})
pd.DataFrame(rows).to_parquet(f"{TMP}/annot/part-chr1.parquet", index=False)
ok = run([f"{S}/41_canonical_filter.py", "--annot", f"{TMP}/annot", "--gtf", f"{TMP}/test.gtf.gz",
          "--out-dir", f"{TMP}/annot_canon"])
check("41 runs", ok)
if ok:
    d = pd.read_parquet(f"{TMP}/annot_canon/part-chr1.parquet")
    ex = d[d.compartment == "exon"]
    check("41 drops exonic bases of non-canonical exons only",
          ex.pos.between(2000, 2099).sum() == 0 and len(ex) == 200 + 200 + 300, f"{len(ex)} exonic bases")
    check("41 keeps intronic bases", (d.compartment == "intron").sum() == 100)
    s = pd.read_csv(f"{TMP}/annot_canon/canonical_summary.tsv", sep="\t").set_index("gene_id")
    check("41 summary counts", s.loc["G1", "bases_before"] == 500 and s.loc["G1", "bases_after"] == 400,
          s.to_dict())

# ------------------------------------------------------------------ 42 starts in LTRs
# 300 lncRNA genes on chr2, each with one transcript; starts sit at consensus position ~60 of a sense LTR12C copy
# (copy covers consensus 1-400) in 240 genes, and at random positions of antisense THE1B copies in 60 genes
glines, rm, lg = [], [], []
for i in range(300):
    base = 100000 + i * 5000
    sense = i < 240
    strand = "+" if i % 2 == 0 else "-"
    cs, ce = 1, 400
    copy_s, copy_e = base, base + 400                      # 0-based [copy_s, copy_e)
    if sense:
        te_strand = strand
        if strand == "+":
            tss0 = copy_s + 59 + int(rng.integers(-3, 4))
        else:
            tss0 = copy_e - 1 - (59 + int(rng.integers(-3, 4)))
        name, cls = "LTR12C", "LTR/ERV1"
    else:
        te_strand = "-" if strand == "+" else "+"
        tss0 = copy_s + int(rng.integers(0, 400))
        name, cls = "THE1B", "LTR/ERVL-MaLR"
    if strand == "+":
        ts, te_ = tss0 + 1, tss0 + 2000
    else:
        ts, te_ = tss0 - 1999 + 1, tss0 + 1
    gid = f"L{i}"
    glines += [gtf_line("chr2", "gene", ts, te_, strand, gid), gtf_line("chr2", "transcript", ts, te_, strand, gid,
                                                                        f"LT{i}", canonical=True),
               gtf_line("chr2", "exon", ts, te_, strand, gid, f"LT{i}")]
    lg.append({"gene_id": gid, "chrom": "chr2", "drop_reason": ""})
    if te_strand == "+":
        rm.append(f"1000 10.0 0.0 0.0 chr2 {copy_s + 1} {copy_e} (1000) + {name} {cls} {cs} {ce} (0) {i + 1}")
    else:
        rm.append(f"1000 10.0 0.0 0.0 chr2 {copy_s + 1} {copy_e} (1000) C {name} {cls} (0) {ce} {cs} {i + 1}")
with gzip.open(f"{TMP}/ltr.gtf.gz", "wt") as fh:
    fh.write("\n".join(glines) + "\n")
with gzip.open(f"{TMP}/rmsk.out.gz", "wt") as fh:
    fh.write("   SW  perc perc perc  query position\nscore div. del. ins. sequence\n\n" + "\n".join(rm) + "\n")
pd.DataFrame(lg).to_csv(f"{TMP}/lnc_genes.tsv.gz", sep="\t", index=False)
ok = run([f"{S}/42_tss_ltr.py", "--gtf", f"{TMP}/ltr.gtf.gz", "--lnc-genes", f"{TMP}/lnc_genes.tsv.gz",
          "--rmsk-out", f"{TMP}/rmsk.out.gz", "--chroms", "chr2", "--n-perm", "200", "--out", f"{TMP}/tss"])
check("42 runs", ok)
if ok:
    s = pd.read_csv(f"{TMP}/tss_summary.tsv", sep="\t").set_index("group")
    t = pd.read_csv(f"{TMP}/tss_tss.tsv", sep="\t")
    check("42 finds every start in an LTR", len(t) == 300, len(t))
    check("42 orientation of planted promoters", s.loc["LTR12", "share_sense"] == 1.0
          and s.loc["THE1", "share_sense"] == 0.0, s["share_sense"].to_dict())
    ws = s.loc["LTR12", "window_start"]
    check("42 recurrent position recovered (consensus ~60)", 10 <= ws <= 60 and s.loc["LTR12", "window_share"] > 0.95
          and s.loc["LTR12", "window_p"] < 0.01,
          f"window {ws}-{s.loc['LTR12', 'window_end']}, share {s.loc['LTR12', 'window_share']:.2f}, "
          f"P {s.loc['LTR12', 'window_p']:.3f}")
    c = t[t.group == "LTR12"].cons_pos
    check("42 consensus positions in both strands", c.between(55, 65).all(), c.describe().to_dict())

# ------------------------------------------------------------------ 43 facts
os.makedirs(f"{TMP}/step1", exist_ok=True)
pd.DataFrame({"copy_id": ["a", "a", "b", "c", "d"], "unit_id": ["a1", "a2", "b1", "c1", "d1"],
              "chrom": ["chr1", "chr1", "chr2", "chrX", "chr3"], "subfamily": ["AluSx", "AluSx", "MIRb", "L2a", "THE1B"],
              "compartment": ["exon", "intron", "intron", "exon", "exon"]}).to_csv(
    f"{TMP}/te_copies.tsv.gz", sep="\t", index=False)
pd.DataFrame({"gene_id": ["G1", "G2", "G3"], "chrom": ["chr1", "chr1", "chrX"],
              "drop_reason": ["", "low_expression", ""]}).to_csv(f"{TMP}/step1/lnc_genes.tsv.gz", sep="\t", index=False)
pd.DataFrame({"chrom": ["chr1", "chr1"], "start": [0, 500], "end": [100, 600], "compartment": ["exon", "intron"],
              "gene_id": ["G1", "G1"]}).to_csv(f"{TMP}/step1/compartments.tsv.gz", sep="\t", index=False)
ok = run([f"{S}/43_facts.py", "--copies", f"{TMP}/te_copies.tsv.gz", "--step1", f"{TMP}/step1", "--out", f"{TMP}/facts"])
check("43 runs", ok)
if ok:
    f = pd.read_csv(f"{TMP}/facts.tsv", sep="\t").set_index("quantity")["value"]
    check("43 counts autosomal Alu/MIR/L2 copies once", int(f["Alu/MIR/L2 copies, autosomes, all compartments"]) == 2,
          f.to_dict())
    check("43 kept genes", int(f["kept lncRNA genes (autosomes)"]) == 1)

# ------------------------------------------------------------------ 37 --well-supported
gl = []
for i in range(300):
    for sf in ("AluS*", "MIR*"):
        for ori in ("sense", "antisense", "pooled"):
            gl.append({"subfamily": sf, "orientation": ori, "gene_id": f"ENSG{i:05d}.3", "n_exon": 50,
                       "n_intron": 200, "gtex": float(i) / 30 if i < 200 else np.nan, "tau": 1 - i / 300,
                       "excess_poly_noncpg": rng.normal(0, 0.01), "diff_poly_noncpg": rng.normal(0, 0.01)})
pd.DataFrame(gl).to_csv(f"{TMP}/genes.tsv", sep="\t", index=False)
pd.DataFrame({"gene_id": [f"ENSG{i:05d}" for i in range(300)],
              "project_class": ["lincRNA" if i % 3 == 0 else "antisense" for i in range(300)]}).to_csv(
    f"{TMP}/gene_class.tsv", sep="\t", index=False)
os.makedirs(f"{TMP}/comp", exist_ok=True)
for g in ("AluS", "MIR"):
    pd.DataFrame({"poly_c_noncpg_intron": [0.2, 0.2], "units_intron": [10, 10]}).to_csv(
        f"{TMP}/comp/te_{g}_paired_bins.tsv", sep="\t", index=False)
ok = run([f"{S}/37_bounds_subsets.py", "--genes", f"{TMP}/genes.tsv", "--comp-dir", f"{TMP}/comp", "--n-boot", "200",
          "--well-supported", "--gene-class", f"{TMP}/gene_class.tsv", "--skip-expression-subsets",
          "--out", f"{TMP}/bwell"])
check("37 --well-supported runs", ok)
if ok:
    sub = pd.read_csv(f"{TMP}/bwell_subsets.tsv", sep="\t").set_index("subset")["genes"]
    check("37 well-supported subsets", sub.get("gtex_any") == 200 and sub.get("lincrna_v29") == 100
          and sub.get("lincrna_v29_gtex") == 67 and "expr_top25" not in sub, sub.to_dict())
    b = pd.read_csv(f"{TMP}/bwell_bounds.tsv", sep="\t")
    check("37 bounds for every subset", set(b["table"]) == {"all", "gtex_any", "lincrna_v29", "lincrna_v29_gtex"})

# ------------------------------------------------------------------ figure S8 from the step-42 output
if M and not FAIL:
    os.makedirs(f"{TMP}/rev6", exist_ok=True)
    for f in ("summary", "hist"):
        os.replace(f"{TMP}/tss_{f}.tsv", f"{TMP}/rev6/tss_ltr_{f}.tsv")
    ok = run([os.path.join(M, "figures", "make_figS8.py"), f"{TMP}/rev6", f"{TMP}/FigS8"])
    check("Fig. S8 is drawn", ok and os.path.getsize(f"{TMP}/FigS8.png") > 10000)

print("ALL PATCH 6 CHECKS PASSED" if not FAIL else f"FAILED: {FAIL}")
sys.exit(1 if FAIL else 0)
