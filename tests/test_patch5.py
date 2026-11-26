#!/usr/bin/env python3
"""Tests of the patch-5 scripts (steps 33-38) on small synthetic inputs. No network, a few minutes.

Usage: python tests/test_patch5.py <scripts dir> <scratch dir>
The scripts directory must hold the frozen pipeline scripts (17_te_architecture.py, ridl_utils.py, ridl_stats.py)
and the revision scripts (23_bounds.py) next to the patch-5 scripts. Ends with 'ALL PATCH 5 CHECKS PASSED'.
"""
import gzip
import importlib.util
import os
import subprocess
import sys

import numpy as np
import pandas as pd

S, TMP = sys.argv[1], sys.argv[2]
os.makedirs(TMP, exist_ok=True)
sys.path.insert(0, S)
FAIL = []


def check(name, ok, detail=""):
    print(f"[{'ok' if ok else 'FAIL'}] {name} {detail}")
    if not ok:
        FAIL.append(name)


def load(fname, modname):
    spec = importlib.util.spec_from_file_location(modname, os.path.join(S, fname))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def run(args):
    r = subprocess.run([sys.executable] + args, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-3000:], r.stderr[-3000:])
    return r.returncode == 0


rng = np.random.default_rng(7)

# ------------------------------------------------------------------ 33: weights and planted signal
m33 = load("33_spikein_scan.py", "scan33")
rows = []
for u in range(400):
    lab = int(u < 120)
    gene = f"g{u % 60}"
    for b in range(10):
        for k in range(15):
            rows.append((f"u{u}", gene, lab, b, True, float(rng.random() < 0.2), rng.integers(0, 4), "sf"))
d = pd.DataFrame(rows, columns=["unit_id", "gene_id", "_lab", "bin", "gn_callable", "gn_poly", "tri", "subfamily_orig"])
d["cpg"] = 0
comp = m33.Comparison(d, 10, ["tri", "cpg", "subfamily_orig"], 50, "global", 5, np.random.default_rng(1))
r0 = comp.run()
r1 = comp.run(None, np.ones(len(comp.agg.units)))
check("33 unit weights of one reproduce the unweighted statistics",
      np.allclose(r0["bin"], r1["bin"], equal_nan=True) and np.isclose(r0["all"], r1["all"]))
mult = {g: 2.0 for g in comp.agg.clusters}
w2 = comp.unit_weights(mult)
r2 = comp.run(None, w2)
check("33 doubling every gene leaves the log ORs unchanged", np.allclose(r0["bin"], r2["bin"], equal_nan=True))
v, n = comp.thinned(np.array([3, 4]), 0.9, 1.0, np.random.default_rng(2))
r3 = comp.run(v)
check("33 a strong planted deficit is negative in the planted bins and near zero elsewhere",
      np.all(r3["bin"][[3, 4]] < -0.5) and np.all(np.abs(r3["bin"][[0, 1, 7, 8]]) < 0.4),
      f"planted {np.round(r3['bin'][[3, 4]], 2)}")

# 33 end to end: scan mode with gene resampling on a small annotation
rows = []
for u in range(300):
    comp_ = "exon" if u < 90 else "intron"
    for ori in ("sense", "antisense"):
        for pos_ in range(1, 201):
            rows.append({"unit_id": f"{ori}{u}", "subfamily": "AluS*", "compartment": comp_, "orientation": ori,
                         "gene_id": f"g{u % 45}", "cons_pos": pos_, "gn_callable": True,
                         "gn_poly": bool(rng.random() < 0.2), "tri": "ACA", "cpg": False, "subfamily_orig": "AluSx",
                         "chrom": "chr1"})
os.makedirs(f"{TMP}/annot33", exist_ok=True)
pd.DataFrame(rows).to_parquet(f"{TMP}/annot33/chr1.parquet", index=False)
ok = run([f"{S}/33_spikein_scan.py", "--annot", f"{TMP}/annot33", "--subfamilies", "AluS*", "--orientation", "sense",
          "--scan", "40", "--scan-step", "2", "--reductions", "0,0.5", "--fractions", "1", "--reps", "3",
          "--n-perm", "300", "--boot-genes", "--min-units", "10", "--out", f"{TMP}/scan33"])
if ok:
    P = pd.read_csv(f"{TMP}/scan33_power.tsv", sep="\t")
    check("33 scan covers the element in steps of two bins", sorted(P["seg_start"].unique()) == [1, 41, 81, 121, 161],
          str(sorted(P["seg_start"].unique())))
    strong = P[(P.reduction == 0.5)]["power_window"].mean()
    null = P[(P.reduction == 0)]["power_window"].mean()
    check("33 planted 50% reduction in all exonic copies is detected far more often than none", strong > 0.6 and null < 0.3,
          f"r=0.5: {strong:.2f}, r=0: {null:.2f}")
else:
    check("33 scan runs", False)

# ------------------------------------------------------------------ 34: usage summary
pos = []
us = []
for st, n in (("pc_site", 90), ("lnc_TE_site", 60), ("lnc_nonTE_site", 60)):
    for i in range(n):
        key = f"chr1:{1000 + 100 * i + hash(st) % 50}-{5000 + 100 * i}+"
        pos.append({"set": st, "site_id": key + "d"})
        pos.append({"set": st, "site_id": key + "a"})
        matched = i % 3 != 0
        reads = float(i * 10) if matched and i % 5 else 0.0
        us.append({"intron_key": key, "matched": matched, "total_reads": reads, "mean_reads": reads / 100,
                   "samples_ge1": int(reads > 0) * 5, "tissues_mean_ge1": int(reads > 0)})
pd.DataFrame(pos).to_csv(f"{TMP}/pos.tsv.gz", sep="\t", index=False)
pd.DataFrame(us).drop_duplicates("intron_key").to_csv(f"{TMP}/usage.tsv.gz", sep="\t", index=False)
ok = run([f"{S}/34_usage_summary.py", "--positions", f"{TMP}/pos.tsv.gz", "--usage", f"{TMP}/usage.tsv.gz",
          "--out", f"{TMP}/us"])
if ok:
    st = pd.read_csv(f"{TMP}/us_sets.tsv", sep="\t").set_index("set")
    tt = pd.read_csv(f"{TMP}/us_tertiles.tsv", sep="\t")
    pc = st.loc["pc_site"]
    check("34 share in the GTEx table", abs(pc["in_gtex_file"] - 60 / 90) < 1e-9, f"{pc['in_gtex_file']:.3f}")
    t3 = tt[(tt.set == "pc_site")]
    check("34 tertiles ordered by reads", list(t3["tertile"]) == ["usage q1", "usage q2", "usage q3"]
          and (t3["reads_min"].diff().dropna() > 0).all())
else:
    check("34 runs", False)

# ------------------------------------------------------------------ 35: exact null vs simulation
m35 = load("35_boundary_null.py", "bn35")
good = True
for _ in range(200):
    e0 = int(rng.integers(1000, 2000)); e1 = e0 + int(rng.integers(30, 400)); L = int(rng.integers(20, 400))
    sites = np.sort(np.unique(np.concatenate([[e0, e1 - 1], rng.integers(e0 - 3000, e1 + 3000, 3)])))
    p = m35.null_prob(L, e0, e1, sites, 10, 10)
    c0, c1 = e0 + 10, e1 - 10
    if c1 <= c0:
        c0, c1 = e0, e1
    xs = np.arange(c0 - L + 1, c1)
    sim = np.mean([np.any((sites >= x - 10) & (sites < x + L + 10)) for x in xs])
    good &= abs(p - sim) < 1e-9
check("35 exact null equals enumeration of all start positions", bool(good))

# 35 end to end on the synthetic gene of test 38 is done below; here a minimal unit table
os.makedirs(f"{TMP}/annot35", exist_ok=True)
pd.DataFrame({"unit_id": ["u1"] * 50 + ["u2"] * 50, "copy_id": ["c1"] * 50 + ["c2"] * 50, "compartment": "exon",
              "gene_id": "G1", "orientation": ["antisense"] * 50 + ["sense"] * 50,
              "unit_class": ["boundary"] * 50 + ["embedded"] * 50, "chrom": "chr1",
              "pos": list(range(1150, 1200)) + list(range(3100, 3150))}).to_parquet(f"{TMP}/annot35/chr1.parquet", index=False)
pd.DataFrame({"copy_id": ["c1", "c2"], "chrom": ["chr1", "chr1"], "start": [1150, 3100], "end": [1250, 3150]}).to_csv(
    f"{TMP}/copies35.tsv.gz", sep="\t", index=False)
pd.DataFrame({"chrom": ["chr1", "chr1"], "pos": [1199, 3000], "gene_id": ["G1", "G1"]}).to_csv(
    f"{TMP}/ss35.tsv.gz", sep="\t", index=False)

# ------------------------------------------------------------------ 36: copy-aware ages
sp = ["panTro6", "rheMac10", "mm39"]
S36 = pd.DataFrame({
    "set": ["lnc_TE_site"] * 4 + ["pc_site"],
    "site_id": list("abcde"), "rep_family": ["Alu", "Alu", "MIR", "MIR", ""],
    "te_orientation": ["antisense"] * 4 + [""], "ancestral": ["created", "preformed", "preformed", "created", ""],
    "site_type": ["acceptor"] * 5,
    "site_panTro6": ["present"] * 5, "copy_panTro6": ["present"] * 4 + [""],
    "site_rheMac10": ["present", "present", "present", "absent", "present"], "copy_rheMac10": ["present"] * 4 + [""],
    "site_mm39": ["present", "unaligned", "present", "present", "present"],
    "copy_mm39": ["absent", "absent", "present", "absent", ""],
    "age_class": ["shared with mm39", "shared with rheMac10", "shared with mm39", "shared with mm39", "shared with mm39"]})
S36.to_csv(f"{TMP}/age_sites.tsv.gz", sep="\t", index=False)
ok = run([f"{S}/36_age_copyaware.py", "--sites", f"{TMP}/age_sites.tsv.gz", "--species", ",".join(sp),
          "--out", f"{TMP}/ageca"])
if ok:
    o = pd.read_csv(f"{TMP}/ageca_sites.tsv.gz", sep="\t").set_index("site_id")["age_class_ca"]
    check("36 Alu site present in mouse without the copy is dated to macaque", o["a"] == "shared with rheMac10", o["a"])
    check("36 MIR site with the copy in mouse keeps its mouse age", o["c"] == "shared with mm39", o["c"])
    check("36 MIR site without the copy in mouse falls back to chimpanzee", o["d"] == "shared with panTro6", o["d"])
    check("36 non-TE sites unchanged", o["e"] == "shared with mm39", o["e"])
else:
    check("36 runs", False)

# ------------------------------------------------------------------ 37: subsets -> 23_bounds
gl = []
for i in range(300):
    for sf in ("AluS*", "MIR*"):
        for ori in ("sense", "antisense", "pooled"):
            gl.append({"subfamily": sf, "orientation": ori, "gene_id": f"G{i}", "n_exon": 50, "n_intron": 200,
                       "gtex": float(i) / 30, "tau": 1 - i / 300,
                       "excess_poly_noncpg": rng.normal(0, 0.01), "diff_poly_noncpg": rng.normal(0, 0.01)})
pd.DataFrame(gl).to_csv(f"{TMP}/genes.tsv", sep="\t", index=False)
os.makedirs(f"{TMP}/comp", exist_ok=True)
for g in ("AluS", "MIR"):
    pd.DataFrame({"poly_c_noncpg_intron": [0.2, 0.2], "units_intron": [10, 10]}).to_csv(
        f"{TMP}/comp/te_{g}_paired_bins.tsv", sep="\t", index=False)
ok = run([f"{S}/37_bounds_subsets.py", "--genes", f"{TMP}/genes.tsv", "--comp-dir", f"{TMP}/comp", "--n-boot", "200",
          "--out", f"{TMP}/bsub"])
if ok:
    b = pd.read_csv(f"{TMP}/bsub_bounds.tsv", sep="\t")
    sub = pd.read_csv(f"{TMP}/bsub_subsets.tsv", sep="\t").set_index("subset")["genes"]
    check("37 subsets have the expected gene numbers", sub["expr_top25"] == 75 and sub["expr_top10"] == 30
          and sub["broad_tau25"] == 75, sub.to_dict())
    check("37 bounds written for every subset", set(b["table"]) == {"all", "expr_top25", "expr_top10", "broad_tau25"})
else:
    check("37 runs", False)

# ------------------------------------------------------------------ 38: local background
# one gene on chr1 (+): exons [1000,1200) and [3000,3300); an Alu covering [990,1100) and an L2 at [2900,3100)
gtf = [("chr1", "t", "gene", 1001, 3300, ".", "+", ".", 'gene_id "G1"; gene_type "lncRNA";'),
       ("chr1", "t", "transcript", 1001, 3300, ".", "+", ".", 'gene_id "G1"; transcript_id "T1"; gene_type "lncRNA";'),
       ("chr1", "t", "exon", 1001, 1200, ".", "+", ".", 'gene_id "G1"; transcript_id "T1"; gene_type "lncRNA";'),
       ("chr1", "t", "exon", 3001, 3300, ".", "+", ".", 'gene_id "G1"; transcript_id "T1"; gene_type "lncRNA";')]
with gzip.open(f"{TMP}/t.gtf.gz", "wt") as fh:
    for r in gtf:
        fh.write("\t".join(map(str, r)) + "\n")
rmsk = ["   SW  perc perc perc  query      position in query           matching       repeat              position in  repeat",
        "score  div. del. ins.  sequence    begin     end    (left)    repeat         class/family         begin  end (left)   ID",
        "",
        " 1000  10.0  1.0  1.0  chr1          991     1100 (1000)   +  AluSx          SINE/Alu                1    110  (190)  1",
        " 1000  10.0  1.0  1.0  chr1         2901     3100 (1000)   C  L2a            LINE/L2              (10)   3000   2801  2"]
with open(f"{TMP}/t.out", "w") as fh:
    fh.write("\n".join(rmsk) + "\n")
pd.DataFrame({"gene_id": ["G1"], "drop_reason": [""]}).to_csv(f"{TMP}/lnc.tsv.gz", sep="\t", index=False)
ok = run([f"{S}/38_feature_local_bg.py", "--gtf", f"{TMP}/t.gtf.gz", "--gene-set", f"lnc={TMP}/lnc.tsv.gz",
          "--rmsk-out", f"{TMP}/t.out", "--chroms", "chr1", "--dist-bins", "20,50,100", "--per-feature", "20",
          "--out", f"{TMP}/fb"])
if ok:
    L = pd.read_csv(f"{TMP}/fb_local.tsv", sep="\t")
    f = L[L.side == "feature"].set_index("feature")["any TE"]
    check("38 features inside TEs (TSS in Alu, acceptor in L2)", f["TSS"] == 1 and f["acceptor"] == 1, f.to_dict())
    tss_ex = L[(L.feature == "TSS") & (L.side == "exon side") & (L.dist_bin == "50-100")]
    check("38 TSS controls 50-100 nt downstream fall in the exon and in the Alu",
          set(tss_ex["falls_in"]) == {"exon"} and float(tss_ex["Alu"].iloc[0]) == 1.0)
    don_in = L[(L.feature == "donor") & (L.side == "other side")]
    check("38 donor controls on the other side fall in the intron", set(don_in["falls_in"]) == {"intron"})
else:
    check("38 runs", False)

ok = run([f"{S}/35_boundary_null.py", "--annot-sub", f"AluS={TMP}/annot35", "--copies", f"{TMP}/copies35.tsv.gz",
          "--splice-sites", f"{TMP}/ss35.tsv.gz", "--gtf", f"{TMP}/t.gtf.gz", "--n-boot", "50", "--out", f"{TMP}/bn"])
if ok:
    b = pd.read_csv(f"{TMP}/bn_units.tsv.gz", sep="\t").set_index("unit_id")
    # c1: 100 nt copy on the 200-nt first exon (sites 1199 and 3000): exact probability by enumeration
    xs = np.arange(1010 - 100 + 1, 1190)
    p1 = np.mean([np.any((np.array([1199, 3000]) >= x - 10) & (np.array([1199, 3000]) < x + 110)) for x in xs])
    check("35 end to end: null of a copy on the first exon", abs(b.loc["u1", "p_null"] - p1) < 1e-9,
          f"{b.loc['u1', 'p_null']:.4f} vs {p1:.4f}")
    check("35 observed classes reproduced from coordinates", (b["boundary_obs"] == b["boundary_recomputed"]).all())
else:
    check("35 runs", False)

print("ALL PATCH 5 CHECKS PASSED" if not FAIL else f"FAILED: {FAIL}")
sys.exit(1 if FAIL else 0)
