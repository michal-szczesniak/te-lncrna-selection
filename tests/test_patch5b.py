#!/usr/bin/env python3
"""Tests of the patch-5b changes (step 36 --rule, step 39) on small synthetic inputs. No network, under a minute.

Usage: python tests/test_patch5b.py <scripts dir> <scratch dir>
Ends with 'ALL PATCH 5B CHECKS PASSED'.
"""
import os
import subprocess
import sys

import numpy as np
import pandas as pd

S, TMP = sys.argv[1], sys.argv[2]
os.makedirs(TMP, exist_ok=True)
FAIL = []


def check(name, ok, detail=""):
    print(f"[{'ok' if ok else 'FAIL'}] {name} {detail}")
    if not ok:
        FAIL.append(name)


def run(args):
    env = dict(os.environ, PYTHONPATH=S + os.pathsep + os.environ.get("PYTHONPATH", ""))
    r = subprocess.run([sys.executable] + args, capture_output=True, text=True, env=env)
    if r.returncode != 0:
        print(r.stdout[-3000:], r.stderr[-3000:])
    return r.returncode == 0


rng = np.random.default_rng(11)

# ------------------------------------------------------------------ 36: rule absent vs strict
sp = ["panTro6", "rheMac10", "mm39"]
S36 = pd.DataFrame({
    "set": ["lnc_TE_site"] * 3 + ["pc_site"],
    "site_id": list("abcd"), "rep_family": ["Alu", "Alu", "MIR", ""],
    "te_orientation": ["antisense"] * 3 + [""], "ancestral": ["created", "preformed", "preformed", ""],
    "site_type": ["acceptor"] * 4,
    "site_panTro6": ["present"] * 4, "copy_panTro6": ["present", "short", "present", ""],
    "site_rheMac10": ["present"] * 4, "copy_rheMac10": ["present", "short", "present", ""],
    "site_mm39": ["present", "absent", "present", "present"], "copy_mm39": ["absent", "short", "short", ""],
    "age_class": ["shared with mm39", "shared with rheMac10", "shared with mm39", "shared with mm39"]})
S36.to_csv(f"{TMP}/age_sites.tsv.gz", sep="\t", index=False)
for rule in ("absent", "strict"):
    ok = run([f"{S}/36_age_copyaware.py", "--sites", f"{TMP}/age_sites.tsv.gz", "--species", ",".join(sp),
              "--rule", rule, "--out", f"{TMP}/ageca_{rule}"])
    check(f"36 runs with --rule {rule}", ok)
if not FAIL:
    a = pd.read_csv(f"{TMP}/ageca_absent_sites.tsv.gz", sep="\t").set_index("site_id")["age_class_ca"]
    s = pd.read_csv(f"{TMP}/ageca_strict_sites.tsv.gz", sep="\t").set_index("site_id")["age_class_ca"]
    check("36 absent: copy absent in mouse vetoes the mouse call", a["a"] == "shared with rheMac10", a["a"])
    check("36 absent: short copy keeps the site calls", a["b"] == "shared with rheMac10" and a["c"] == "shared with mm39",
          f"{a['b']} / {a['c']}")
    check("36 strict: short copy vetoes every call", s["b"] == "unresolved" and s["c"] == "shared with rheMac10",
          f"{s['b']} / {s['c']}")
    check("36 non-TE site unchanged under both rules", a["d"] == s["d"] == "shared with mm39")
    r = pd.read_csv(f"{TMP}/ageca_strict_reclassified.tsv", sep="\t")
    check("36 reclassified table splits absent and short copies",
          {"copy_absent", "copy_short"} <= set(r.columns) and r["copy_short"].sum() == 3 and r["copy_absent"].sum() == 1,
          r.to_dict("records"))

# ------------------------------------------------------------------ 39: pooled conserved stretches
def units(group, subs, p_exon, p_intron, n_genes=800, start=0):
    rows = []
    for gi in range(start, start + n_genes):
        for comp, p, k in (("exon", p_exon, 3), ("intron", p_intron, 6)):
            for j in range(k):
                L = int(rng.integers(40, 300)) if comp == "intron" else int(rng.integers(30, 160))
                hit = rng.random() < p * (0.5 + L / 300)       # longer units more often conserved; exons shorter
                rows.append({"subfamily": group, "subfamily_orig": subs[j % len(subs)], "compartment": comp,
                             "orientation": "sense", "unit_id": f"{group}_{gi}_{comp}_{j}", "gene_id": f"G{gi}",
                             "chrom": "chr1", "n_bases": L,
                             "best_phylop_p447": 2.0 if hit else 0.0, "best_phylop_prim": 0.0})
    return rows


te = pd.DataFrame(units("AluS*", ["AluSx", "AluSz"], 0.12, 0.06) + units("MIR*", ["MIRb"], 0.06, 0.06))
nt = pd.DataFrame(units("nonTE", ["nonTE"], 0.12, 0.06))
nt.loc[nt.index[:30], "compartment"] = "intergenic"
te.to_csv(f"{TMP}/te_units.tsv.gz", sep="\t", index=False)
nt.to_csv(f"{TMP}/nt_units.tsv.gz", sep="\t", index=False)
ok = run([f"{S}/39_conserved_pooled.py", "--variant", "v1", f"{TMP}/te_units.tsv.gz", f"{TMP}/nt_units.tsv.gz",
          "--groups", "AluS*,MIR*,L2abc", "--thresholds", "1.5", "--n-boot", "200", "--out", f"{TMP}/pooled"])
check("39 runs (a missing group is skipped)", ok)
if ok:
    R = pd.read_csv(f"{TMP}/pooled.tsv", sep="\t")
    R = R[R["track"] == "p447"].set_index(["group", "design"])
    al, mi, nn = R.loc[("AluS*", "all genes")], R.loc[("MIR*", "all genes")], R.loc[("nonTE", "all genes")]
    check("39 exonic excess detected for the planted group", al["ci_lo"] > 1.2 and al["p_cmh"] < 0.01,
          f"OR {al['mh_or']:.2f} [{al['ci_lo']:.2f},{al['ci_hi']:.2f}]")
    check("39 no excess for the null group", mi["ci_lo"] < 1 < mi["ci_hi"],
          f"OR {mi['mh_or']:.2f} [{mi['ci_lo']:.2f},{mi['ci_hi']:.2f}]")
    check("39 TE-specific ratio covers 1 when the baseline has the same excess",
          al["ratio_lo"] < 1 < al["ratio_hi"], f"{al['ratio_vs_nonte']:.2f} [{al['ratio_lo']:.2f},{al['ratio_hi']:.2f}]")
    check("39 TE-specific ratio below 1 for the null group", mi["ratio_hi"] < 1,
          f"{mi['ratio_vs_nonte']:.2f} [{mi['ratio_lo']:.2f},{mi['ratio_hi']:.2f}]")
    check("39 intergenic units ignored", nn["units_exon"] + nn["units_intron"] == len(nt) - 30)
    sg = R.loc[("AluS*", "same genes")]
    check("39 same-genes design keeps genes with both compartments", sg["genes"] == 800, sg["genes"])

# ------------------------------------------------------------------ 40: unit table without subfamily_orig
os.makedirs(f"{TMP}/annot40", exist_ok=True)
rows = []
for k, (comp, gene) in enumerate([("exon", "G1"), ("intron", "G1"), ("exon", "G2")]):
    for p in range(60):
        rows.append({"copy_id": f"w{k}", "pos": 1000 * k + p, "cons_pos": p + 1, "unit_id": f"w{k}", "subfamily": "nonTE",
                     "compartment": comp, "orientation": "sense", "gene_id": gene, "chrom": "chr1",
                     "phylop_p447": 2.0 if (k == 0 and 10 <= p < 30) else 0.0, "phylop_prim": 0.0})
pd.DataFrame(rows).to_parquet(f"{TMP}/annot40/part-chr1.parquet", index=False)
ok = run([f"{S}/40_conserved_units.py", "--annot", f"{TMP}/annot40", "--out", f"{TMP}/u40"])
check("40 runs on an annotation without subfamily_orig", ok)
if ok:
    u = pd.read_csv(f"{TMP}/u40_units.tsv.gz", sep="\t").set_index("unit_id")
    check("40 fills subfamily_orig and finds the conserved stretch",
          (u["subfamily_orig"] == "nonTE").all() and abs(u.loc["w0", "best_phylop_p447"] - 2.0) < 1e-6
          and abs(u.loc["w1", "best_phylop_p447"]) < 1e-6, u[["best_phylop_p447"]].to_dict())

print("ALL PATCH 5B CHECKS PASSED" if not FAIL else f"FAILED: {FAIL}")
sys.exit(1 if FAIL else 0)
