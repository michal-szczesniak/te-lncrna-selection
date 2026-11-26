#!/usr/bin/env python3
"""Generate a small synthetic data set in real-world formats to test the pilot scripts.

Planted truth:
  * AluSx copies in lncRNA EXONS, SENSE orientation: consensus 101-130 has fewer gnomAD variants
    (x0.5), more singletons, higher phyloP (+1.5). Nothing elsewhere / in antisense / in introns.
  * MIRb copies in EXONS, SENSE: same kind of signal at consensus 51-80.
  * AluSz = AluSx with a 5-nt deletion (46-50), a 3-nt insertion after 200 and 3% substitutions;
    its exonic sense copies carry the signal at AluSz 96-125 (= AluSx 101-130 after grouping).
  * one AluSx copy per multi-exon gene straddles an exon/intron boundary (no planted signal).
Writes: genome.fa(+fai), chrom.sizes, gencode.gtf.gz, gtex.gct.gz, sd_mask.bed, rmsk.txt.gz (UCSC),
        rmsk.out.gz + rmsk.align.gz (RepeatMasker), consensus.fa, phylop.bw, gnomad.vcf.gz(+tbi),
        gnomad_chr1/2.vcf.gz (per chromosome, full headers),
        an.tsv.gz(+tbi), controls.bed, truth_mapping.tsv.gz, test.maf
"""
import gzip
import os
import sys

import numpy as np
import pysam
import pyBigWig

OUT = sys.argv[1] if len(sys.argv) > 1 else "synthetic"
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(42)
BASES = np.array(list("ACGT"))
COMP = str.maketrans("ACGTN-", "TGCAN-")


def rc(s):
    return s.translate(COMP)[::-1]


def rand_seq(n, gc=0.41):
    p = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
    return "".join(rng.choice(BASES, size=n, p=p))


CONS = {"AluSx": ("SINE", "Alu", rand_seq(300, 0.55)), "MIRb": ("SINE", "MIR", rand_seq(260, 0.45))}
PLANT = {("AluSx", "exon", "sense"): (101, 130), ("MIRb", "exon", "sense"): (51, 80),
         ("AluSz", "exon", "sense"): (96, 125)}
_sx = CONS["AluSx"][2]
_sz = _sx[:45] + _sx[50:200] + "GAT" + _sx[200:]
_sz = "".join(b if rng.random() >= 0.03 else str(rng.choice([x for x in "ACGT" if x != b])) for b in _sz)
CONS["AluSz"] = ("SINE", "Alu", _sz)
CHROMS = {"chr1": 6_000_000, "chr2": 6_000_000}
genome = {c: list(rand_seq(n)) for c, n in CHROMS.items()}

# ------------------------------------------------------------------ gene models
genes, gtf, gtex_rows = [], [], []
gid_n = 0
for chrom, size in CHROMS.items():
    pos = 30_000
    while pos < size - 120_000:
        gid_n += 1
        strand = "+" if rng.random() < 0.5 else "-"
        mono = rng.random() < 0.2
        n_ex = 1 if mono else int(rng.integers(3, 5))
        exons, p = [], pos
        for k in range(n_ex):
            L = int(rng.integers(3000, 5000)) if mono else int(rng.integers(1500, 3000))
            exons.append((p, p + L))
            p += L
            if k < n_ex - 1:
                p += int(rng.integers(3000, 8000))
        gid = f"ENSG{gid_n:011d}.1"
        genes.append(dict(chrom=chrom, start=exons[0][0], end=exons[-1][1], strand=strand, gid=gid,
                          exons=exons, gtype="lncRNA"))
        pos = exons[-1][1] + int(rng.integers(25_000, 40_000))
# one protein-coding gene overlapping exon 1 of the 3rd lncRNA (that lncRNA must be dropped)
g3 = genes[2]
genes.append(dict(chrom=g3["chrom"], start=g3["exons"][0][0] + 200, end=g3["exons"][0][0] + 900,
                  strand="+" if g3["strand"] == "-" else "-", gid="ENSG99999999999.1",
                  exons=[(g3["exons"][0][0] + 200, g3["exons"][0][0] + 900)], gtype="protein_coding"))
for g in genes:
    attr_g = f'gene_id "{g["gid"]}"; gene_type "{g["gtype"]}"; gene_name "G{g["gid"][-6:-2]}";'
    tid = g["gid"].replace("ENSG", "ENST")
    attr_t = attr_g + f' transcript_id "{tid}"; tag "Ensembl_canonical";'
    gtf.append(f'{g["chrom"]}\tHAVANA\tgene\t{g["start"] + 1}\t{g["end"]}\t.\t{g["strand"]}\t.\t{attr_g}')
    gtf.append(f'{g["chrom"]}\tHAVANA\ttranscript\t{g["start"] + 1}\t{g["end"]}\t.\t{g["strand"]}\t.\t{attr_t}')
    for k, (s, e) in enumerate(g["exons"]):
        gtf.append(f'{g["chrom"]}\tHAVANA\texon\t{s + 1}\t{e}\t.\t{g["strand"]}\t.\t{attr_t} exon_number {k + 1};')
    tpm = 0.2 if g["gid"].endswith("05.1") else float(rng.uniform(1.5, 40))
    gtex_rows.append(f'{g["gid"]}\tG\t{tpm:.3f}\t{tpm / 2:.3f}')
with gzip.open(f"{OUT}/gencode.gtf.gz", "wt") as fh:
    fh.write("##description: synthetic\n" + "\n".join(gtf) + "\n")
with gzip.open(f"{OUT}/gtex.gct.gz", "wt") as fh:
    fh.write(f"#1.2\n{len(gtex_rows)}\t2\nName\tDescription\tTissueA\tTissueB\n" + "\n".join(gtex_rows) + "\n")

# ------------------------------------------------------------------ TE placement
occupied = {c: [] for c in CHROMS}


def free(chrom, s, e):
    return all(e <= a or s >= b for a, b in occupied[chrom])


def make_copy(name, cs, ce):
    cons = CONS[name][2]
    ca, qa, truth = [], [], []
    for cpos in range(cs, ce + 1):
        base = cons[cpos - 1]
        r = rng.random()
        if cs < cpos < ce and r < 0.01:  # deletion in copy
            ca.append(base); qa.append("-"); continue
        if rng.random() < 0.10:
            base_q = rng.choice([b for b in "ACGT" if b != base])
        else:
            base_q = base
        ca.append(base); qa.append(base_q); truth.append(cpos)
        if cpos < ce and rng.random() < 0.01:  # insertion in copy
            ca.append("-"); qa.append(rng.choice(BASES)); truth.append(None)
    return "".join(ca), "".join(qa), truth


copies = []


def place(chrom, lo, hi, name, comp, gene=None, full=True, at=None):
    L = len(CONS[name][2])
    cs, ce = (1, L) if full or rng.random() < 0.7 else (int(rng.integers(1, 60)), L)
    ca, qa, truth = make_copy(name, cs, ce)
    seq = qa.replace("-", "")
    if at is not None:
        s = at
        if not free(chrom, s - 5, s + len(seq) + 5):
            return
    else:
        for _ in range(30):
            if hi - lo - len(seq) <= 0:
                return
            s = int(rng.integers(lo, hi - len(seq)))
            if free(chrom, s - 5, s + len(seq) + 5):
                break
        else:
            return
    te_strand = "+" if rng.random() < 0.5 else "-"
    gseq = seq if te_strand == "+" else rc(seq)
    genome[chrom][s:s + len(seq)] = list(gseq)
    occupied[chrom].append((s, s + len(seq)))
    ori = "." if gene is None else ("sense" if te_strand == gene["strand"] else "antisense")
    copies.append(dict(chrom=chrom, start=s, end=s + len(seq), te_strand=te_strand, name=name, cs=cs, ce=ce,
                       L=L, ca=ca, qa=qa, truth=truth, comp=comp, ori=ori,
                       gid=None if gene is None else gene["gid"]))


for g in genes:
    if g["gtype"] != "lncRNA":
        continue
    tss = g["start"] if g["strand"] == "+" else g["end"]
    for k, (s, e) in enumerate(g["exons"]):
        lo, hi = s + 150, e - 150
        if g["strand"] == "+" and k == 0:
            lo = max(lo, tss + 700)
        if g["strand"] == "-" and k == len(g["exons"]) - 1:
            hi = min(hi, tss - 700)
        for name in ("AluSx", "MIRb", "AluSx", "AluSz"):
            place(g["chrom"], lo, hi, name, "exon", g)
    if len(g["exons"]) > 1:  # straddles the end of exon 1 (exonic part ~120 nt, rest intronic)
        place(g["chrom"], 0, 0, "AluSx", "straddle", g, at=g["exons"][0][1] - 120)
    for (s1, e1), (s2, e2) in zip(g["exons"][:-1], g["exons"][1:]):
        lo, hi = e1 + 400, s2 - 120
        for name in ("AluSx", "AluSx", "MIRb", "AluSz"):
            place(g["chrom"], lo, hi, name, "intron", g)
for chrom, size in CHROMS.items():
    spans = sorted((g["start"], g["end"]) for g in genes if g["chrom"] == chrom)
    gaps = [(a[1] + 12_000, b[0] - 12_000) for a, b in zip(spans[:-1], spans[1:]) if b[0] - a[1] > 30_000]
    for i in range(400):
        lo, hi = gaps[int(rng.integers(len(gaps)))]
        place(chrom, lo, hi, ("AluSx", "MIRb", "AluSz")[i % 3], "intergenic")
print(f"placed {len(copies)} copies", file=sys.stderr)

# ------------------------------------------------------------------ genome, sizes, consensus
with open(f"{OUT}/genome.fa", "w") as fh:
    for c, s in genome.items():
        seq = "".join(s)
        fh.write(f">{c}\n" + "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\n")
pysam.faidx(f"{OUT}/genome.fa")
with open(f"{OUT}/chrom.sizes", "w") as fh:
    for c, n in CHROMS.items():
        fh.write(f"{c}\t{n}\n")
with open(f"{OUT}/consensus.fa", "w") as fh:
    for n, (cl, fam, s) in CONS.items():
        fh.write(f">{n}#{cl}/{fam}\n{s}\n")
with open(f"{OUT}/sd_mask.bed", "w") as fh:  # mask one intergenic stretch
    fh.write("chr2\t5000000\t5100000\n")

# ------------------------------------------------------------------ RepeatMasker formats
copies.sort(key=lambda c: (c["chrom"], c["start"]))
ucsc, outl, alig, truth_rows = [], [], [], []
for i, c in enumerate(copies):
    cl, fam, _ = CONS[c["name"]]
    left = c["L"] - c["ce"]
    minus = c["te_strand"] == "-"
    rs, rl = (-left, c["cs"]) if minus else (c["cs"], -left)
    ucsc.append(f'0\t2000\t100\t5\t5\t{c["chrom"]}\t{c["start"]}\t{c["end"]}\t{c["te_strand"]}\t{c["name"]}'
                f'\t{cl}\t{fam}\t{rs}\t{c["ce"]}\t{rl}\t0')
    qleft = CHROMS[c["chrom"]] - c["end"]
    rep = f'({left}) {c["ce"]} {c["cs"]}' if minus else f'{c["cs"]} {c["ce"]} ({left})'
    outl.append(f'2000 10.0 1.0 1.0 {c["chrom"]} {c["start"] + 1} {c["end"]} ({qleft}) '
                f'{"C" if minus else "+"} {c["name"]} {cl}/{fam} {rep} {i + 1}')
    q, s = (rc(c["qa"]), rc(c["ca"])) if minus else (c["qa"], c["ca"])
    blk = [f'2000 10.00 1.00 1.00 {c["chrom"]} {c["start"] + 1} {c["end"]} ({qleft}) '
           f'{"C " if minus else ""}{c["name"]}#{cl}/{fam} {rep} m_b1s001i{i} {i + 1}', ""]
    qp, sp = c["start"] + 1, (c["ce"] if minus else c["cs"])
    step = -1 if minus else 1
    for k in range(0, len(q), 50):
        qs_, ss_ = q[k:k + 50], s[k:k + 50]
        nq, ns = len(qs_.replace("-", "")), len(ss_.replace("-", ""))
        q_end, s_end = qp + nq - 1, sp + step * (ns - 1)
        blk.append(f'  {c["chrom"]:<14}{qp:>10} {qs_} {q_end}')
        blk.append(" " * 25 + "".join(" " if a == b else ("-" if "-" in (a, b) else "v") for a, b in zip(qs_, ss_)))
        blk.append(f'{"C " if minus else "  "}{c["name"] + "#" + cl + "/" + fam:<14}{sp:>8} {ss_} {s_end}')
        blk.append("")
        qp, sp = q_end + 1, s_end + step
    blk += ["Matrix = 25p41g.matrix", "Kimura (with divCpGMod) = 10.00", "Transitions / transversions = 1.00 (1/1)",
            "Gap_init rate = 0.01 (1 / 100), avg. gap size = 1.00 (1 / 1)", "", ""]
    alig.extend(blk)
    seq = c["qa"].replace("-", "")
    for k, cp in enumerate(c["truth"]):
        if cp is None:
            continue
        gpos = c["start"] + k if not minus else c["end"] - 1 - k
        truth_rows.append(f"{c['chrom']}\t{gpos}\t{cp}")
with gzip.open(f"{OUT}/rmsk.txt.gz", "wt") as fh:
    fh.write("\n".join(ucsc) + "\n")
with gzip.open(f"{OUT}/rmsk.out.gz", "wt") as fh:
    fh.write("   SW  perc perc perc  query  position in query  matching  repeat  position in repeat\n"
             "score  div. del. ins.  sequence  begin end (left)  repeat  class/family  begin end (left)  ID\n\n")
    fh.write("\n".join(outl) + "\n")
with gzip.open(f"{OUT}/rmsk.align.gz", "wt") as fh:
    fh.write("\n".join(alig) + "\n")
with gzip.open(f"{OUT}/truth_mapping.tsv.gz", "wt") as fh:
    fh.write("chrom\tpos\tcons_pos_true\n" + "\n".join(truth_rows) + "\n")

# ------------------------------------------------------------------ phyloP, gnomAD, AN
bw = pyBigWig.open(f"{OUT}/phylop.bw", "w")
bw.addHeader(list(CHROMS.items()))
vcf_rows, an_rows = [], []
truth_by = {}
for r in truth_rows:
    ch, p, cp = r.split("\t")
    truth_by[(ch, int(p))] = int(cp)
for c in copies:
    ch = c["chrom"]
    gseq = "".join(genome[ch][c["start"] - 1:c["end"] + 1])
    lo, hi = PLANT.get((c["name"], c["comp"], c["ori"]), (0, -1))
    vals = rng.normal(0, 1, c["end"] - c["start"]).astype(float)
    for k, p in enumerate(range(c["start"], c["end"])):
        cp = truth_by.get((ch, p))
        planted = cp is not None and lo <= cp <= hi
        if planted:
            vals[k] += 1.5
        b = gseq[k + 1]
        cpg = (b == "C" and gseq[k + 2] == "G") or (b == "G" and gseq[k] == "C")
        prob = (0.45 if cpg else 0.15) * (0.5 if planted else 1.0)
        an = 50_000 if rng.random() < 0.02 else 152_000
        an_rows.append((ch, p + 1, an))
        if rng.random() < prob:
            ac = 1 if rng.random() < (0.75 if planted else 0.5) else int(rng.integers(2, 500))
            alt = rng.choice([x for x in "ACGT" if x != b])
            vcf_rows.append((ch, p + 1, b, alt, ac, an))
    bw.addEntries(ch, c["start"], values=list(vals), span=1, step=1)
bw.close()
vcf_rows.sort()
an_rows.sort()
with open(f"{OUT}/gnomad.vcf", "w") as fh:
    fh.write("##fileformat=VCFv4.2\n")
    for c, n in CHROMS.items():
        fh.write(f"##contig=<ID={c},length={n}>\n")
    fh.write('##INFO=<ID=AC,Number=A,Type=Integer,Description="AC">\n'
             '##INFO=<ID=AN,Number=1,Type=Integer,Description="AN">\n'
             '##FILTER=<ID=PASS,Description="All filters passed">\n'
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for ch, p, ref, alt, ac, an in vcf_rows:
        fh.write(f"{ch}\t{p}\t.\t{ref}\t{alt}\t.\tPASS\tAC={ac};AN={an}\n")
pysam.tabix_index(f"{OUT}/gnomad.vcf", preset="vcf", force=True)
# also one file per chromosome whose header lists ALL chromosomes, like the real gnomAD files
for c in CHROMS:
    with open(f"{OUT}/gnomad_{c}.vcf", "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        for cc, n in CHROMS.items():
            fh.write(f"##contig=<ID={cc},length={n}>\n")
        fh.write('##INFO=<ID=AC,Number=A,Type=Integer,Description="AC">\n'
                 '##INFO=<ID=AN,Number=1,Type=Integer,Description="AN">\n'
                 '##FILTER=<ID=PASS,Description="All filters passed">\n'
                 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for ch, p, ref, alt, ac, an in vcf_rows:
            if ch == c:
                fh.write(f"{ch}\t{p}\t.\t{ref}\t{alt}\t.\tPASS\tAC={ac};AN={an}\n")
    pysam.tabix_index(f"{OUT}/gnomad_{c}.vcf", preset="vcf", force=True)

with open(f"{OUT}/an.tsv", "w") as fh:
    for ch, p, an in an_rows:
        fh.write(f"{ch}\t{p}\t{an}\n")
pysam.tabix_index(f"{OUT}/an.tsv", seq_col=0, start_col=1, end_col=1, zerobased=False, force=True)

# ------------------------------------------------------------------ positive controls + tiny MAF
dropped = {genes[2]["gid"]} | {g["gid"] for g in genes if g["gid"].endswith("05.1")}
ex_sense = [c for c in copies if c["name"] == "AluSx" and c["comp"] == "exon" and c["ori"] == "sense"
            and c["cs"] == 1 and c["gid"] not in dropped]
ex_anti = [c for c in copies if c["name"] == "AluSx" and c["comp"] == "exon" and c["ori"] == "antisense"
           and c["gid"] not in dropped]
with open(f"{OUT}/controls.bed", "w") as fh:
    for nm, c in (("planted_sense_copy", ex_sense[0]), ("antisense_copy", ex_anti[0])):
        fh.write(f'{c["chrom"]}\t{c["start"]}\t{c["end"]}\t{nm}\n')
c = copies[0]
ref = "".join(genome[c["chrom"]][c["start"]:c["start"] + 40])
with open(f"{OUT}/test.maf", "w") as fh:
    fh.write("##maf version=1\n\na score=0\n")
    fh.write(f's hg38.{c["chrom"]} {c["start"]} 40 + {CHROMS[c["chrom"]]} {ref[:20]}--{ref[20:]}\n')
    fh.write(f's panTro6.chr1 100 42 + 1000000 {ref[:20]}AA{ref[20:]}\n')
    fh.write(f's mm39.chr5 100 30 + 1000000 {"-" * 10}{ref[10:20]}--{ref[20:]}\n\n')
print("done", file=sys.stderr)
