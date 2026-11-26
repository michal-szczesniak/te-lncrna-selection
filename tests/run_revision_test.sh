#!/usr/bin/env bash
# End-to-end test of the revision stages (s11, s12) on synthetic data, run with bash instead of sbatch.
# Builds the synthetic data set of the published pipeline, a fake project tree with the inputs the revision
# stages expect, two fake 'species' and a GTEx-like junction table, then runs every part and checks outputs.
# Usage: bash tests/run_revision_test.sh <bundle dir with scripts/, pipeline/, tests/> [work dir]
set -euo pipefail
B=$(cd "$1" && pwd)
T=${2:-$(mktemp -d)}
mkdir -p "$T"
T=$(cd "$T" && pwd)
echo "test directory: $T"

# 1. synthetic data and the main pipeline steps used by the revision (steps 1-4, 12, 14, 15, 18, 19)
bash $B/scripts/tests/run_synthetic_test.sh $T/base > $T/base.log 2>&1
SYN=$T/base/syn; BW=$T/base/work
python3 $B/tests/make_splice_test.py $SYN $BW $T/ss > $T/ss.log
P=$T/proj; mkdir -p $P/data $P/work $P/logs
ln -sfn $B/scripts $P/scripts; ln -sfn $B/pipeline $P/pipeline
D=$P/data
ln -sf $T/ss/genome_ss.fa $D/hg38.fa; ln -sf $T/ss/genome_ss.fa.fai $D/hg38.fa.fai
ln -sf $SYN/rmsk.out.gz $D/hg38.fa.out.gz; ln -sf $SYN/rmsk.align.gz $D/hg38.fa.align.gz
ln -sf $SYN/phylop.bw $D/hg38.phyloP447way.bw; ln -sf $SYN/phylop.bw $D/hg38.phyloP447wayPrimates.bw
ln -sf $SYN/an.tsv.gz $D/gnomad_AN_v4.1.tsv.bgz; ln -sf $SYN/an.tsv.gz.tbi $D/gnomad_AN_v4.1.tsv.bgz.tbi
mkdir -p $D/rev
for sp in sp1 sp2; do
    cap="$(tr '[:lower:]' '[:upper:]' <<< ${sp:0:1})${sp:1}"
    gzip -c $T/ss/hg38To$cap.chain > $D/rev/hg38To$cap.over.chain.gz
    ln -sf $T/ss/$sp.fa $D/rev/$sp.fa; ln -sf $T/ss/$sp.fa.fai $D/rev/$sp.fa.fai
done
ln -sf $T/ss/junctions.gct.gz $D/rev/GTEx_Analysis_v10_STARv2.7.10a_junctions.gct.gz
ln -sf $T/ss/attributes.txt $D/rev/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt
W=$P/work
cp -r $BW/step1 $W/step1_all; cp -r $BW/step2 $W/step2_v2
mkdir -p $W/step3 $W/gnomad
cp $BW/map_align.parquet $BW/map_grouped.parquet $BW/copies_grouped.tsv.gz $BW/consensus_rec.fa $BW/consensus_grouped.fa $W/step3/
for c in chr1 chr2; do cp $SYN/gnomad_$c.vcf.gz $W/gnomad/gn.$c.vcf.gz; cp $SYN/gnomad_$c.vcf.gz.tbi $W/gnomad/gn.$c.vcf.gz.tbi; done
echo 'conda () { :; }' > $T/conda_stub.sh

export PROJECT=$P PIPE=$P/pipeline CONDA_SH=$T/conda_stub.sh CONDA_ENV=none
export REV_GROUPS_STR="AluS*,MIRb" REV_EXTRA_STR="AluS*" REV_ALU_STR="AluS*,AluS*" REV_OLD_LIST="MIRb" REV_ALU_LIST="AluS*"
export SPECIES_STR="sp1,sp2" N_PERM=100 SPIKE_REPS=2 SPIKE_DEFAULT="101-120,all" N_CHR_GNOMAD=2 MIN_EXONIC_GROUP=50
export TE_SELECT_LTR="--rmsk $D/hg38.fa.out.gz --rmsk-format rmout --subfamily-regex ^AluS --assign perbase --max-intronic 300 --max-intergenic 100"
export TE_SELECT_L1="--rmsk $D/hg38.fa.out.gz --rmsk-format rmout --subfamily-regex ^MIR --assign perbase --max-intronic 300 --max-intergenic 100"
export FAM_LTR_STR="AluS_g:^AluS" FAM_L1_STR="MIR_g:^MIR"
cd $P
source $P/pipeline/config.sh   # functions annotate, gnomad_args and variables of the main pipeline

# annotation with both phyloP tracks, as on the cluster
annotate $W/step3/map_grouped.parquet $W/step3/copies_grouped.tsv.gz $W/step1_all $W/gnomad $W/annot $W/step3/consensus_grouped.fa 2> $T/annot.log
python3 $S/28_annot_subset.py --annot $W/annot --subfamilies MIRb --copies $W/step3/copies_grouped.tsv.gz \
    --splice-sites $W/step1_all/splice_sites.tsv.gz --relabel nonTE --out-dir $W/annot_nonte 2>> $T/annot.log
python3 - $W << 'PY'
import sys, numpy as np, pandas as pd
W = sys.argv[1]
g = pd.read_csv(f"{W}/step1_all/lnc_genes.tsv.gz", sep="\t"); g["gene_id"] = g["gene_id"].str.split(".").str[0]
r = np.random.default_rng(1)
for p in ["stab", "nuc", "splice", "synth", "tau"]:
    g[p] = r.normal(size=len(g))
g[["gene_id", "stab", "nuc", "splice", "synth", "tau"]].to_csv(f"{W}/gene_props.tsv", sep="\t", index=False)
PY
mkdir -p $W/trend $W/comp $W/side $W/splice
python3 $S/12_gene_trend.py --annot $W/annot_nonte --subfamilies nonTE --gene-table $W/gene_props.tsv --n-boot 20 --out $W/trend/nonte 2> /dev/null
python3 $S/12_gene_trend.py --annot $W/annot --subfamilies 'AluS*,MIRb' --gene-table $W/gene_props.tsv --pool-orientations \
    --n-boot 20 --baseline $W/trend/nonte_genes.tsv --out $W/trend/te_vs_nonte 2> /dev/null
python3 $S/14_composition.py --annot $W/annot --subfamilies 'AluS*' --paired-genes --out $W/comp/te_AluS_paired 2> /dev/null
python3 $S/14_composition.py --annot $W/annot --subfamilies MIRb --paired-genes --out $W/comp/te_MIRb_paired 2> /dev/null
python3 $S/15_conserved_fragments.py --annot $W/annot --out $W/side/conserved 2> /dev/null
python3 $S/18_splice_prepare.py --gtf $SYN/gencode.gtf.gz --lnc-genes $W/step1_all/lnc_genes.tsv.gz --pc-genes $W/step1_all/lnc_genes.tsv.gz \
    --rmsk-out $D/hg38.fa.out.gz --genome-fasta $D/hg38.fa --copies $W/step3/copies_grouped.tsv.gz --map $W/step3/map_grouped.parquet \
    --te-regions $W/step2_v2/gnomad_regions.bed --chroms chr1,chr2 --max-introns 300 --out-dir $W/splice 2> /dev/null
python3 $S/19_splice_analyze.py --positions $W/splice/positions.tsv.gz $(gnomad_args $W/gnomad) --an-tsv $D/gnomad_AN_v4.1.tsv.bgz \
    $PHYLOP --genome-fasta $D/hg38.fa --n-boot 20 --out $W/splice/res 2> /dev/null
echo "inputs ready"

# 2. the revision stages
run () { echo "== $*"; SLURM_ARRAY_TASK_ID=${TASK:-0} bash $P/pipeline/$1 "${@:2}" > $T/log_$(echo "$*" | tr ' /' '__').txt 2>&1 || { echo "FAILED: $*"; tail -20 $T/log_$(echo "$*" | tr ' /' '__').txt; exit 1; }; }
run s11_revision.sbatch subsets
for i in 1 2 5 6 9 11 12; do TASK=$i run s11_revision.sbatch profiles; done
run s11_revision.sbatch trend_sub
run s11_revision.sbatch bounds
for i in 1 6; do TASK=$i run s11_revision.sbatch power; done
run s11_revision.sbatch posctrl
run s11_revision.sbatch posctrl_q
run s11_revision.sbatch dropout
run s11_revision.sbatch conserved
run s11_revision.sbatch splice_ctrl
run s11_revision.sbatch age
run s11_revision.sbatch age_base
run s11_revision.sbatch usage
run s11_revision.sbatch strata
run s12_ltr_l1.sbatch prepare
run s12_ltr_l1.sbatch map
for X in ltr l1; do mkdir -p $W/gnomad_$X; cp $W/gnomad/gn.chr* $W/gnomad_$X/; done
run s12_ltr_l1.sbatch annotate
for i in 1 2; do TASK=$i run s12_ltr_l1.sbatch analyze; done
run s12_ltr_l1.sbatch map_lin
run s12_ltr_l1.sbatch annotate_lin
run s12_ltr_l1.sbatch analyze_lin
run s11_revision.sbatch report

# 3. checks
python3 - $W/rev $P/results/revision << 'PY'
import os, sys
import pandas as pd
RW, O = sys.argv[1], sys.argv[2]
need = ["bounds_bounds.tsv", "bounds_functional.tsv", "power/AluS_sense_power.tsv", "posctrl/old_p447_t1.0_posctrl.tsv",
        "dropout/te_decomposition.tsv", "conserved_p447_paired.tsv", "splice_intron/res_stats.tsv",
        "splice_samegene/res_stats.tsv", "splice_age_summary.tsv", "junction_usage.tsv.gz",
        "splice_strata_total_reads_strata.tsv", "prof_sub/AluS_embedded/whole_element.tsv",
        "prof_sub/AluS_all/whole_element.tsv", "ltr/bounds_class_bounds.tsv", "l1/prof/class_paired/whole_element.tsv",
        "l1lin/bounds_class_bounds.tsv", "l1lin/prof/class_paired/whole_element.tsv",
        "posctrl/alu_prim_q0.99_posctrl.tsv", "posctrl/old_p447_q0.99_posctrl.tsv", "splice_age_base_baseline.tsv"]
miss = [f for f in need if not os.path.exists(os.path.join(RW, f))]
assert not miss, f"missing outputs: {miss}"
w = pd.read_csv(os.path.join(RW, "prof_sub/AluS_all/whole_element.tsv"), sep="\t")
assert "phylop_consmatch" in set(w["layer"]), "consensus-match layer missing"
p = pd.read_csv(os.path.join(RW, "power/AluS_sense_power.tsv"), sep="\t")
assert (p.loc[(p.segment == "all") & (p.reduction == 0.5) & (p.fraction_units == 1), "power_whole"] == 1).all()
pc = pd.read_csv(os.path.join(RW, "posctrl/old_p447_t1.0_posctrl.tsv"), sep="\t")
assert (pc["mh_log_or"].dropna() < 0).all(), "planted conservation should come with fewer variants"
ag = pd.read_csv(os.path.join(RW, "splice_age_summary.tsv"), sep="\t")
x = ag[(ag.set == "lnc_nonTE_site") & (ag.subset == "all") & (ag.measure == "present_in_sp1 (of aligned)")]["share"].iloc[0]
assert 0.35 < x < 0.65, f"sp1 destroys half of the sites; present share {x:.2f}"
assert len(os.listdir(O)) > 20
uc = pd.read_csv(os.path.join(O, "unit_classes.tsv"), sep="\t")
assert {"boundary", "embedded"} <= set(uc["unit_class"]), "unit classes missing"
q = pd.read_csv(os.path.join(RW, "posctrl/alu_prim_q0.99_posctrl.tsv"), sep="\t")
assert q["threshold_quantile"].eq(0.99).all() and q["share_bases_in_stretch"].max() > 0
ll = pd.read_csv(os.path.join(RW, "l1lin/prof/class_paired/whole_element.tsv"), sep="\t")
ab = pd.read_csv(os.path.join(RW, "splice_age_base_baseline.tsv"), sep="\t")
assert {"lnc_TE_site", "lnc_nonTE_site", "pc_site"} <= set(ab["set"]), "baseline sets missing"
st = pd.read_csv(os.path.join(RW, "splice_strata_total_reads_strata.tsv"), sep="\t")
assert st["comparison"].str.contains("strata x CpG").any(), "CpG-stratified rows missing"
assert len(ll) > 0
print("ALL REVISION CHECKS PASSED")
PY
