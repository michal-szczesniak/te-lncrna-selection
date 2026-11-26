#!/bin/bash
#SBATCH --job-name=ridl_s10_report
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%j.log
# Table S6 summary, synthetic-data runs for Fig. S1, Figures 1-6 and S1-S4.
source "$PIPE/config.sh"
mkdir -p results/figures results/supplementary results/tables

stamp "Table S6: main analysis and robustness runs"
python $S/compare_checks.py | tee results/tables/robustness_runs.txt
cp $W/attrition.tsv $W/attrition_genes.tsv results/tables/ 2>/dev/null || true

stamp "synthetic data with planted constraint and with a uniform exonic deficit (Fig. S1)"
rm -rf $W/synthetic $W/synthetic_shift
bash $S/tests/run_synthetic_test.sh $W/synthetic
TMP=$(mktemp -d)
cp -r $S $TMP/scripts
cp figures/make_synthetic_data_shift.py $TMP/scripts/tests/make_synthetic_data.py
# the assertions at the end of the test expect no deficit outside the planted segments, so they fail here by design
SHIFT=0.75 bash $TMP/scripts/tests/run_synthetic_test.sh $W/synthetic_shift || stamp "shifted run: assertions failed as expected"
test -s $W/synthetic_shift/work/prof_main/profiles.tsv
rm -rf $TMP

stamp "figures"
python figures/make_figures.py . results/figures
python figures/make_supp_figures.py . $W/synthetic $W/synthetic_shift results/supplementary
ls -l results/figures results/supplementary
stamp "done"
