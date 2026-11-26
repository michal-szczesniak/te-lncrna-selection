#!/bin/bash
#SBATCH --job-name=ridl_s13_flank_gc
#SBATCH --array=1-8
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%A_%a.log
# GC matching on the flanks of each unit (step 32), as a robustness check of the GC-matched conservation in
# Fig. 3a; same genes (paired), autosomes. Outputs: work/rev/flankgc/<set>_paired_{summary,bins}.tsv and
# _units.tsv.gz; the summaries are copied to results/revision/flankgc_<set>_paired_summary.tsv.
#   1 repeat-free lncRNA windows   2 protein-coding windows   3-6 Alu S, Alu J, MIR, L2 a/b/c   7 LTR   8 L1
# Submit (login node is fine):  export PIPE=$PWD/pipeline; sbatch --array=1-8 pipeline/s13_flank_gc.sbatch
source "$PIPE/config_rev.sh"
t=${SLURM_ARRAY_TASK_ID:-1}
O=$RW/flankgc
mkdir -p $O results/revision
FLANK=${FLANK:-1000}
case $t in
1) N=lnc;    A=(--annot $W/annot_nonte --subfamilies nonTE) ;;
2) N=pc;     A=(--annot $W/annot_pc --subfamilies nonTE) ;;
3) N=AluS;   A=(--annot $W/annot --subfamilies 'AluS*' --query "$AUTO") ;;
4) N=AluJ;   A=(--annot $W/annot --subfamilies 'AluJ*' --query "$AUTO") ;;
5) N=MIR;    A=(--annot $W/annot --subfamilies 'MIR*' --query "$AUTO") ;;
6) N=L2abc;  A=(--annot $W/annot --subfamilies 'L2abc' --query "$AUTO") ;;
7) N=LTR;    A=(--annot $RW/annot_ltr_class --subfamilies LTR) ;;
8) N=L1;     A=(--annot $RW/annot_l1lin_class --subfamilies L1) ;;
*) echo "task $t: nothing to do"; exit 1 ;;
esac
if [ ! -e "${A[1]}" ]; then echo "missing annotation ${A[1]} for $N"; exit 1; fi
stamp "flank GC ($FLANK nt), set $N"
python $S/32_flank_gc.py "${A[@]}" --paired-genes --genome-fasta $D/hg38.fa --flank $FLANK \
    --out $O/${N}_paired
cp $O/${N}_paired_summary.tsv results/revision/flankgc_${N}_paired_summary.tsv
stamp "task $t done"
