#!/bin/bash
# Submits the revision analyses (stages s00r, s03r, s11, s12) as a chain of Slurm jobs. Needs a finished main run
# (work/ of the published pipeline). Run from anywhere:
#     bash pipeline/submit_revision.sh              # submit
#     DRY=1 bash pipeline/submit_revision.sh        # print the sbatch commands only
#     SKIP="dl power" bash pipeline/submit_revision.sh   # resume: skip finished stages (keys below)
set -eo pipefail
export PIPE=$(cd "$(dirname "$0")" && pwd)
[ -f "$PIPE/local.sh" ] && source "$PIPE/local.sh"
PROJECT=${PROJECT:-$(cd "$PIPE/.." && pwd)}
export PROJECT
cd "$PROJECT"
mkdir -p logs

sub () {  # sub <key> <dependency ids joined by ':'> <sbatch args...>; prints the job id or '-' if skipped
    local key=$1 dep=$2; shift 2
    if [[ " $SKIP " == *" $key "* ]]; then echo "skip $key" >&2; echo "-"; return; fi
    dep=$(echo "$dep" | tr ':' '\n' | { grep -v '^-$' || true; } | paste -sd: -)
    local d=""
    [ -n "$dep" ] && d="--dependency=afterok:$dep"
    if [ -n "$DRY" ]; then
        echo "sbatch $d $*" >&2; echo "$key"
    else
        sbatch --parsable $d "$@"
    fi
}
P=$PIPE
DL=$(sub dl - $P/s00r_download.sbatch)
SUB=$(sub subsets - $P/s11_revision.sbatch subsets)
PRO=$(sub profiles $SUB --array=1-15 $P/s11_revision.sbatch profiles)
TRS=$(sub trend_sub $SUB $P/s11_revision.sbatch trend_sub)
BND=$(sub bounds $TRS $P/s11_revision.sbatch bounds)
POW=$(sub power - --array=1-8 $P/s11_revision.sbatch power)
POS=$(sub posctrl - $P/s11_revision.sbatch posctrl)
PQ=$(sub posctrl_q - $P/s11_revision.sbatch posctrl_q)
DRP=$(sub dropout - $P/s11_revision.sbatch dropout)
CON=$(sub conserved - $P/s11_revision.sbatch conserved)
SPC=$(sub splice_ctrl - $P/s11_revision.sbatch splice_ctrl)
AGE=$(sub age $DL $P/s11_revision.sbatch age)
AGB=$(sub age_base $AGE $P/s11_revision.sbatch age_base)
USE=$(sub usage $DL $P/s11_revision.sbatch usage)
STR=$(sub strata $SPC:$AGE:$USE $P/s11_revision.sbatch strata)
LP=$(sub l_prep - $P/s12_ltr_l1.sbatch prepare)
LM=$(sub l_map $LP $P/s12_ltr_l1.sbatch map)
GL=$(sub g_ltr $LP --array=1-23 $P/s03r_gnomad.sbatch ltr)
G1=$(sub g_l1 $LP --array=1-23 $P/s03r_gnomad.sbatch l1)
LA=$(sub l_annot $LM:$GL:$G1 $P/s12_ltr_l1.sbatch annotate)
LX=$(sub l_analyze $LA --array=1-2 $P/s12_ltr_l1.sbatch analyze)
LML=$(sub l_maplin $LP $P/s12_ltr_l1.sbatch map_lin)
LAL=$(sub l_annlin $LML:$G1 $P/s12_ltr_l1.sbatch annotate_lin)
LXL=$(sub l_anlin $LAL $P/s12_ltr_l1.sbatch analyze_lin)
REP=$(sub report $PRO:$BND:$POW:$POS:$PQ:$DRP:$CON:$STR:$AGB:$LX:$LXL $P/s11_revision.sbatch report)
if [ -n "$DRY" ]; then echo "dry run only; nothing submitted"; else echo "submitted; final report job: $REP"; fi
