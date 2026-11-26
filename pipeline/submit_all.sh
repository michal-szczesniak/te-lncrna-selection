#!/bin/bash
# Submits the whole analysis as a chain of Slurm jobs. Run from anywhere:
#     bash pipeline/submit_all.sh            # submit
#     DRY=1 bash pipeline/submit_all.sh      # print the sbatch commands only
#     DOWNLOAD=1 bash pipeline/submit_all.sh # start with stage 0 (data download)
#     SKIP="s01 g_te" bash pipeline/submit_all.sh   # resume: skip stages that already finished (keys below)
# Settings are in pipeline/config.sh; machine-specific paths in pipeline/local.sh (see pipeline/local.sh.example).
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

S0=-
[ -n "$DOWNLOAD" ] && S0=$(sub s00 - $P/s00_download.sbatch)
S1=$(sub s01 $S0 $P/s01_prepare.sbatch)
S2=$(sub s02 $S1 $P/s02_map.sbatch)
GTE=$(sub g_te $S1 --array=1-23 $P/s03_gnomad.sbatch te)
GNT=$(sub g_nonte $S1 --array=1-23 $P/s03_gnomad.sbatch nonte)
GPC=$(sub g_pc $S1 --array=1-22 $P/s03_gnomad.sbatch pc)
ATE=$(sub a_te $S2:$GTE $P/s04_annotate.sbatch te)
ANT=$(sub a_nonte $GNT $P/s04_annotate.sbatch nonte)
APC=$(sub a_pc $GPC $P/s04_annotate.sbatch pc)
PRO=$(sub s05 $ATE:$ANT $P/s05_profiles.sbatch)
S6A=$(sub s06a $ATE:$ANT:$APC $P/s06_summaries.sbatch a)
S6B=$(sub s06b $ATE $P/s06_summaries.sbatch b)
SP1=$(sub s07p $S2 $P/s07_splice.sbatch prepare)
GSS=$(sub g_ss $SP1 --array=1-23 $P/s03_gnomad.sbatch ss)
SP2=$(sub s07a $GSS:$GTE $P/s07_splice.sbatch analyze)
T1=$(sub s08p $S1 $P/s08_tss100.sbatch prepare)
GT=$(sub g_tss100 $T1 --array=1-23 $P/s03_gnomad.sbatch tss100)
T2=$(sub s08a $GT:$S6A $P/s08_tss100.sbatch analyze)
E1=$(sub s09p $S1 $P/s09_ext.sbatch prepare)
GE1=$(sub g_ext_te $E1 --array=1-23 $P/s03_gnomad.sbatch ext_te)
GE2=$(sub g_ext_nonte $E1 --array=1-23 $P/s03_gnomad.sbatch ext_nonte)
E2=$(sub s09a $GE1:$GE2 $P/s09_ext.sbatch annotate)
E3=$(sub s09prof $E2 --array=1-8 $P/s09_ext.sbatch profiles)
E4=$(sub s09s $E2 $P/s09_ext.sbatch summary)
REP=$(sub s10 $PRO:$S6A:$S6B:$SP2:$T2:$E3:$E4 $P/s10_report.sbatch)
if [ -n "$DRY" ]; then echo "dry run only; nothing submitted"; else echo "submitted; final report job: $REP"; fi
