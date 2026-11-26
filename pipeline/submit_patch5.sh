#!/bin/bash
# Submits the patch-5 analyses (stage s14) as a chain of Slurm jobs. Needs the finished main run and the revision
# run (work/rev/: annot_sub, junction_usage, splice_age_sites). Run from the project directory or anywhere:
#     DRY=1 bash pipeline/submit_patch5.sh        # print the sbatch commands only
#     bash pipeline/submit_patch5.sh              # submit
#     SKIP="lprep lgn" bash pipeline/submit_patch5.sh   # resume: skip finished stages (keys below)
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
P=$PIPE/s14_patch5.sbatch
LP=$(sub lprep - --array=1-2 $P loose_prep)
LG=$(sub lgn $LP --array=1-92 --mem=4G --time=12:00:00 --cpus-per-task=1 $P loose_gnomad)
LA=$(sub lannot $LG --array=1-2 $P loose_annot)
LX=$(sub lanalyze $LA --array=1-3 $P loose_analyze)
SC=$(sub scan - --array=1-8 $P scan)
US=$(sub usage - --mem=16G --time=04:00:00 $P usage_sum)
BN=$(sub bnull - --time=12:00:00 $P boundary_null)
AC=$(sub ageca - --mem=16G --time=04:00:00 $P age_ca)
BS=$(sub bsub - --mem=16G --time=08:00:00 $P bounds_sub)
FB=$(sub featbg - --time=12:00:00 $P feature_bg)
RP=$(sub report5 $LX:$SC:$US:$BN:$AC:$BS:$FB --mem=4G --time=01:00:00 $P report5)
if [ -n "$DRY" ]; then echo "dry run only; nothing submitted"; else echo "submitted; final report job: $RP"; fi
