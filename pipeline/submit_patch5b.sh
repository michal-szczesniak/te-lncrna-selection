#!/bin/bash
# Submits patch 5b (stage s15) as a chain of Slurm jobs. Needs the finished patch-5 run (work/rev/p5/loose).
#     DRY=1 bash pipeline/submit_patch5b.sh        # print the sbatch commands only
#     bash pipeline/submit_patch5b.sh              # submit
#     SKIP="scan" bash pipeline/submit_patch5b.sh  # resume: skip finished stages (keys below)
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
P=$PIPE/s15_patch5b.sbatch
SC=$(sub scan - --array=1-8 $P scan)
AC=$(sub ageca - --mem=16G --time=04:00:00 $P age_ca)
CN=$(sub consnonte - --array=1-3 --mem=32G --time=08:00:00 $P cons_nonte)
CP=$(sub conspool $CN --mem=32G --time=08:00:00 $P cons_pooled)
RP=$(sub report5b $SC:$AC:$CP --mem=4G --time=01:00:00 $P report5b)
if [ -n "$DRY" ]; then echo "dry run only; nothing submitted"; else echo "submitted; final report job: $RP"; fi
