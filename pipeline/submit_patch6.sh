#!/bin/bash
# Submits patch 6 (stage s16) as a chain of Slurm jobs, ending with the final figures and tables.
#     DRY=1 bash pipeline/submit_patch6.sh        # print the sbatch commands only
#     bash pipeline/submit_patch6.sh              # submit
#     SKIP="canon tss" bash pipeline/submit_patch6.sh   # resume: skip finished stages (keys below)
#     ONLY_RENDER=1 bash pipeline/submit_patch6.sh      # only redraw figures and tables from existing results
#     SYNTH=1 ONLY_RENDER=1 bash pipeline/submit_patch6.sh   # first redo the synthetic runs for Fig. S1
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
P=$PIPE/s16_patch6.sbatch
if [ -n "$ONLY_RENDER" ]; then
    SY=-
    [ -n "$SYNTH" ] && SY=$(sub synthetic - --mem=16G --time=04:00:00 $P synthetic)
    RD=$(sub render "$SY" --mem=16G --time=02:00:00 $P render)
    [ -n "$DRY" ] && echo "dry run only; nothing submitted" || echo "submitted: $RD"
    exit 0
fi
CF=$(sub canon - --array=1-2 --mem=48G --time=06:00:00 $P canon_filter)
CS=$(sub canonsum $CF --time=12:00:00 $P canon_summary)
CB=$(sub canonbounds $CS --mem=16G --time=04:00:00 $P canon_bounds)
BW=$(sub well - --mem=16G --time=04:00:00 $P bounds_well)
TL=$(sub tss - --mem=24G --time=04:00:00 $P tss_ltr)
FA=$(sub facts - --mem=24G --time=02:00:00 $P facts)
RP=$(sub report6 $CB:$BW:$TL:$FA --mem=4G --time=01:00:00 $P report6)
RD=$(sub render $RP --mem=16G --time=02:00:00 $P render)
if [ -n "$DRY" ]; then echo "dry run only; nothing submitted"; else echo "submitted; final job: $RD"; fi
