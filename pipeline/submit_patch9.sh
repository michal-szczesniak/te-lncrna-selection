#!/bin/bash
# Submits patch 9 (stage s17): candidate loci for the illustrative figure.
#     DRY=1 bash pipeline/submit_patch9.sh    # print the sbatch command only
#     bash pipeline/submit_patch9.sh          # submit
set -eo pipefail
export PIPE=$(cd "$(dirname "$0")" && pwd)
[ -f "$PIPE/local.sh" ] && source "$PIPE/local.sh"
PROJECT=${PROJECT:-$(cd "$PIPE/.." && pwd)}
export PROJECT
cd "$PROJECT"
mkdir -p logs
if [ -n "$DRY" ]; then
    echo "sbatch $PIPE/s17_patch9.sbatch locus"
else
    echo "submitted: $(sbatch --parsable $PIPE/s17_patch9.sbatch locus)"
fi
