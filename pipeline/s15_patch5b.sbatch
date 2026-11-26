#!/bin/bash
#SBATCH --job-name=ridl_s15_patch5b
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%A_%a.log
# Patch 5b (outputs in work/rev/p5b/ and results/revision5b/). Usage: sbatch s15_patch5b.sbatch <part>
#   --array=1-8 scan      spike-in power along the element with genes resampled (step 33), now with the permutation
#                         count and replicates of the main analysis ($N_PERM, $SPIKE_REPS); patch 5 used 300
#                         permutations, too few for BH q < 0.1 on a single 20-nt bin
#   age_ca                splice-site ages with the copy veto (step 36): --rule absent and --rule strict
#   cons_nonte --array=1-3  step-15 unit tables of repeat-free windows (step 40; their annotation has no
#                         subfamily_orig column, which step 15 requires): 1 tss0, 2 open, 3 main
#   cons_pooled           exon vs intron conserved stretches pooled per TE group, against repeat-free sequence,
#                         main / tss0 / open (step 39)
#   report5b              tables to results/revision5b/ and results/revision5b_tables.tgz
source "$PIPE/config_rev.sh"
tag () { echo "${1//\*/}"; }
t=${SLURM_ARRAY_TASK_ID:-0}
P5=$RW/p5
P5B=$RW/p5b
mkdir -p $P5B
VARIANTS=(tss0 open main)
case "$1" in
scan)
    G=${G4[$(((t - 1) % 4))]}; ORI=$([ $t -le 4 ] && echo sense || echo antisense); T=$(tag "$G")
    mkdir -p $P5B/power_scan
    python $S/33_spikein_scan.py --annot $W/annot --subfamilies "$G" --orientation $ORI --scan 20 \
        --reductions 0,0.1,0.25,0.5 --fractions 1,0.5 --reps $SPIKE_REPS --n-perm $N_PERM --boot-genes --query "$AUTO" \
        --out $P5B/power_scan/${T}_$ORI ;;
age_ca)
    for R in absent strict; do
        O=$P5B/splice_age_ca; [ $R = strict ] && O=${O}_strict
        python $S/36_age_copyaware.py --sites $RW/splice_age_sites.tsv.gz --species "$(IFS=,; echo "${SPECIES[*]}")" \
            --rule $R --out $O
    done ;;
cons_nonte)
    V=${VARIANTS[$((t - 1))]}
    A=$W/annot_nonte_$V; [ "$V" = main ] && A=$W/annot_nonte
    mkdir -p $P5B/cons
    python $S/40_conserved_units.py --annot $A --out $P5B/cons/${V}_nonte_conserved ;;
cons_pooled)
    TEU=$W/side/conserved_units.tsv.gz
    if [ ! -s $TEU ]; then
        mkdir -p $P5B/cons
        python $S/15_conserved_fragments.py --annot $W/annot --out $P5B/cons/main_te_conserved
        TEU=$P5B/cons/main_te_conserved_units.tsv.gz
    fi
    python $S/39_conserved_pooled.py --variant main $TEU $P5B/cons/main_nonte_conserved_units.tsv.gz \
        --variant tss0 $P5/loose/tss0_conserved_units.tsv.gz $P5B/cons/tss0_nonte_conserved_units.tsv.gz \
        --variant open $P5/loose/open_conserved_units.tsv.gz $P5B/cons/open_nonte_conserved_units.tsv.gz \
        --groups "$(IFS=,; echo "${G4[*]},${GX[*]}")" --n-boot 1000 --out $P5B/conserved_pooled ;;
report5b)
    O=results/revision5b; mkdir -p $O
    for f in $P5B/power_scan/*_power.tsv $P5B/splice_age_ca*_summary.tsv $P5B/splice_age_ca*_reclassified.tsv \
             $P5B/cons/*_summary.tsv $P5B/conserved_pooled.tsv; do
        [ -s "$f" ] && cp "$f" $O/
    done
    (cd $O && md5sum * > MD5SUMS && ls -l)
    tar czf results/revision5b_tables.tgz -C results revision5b
    ls -l results/revision5b_tables.tgz ;;
*) echo "unknown part: $1"; exit 1 ;;
esac
stamp "done $1 ${t}"
