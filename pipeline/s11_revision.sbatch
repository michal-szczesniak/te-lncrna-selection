#!/bin/bash
#SBATCH --job-name=ridl_s11_revision
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%A_%a.log
# Revision analyses on the main data set (outputs in work/rev/). Usage: sbatch s11_revision.sbatch <part>
#   subsets                      unit classes, Alu structure, consensus-match layer (step 28)
#   --array=1-15 profiles        same-gene profiles: embedded exonic copies (1-4), >= 50 nt from splice sites (5-8),
#                                Alu dimers (9-10), all copies with the consensus-match layer (11-15)
#   trend_sub                    gene-level TE-specific difference for embedded copies (step 12)
#   bounds                       bounds with 2,000 resamples, conservation layers, orientation contrast (step 23)
#   --array=1-8 power            spike-in power, 4 groups x 2 orientations (step 22)
#   posctrl                      TE positive control: conserved stretches and human variants (step 24)
#   posctrl_q                    the same with quantile thresholds (top 1%/5% of windows; Alu with primate phyloP)
#   dropout                      decomposition of the all-genes deficit, genes dropping out (step 29)
#   conserved                    conserved stretches exon vs intron, same genes and GC-matched (step 30)
#   splice_ctrl                  location-matched controls for TE-derived splice sites (steps 18b, 19)
#   age                          ancestral state and age of splice sites (step 25; needs s00r)
#   age_base                     neutral baseline for the age: unused copies, intronic dinucleotides (step 25b; needs age)
#   usage                        GTEx v10 junction usage (step 26; needs s00r)
#   strata                       splice-site constraint by usage, age and control set (step 27)
#   report                       copies the main tables to results/revision/
source "$PIPE/config_rev.sh"
tag () { echo "${1//\*/}"; }
t=${SLURM_ARRAY_TASK_ID:-0}
case "$1" in
subsets)
    for G in "${G4[@]}" "${GX[@]}"; do
        python $S/28_annot_subset.py --annot $W/annot --subfamilies "$G" --copies $W/step3/copies_grouped.tsv.gz \
            --splice-sites $W/step1_all/splice_sites.tsv.gz --query "$AUTO" --out-dir $RW/annot_sub/$(tag "$G")
    done ;;
profiles)
    if [ $t -le 4 ]; then G=${G4[$((t - 1))]}; Q="unit_class != 'boundary'"; N=embedded
    elif [ $t -le 8 ]; then G=${G4[$((t - 5))]}; Q="dist_splice >= 50"; N=far50
    elif [ $t -le 10 ]; then G=${GALU[$((t - 9))]}; Q="alu_structure == 'dimer'"; N=dimer
    else G6=("${G4[@]}" "${GX[0]}"); G=${G6[$((t - 11))]}; Q=""; N=all; fi
    [ -z "$G" ] && { stamp "task $t: no group"; exit 0; }
    T=$(tag "$G")
    QA=(); [ -n "$Q" ] && QA=(--query "$Q")
    python $S/06_profiles.py --annot $RW/annot_sub/$T --subfamilies "$G" --bin-size 20 --paired-genes "${QA[@]}" \
        $STRATA_TE --n-perm $N_PERM --n-boot 200 --no-plots --out-dir $RW/prof_sub/${T}_$N ;;
trend_sub)
    mkdir -p $RW/trend_sub
    for G in "${G4[@]}"; do
        T=$(tag "$G")
        python $S/12_gene_trend.py --annot $RW/annot_sub/$T --subfamilies "$G" --gene-table $W/gene_props.tsv \
            --query "unit_class != 'boundary'" --pool-orientations --n-boot 200 --baseline $W/trend/nonte_genes.tsv \
            --out $RW/trend_sub/${T}_embedded
    done ;;
bounds)
    A=(--genes main=$W/trend/te_vs_nonte_genes.tsv)
    [ -s $W/trend_tss100/te_vs_nonte_genes.tsv ] && A+=(--genes tss100=$W/trend_tss100/te_vs_nonte_genes.tsv)
    [ -s $W/trend_ext/te_vs_nonte_genes.tsv ] && A+=(--genes extended=$W/trend_ext/te_vs_nonte_genes.tsv)
    for G in "${G4[@]}"; do
        f=$RW/trend_sub/$(tag "$G")_embedded_genes.tsv
        [ -s $f ] && A+=(--genes embedded_$(tag "$G")=$f)
    done
    python $S/23_bounds.py "${A[@]}" --comp-dir $W/comp --n-boot 2000 --out $RW/bounds ;;
power)
    G=${G4[$(((t - 1) % 4))]}; ORI=$([ $t -le 4 ] && echo sense || echo antisense); T=$(tag "$G")
    mkdir -p $RW/power
    SEG=${SPIKE_SEGS[$G]:-$SPIKE_DEFAULT}
    python $S/22_spikein_power.py --annot $W/annot --subfamilies "$G" --orientation $ORI --segments "$SEG" \
        --reductions 0,0.1,0.25,0.5 --fractions 1,0.5,0.25 --reps $SPIKE_REPS --n-perm $N_PERM --query "$AUTO" \
        --out $RW/power/${T}_$ORI ;;
posctrl)
    mkdir -p $RW/posctrl
    for thr in 1.0 1.5; do
        python $S/24_te_positive_control.py --annot $W/annot --subfamilies "$OLD_LIST" --track p447 \
            --threshold $thr --query "$AUTO" --out $RW/posctrl/old_p447_t$thr
        python $S/24_te_positive_control.py --annot $W/annot --subfamilies "$ALU_LIST" --track prim \
            --threshold $thr --query "$AUTO" --out $RW/posctrl/alu_prim_t$thr
    done ;;
posctrl_q)
    mkdir -p $RW/posctrl
    for q in 0.99 0.95; do
        python $S/24_te_positive_control.py --annot $W/annot --subfamilies "$ALU_LIST" --track prim \
            --threshold-quantile $q --query "$AUTO" --out $RW/posctrl/alu_prim_q$q
    done
    python $S/24_te_positive_control.py --annot $W/annot --subfamilies "$OLD_LIST" --track p447 \
        --threshold-quantile 0.99 --query "$AUTO" --out $RW/posctrl/old_p447_q0.99 ;;
dropout)
    mkdir -p $RW/dropout
    python $S/29_gene_dropout.py --annot $W/annot --subfamilies "$(IFS=,; echo "${G4[*]},${GX[*]}")" --query "$AUTO" \
        --lnc-genes $W/step1_all/lnc_genes.tsv.gz --gene-table $W/gene_props.tsv --n-boot 2000 --out $RW/dropout/te
    python $S/29_gene_dropout.py --annot $W/annot_nonte --subfamilies nonTE --pop-strata tri,cpg --query "$AUTO" \
        --lnc-genes $W/step1_all/lnc_genes.tsv.gz --gene-table $W/gene_props.tsv --n-boot 2000 --out $RW/dropout/nonte ;;
conserved)
    python $S/30_conserved_paired.py --units $W/side/conserved_units.tsv.gz --genome-fasta $D/hg38.fa --track p447 \
        --threshold 1.5 --families '^(?:MIR|L2)' --out $RW/conserved_p447
    python $S/30_conserved_paired.py --units $W/side/conserved_units.tsv.gz --genome-fasta $D/hg38.fa --track prim \
        --threshold 1.0 --families '^Alu' --out $RW/conserved_alu_prim ;;
splice_ctrl)
    for P in intron samegene; do
        python $S/18b_splice_controls.py --positions $W/splice/positions.tsv.gz --copies $W/step3/copies_grouped.tsv.gz \
            --map $W/step3/map_grouped.parquet --genome-fasta $D/hg38.fa --pool $P --out-dir $RW/splice_$P
        python $S/19_splice_analyze.py --positions $RW/splice_$P/positions.tsv.gz $(gnomad_args $W/gnomad) \
            --an-tsv $D/gnomad_AN_v4.1.tsv.bgz $PHYLOP --genome-fasta $D/hg38.fa --n-boot 1000 --out $RW/splice_$P/res
    done ;;
age)
    python $S/25_splice_age.py --positions $W/splice/res_positions.tsv.gz --map-align $W/step3/map_align.parquet \
        --consensus $W/step3/consensus_rec.fa --copies $W/step3/copies_grouped.tsv.gz --genome-fasta $D/hg38.fa \
        $(species_args) --out $RW/splice_age ;;
age_base)
    python $S/25b_age_baseline.py --sites $RW/splice_age_sites.tsv.gz --positions $W/splice/res_positions.tsv.gz \
        --map $W/step3/map_grouped.parquet --copies $W/step3/copies_grouped.tsv.gz --rmsk-out $D/hg38.fa.out.gz \
        --genome-fasta $D/hg38.fa $(species_args) --out $RW/splice_age_base ;;
usage)
    AT=(); [ -s $GTEX_ATTR ] && AT=(--attributes $GTEX_ATTR)
    python $S/26_junction_usage.py --gct $GTEX_JUNC "${AT[@]}" --positions $W/splice/res_positions.tsv.gz \
        --out $RW/junction_usage.tsv.gz ;;
strata)
    ALT=(--positions-alt intronic=$RW/splice_intron/res_positions.tsv.gz --positions-alt samegene=$RW/splice_samegene/res_positions.tsv.gz)
    for M in total_reads donor_share acceptor_share; do
        python $S/27_splice_strata.py --positions $W/splice/res_positions.tsv.gz --usage $RW/junction_usage.tsv.gz \
            --age $RW/splice_age_sites.tsv.gz "${ALT[@]}" --usage-metric $M --n-boot 1000 --out $RW/splice_strata_$M
    done ;;
report)
    O=results/revision; mkdir -p $O
    for f in $RW/bounds_*.tsv $RW/power/*_power.tsv $RW/posctrl/*.tsv $RW/dropout/*.tsv $RW/conserved_*_paired.tsv \
             $RW/splice_age_summary.tsv $RW/splice_age_base_baseline.tsv $RW/splice_strata_*_strata.tsv $RW/splice_intron/res_stats.tsv \
             $RW/splice_samegene/res_stats.tsv $RW/splice_*/summary.txt $RW/ltr/bounds_*.tsv $RW/l1/bounds_*.tsv \
             $RW/ltr/dropout_*.tsv $RW/l1/dropout_*.tsv $RW/l1lin/bounds_*.tsv $RW/l1lin/dropout_*.tsv \
             $RW/ltr/comp/*_summary.tsv $RW/l1/comp/*_summary.tsv $RW/l1lin/comp/*_summary.tsv \
             $RW/ltr/trend/*_trend.tsv $RW/l1/trend/*_trend.tsv $RW/l1lin/trend/*_trend.tsv; do
        [ -s "$f" ] && cp "$f" $O/$(echo ${f#$RW/} | tr '/' '_')
    done
    for d in $RW/prof_sub/* $RW/ltr/prof/* $RW/l1/prof/* $RW/l1lin/prof/*; do
        [ -s $d/whole_element.tsv ] && cp $d/whole_element.tsv $O/whole_$(echo ${d#$RW/} | tr '/' '_').tsv
        [ -s $d/profiles.tsv ] && cp $d/profiles.tsv $O/profiles_$(echo ${d#$RW/} | tr '/' '_').tsv
    done
    cp $W/step3_ltr/groups.tsv $O/groups_ltr.tsv 2>/dev/null || true
    cp $W/step3_l1/groups.tsv $O/groups_l1.tsv 2>/dev/null || true
    cp $W/step3_ltr/groups_auto.tsv $O/groups_auto_ltr.tsv 2>/dev/null || true
    cp $W/step3_l1/groups_auto.tsv $O/groups_auto_l1.tsv 2>/dev/null || true
    # exonic units by orientation and class (copy contains a splice site of its gene = 'boundary'), all genes
    python - "$S" $RW/annot_sub $O/unit_classes.tsv << 'PY'
import glob, os, sys
sys.path.insert(0, sys.argv[1])
import pandas as pd
from ridl_stats import load_annot
rows = []
for d in sorted(glob.glob(os.path.join(sys.argv[2], "*"))):
    if not glob.glob(os.path.join(d, "*.parquet")):
        continue
    a = load_annot(d, columns=["unit_id", "subfamily", "compartment", "orientation", "unit_class", "gene_id"])
    u = a[a["compartment"] == "exon"].drop_duplicates("unit_id")
    t = u.groupby(["subfamily", "orientation", "unit_class"]).agg(units=("unit_id", "size"), genes=("gene_id", "nunique"))
    rows.append(t.reset_index())
if rows:
    out = pd.concat(rows, ignore_index=True)
    out["share_of_orientation"] = out["units"] / out.groupby(["subfamily", "orientation"])["units"].transform("sum")
    out.to_csv(sys.argv[3], sep="\t", index=False)
    print(out.round(3).to_string(index=False))
PY
    (cd $O && md5sum * > MD5SUMS && ls -l)
    tar czf results/revision_tables.tgz -C results revision
    ls -l results/revision_tables.tgz ;;
*) echo "unknown part: $1"; exit 1 ;;
esac
stamp "done $1 ${t}"
