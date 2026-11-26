#!/bin/bash
#SBATCH --job-name=ridl_s12_ltr_l1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%A_%a.log
# LTR and L1 elements through the same pipeline as Alu, MIR and L2 (outputs: work/step2_*, step3_*, annot_*,
# work/rev/ltr, work/rev/l1). Usage: sbatch s12_ltr_l1.sbatch prepare|map|annotate ; --array=1-2 analyze
#   prepare   copies of LTR (THE1, MST, MLT1, LTR12, LTR7, MLT2) and old L1 (L1M*) subfamilies (step 2)
#   map       exact consensus mapping (3' model for L1) and grouping of close subfamilies (steps 3, 31, 3b)
#   annotate  per-base annotation (step 4; needs s03r ltr and l1) and class-level relabelled copies (step 28)
#   analyze   task 1 LTR, task 2 L1: class-level whole-element tests (all and same genes), GC-matched
#             conservation, gene-level TE-specific difference and bounds, decomposition; position profiles and
#             bounds for groups with >= 200 exonic copies in the annotation
#   map_lin, annotate_lin, analyze_lin   L1 with linear consensus coordinates (all copies). Most old L1M
#             subfamilies have only 5'/ORF2 models in the RepeatMasker alignments, so the exact (3' model) mapping
#             keeps ~45% of L1 copies; class-level whole-element tests do not need consensus positions, so they are
#             repeated on all copies (outputs work/step3_l1lin, work/annot_l1lin, work/rev/l1lin)
source "$PIPE/config_rev.sh"
class_level () {  # $1 class label (LTR/L1), $2 class-relabelled annotation, $3 output dir
    local CL=$1 A=$2 O=$3 P PA
    mkdir -p $O/prof $O/comp $O/trend
    stamp "$CL: class-level tests (all subfamilies pooled; strata include the original subfamily)"
    for P in all paired; do
        PA=(); [ $P = paired ] && PA=(--paired-genes)
        python $S/06_profiles.py --annot $A --subfamilies $CL --bin-size 100000 "${PA[@]}" $STRATA_TE \
            --n-perm $N_PERM --n-boot 1000 --no-plots --out-dir $O/prof/class_$P
    done
    python $S/14_composition.py --annot $A --subfamilies $CL --paired-genes --out $O/comp/te_${CL}_paired
    python $S/12_gene_trend.py --annot $A --subfamilies $CL --gene-table $W/gene_props.tsv \
        --pool-orientations --n-boot 200 --baseline $W/trend/nonte_genes.tsv --out $O/trend/class
    python $S/23_bounds.py --genes class=$O/trend/class_genes.tsv --comp-dir $O/comp --n-boot 2000 --out $O/bounds_class
    python $S/29_gene_dropout.py --annot $A --subfamilies $CL --lnc-genes $W/step1_all/lnc_genes.tsv.gz \
        --gene-table $W/gene_props.tsv --n-boot 1000 --out $O/dropout_class
}
case "$1" in
prepare)
    python $S/02_select_te_copies.py $TE_SELECT_LTR --compartments $W/step1_all/compartments.tsv.gz --out-dir $W/step2_ltr
    python $S/02_select_te_copies.py $TE_SELECT_L1 --compartments $W/step1_all/compartments.tsv.gz --out-dir $W/step2_l1 ;;
map)
    for X in ltr l1; do
        mkdir -p $W/step3_$X
        python $S/03_map_to_consensus.py --copies $W/step2_$X/te_copies.tsv.gz --mode align --align-file $D/hg38.fa.align.gz \
            --consensus-out $W/step3_$X/consensus_rec.fa --out $W/step3_$X/map_align.parquet
        if [ $X = ltr ]; then mapfile -t F < <(fam_args "$FAM_LTR_STR"); else mapfile -t F < <(fam_args "$FAM_L1_STR"); fi
        mapfile -t GR < <(python $S/31_make_groups.py --copies $W/step2_$X/te_copies.tsv.gz \
            --consensus $W/step3_$X/consensus_rec.fa "${F[@]}" --report $W/step3_$X/groups_auto.tsv)
        printf '  %s\n' "${GR[@]}"
        python $S/03b_group_subfamilies.py --bases $W/step3_$X/map_align.parquet --copies $W/step2_$X/te_copies.tsv.gz \
            --consensus $W/step3_$X/consensus_rec.fa "${GR[@]}" --out-bases $W/step3_$X/map_grouped.parquet \
            --out-copies $W/step3_$X/copies_grouped.tsv.gz --out-consensus $W/step3_$X/consensus_grouped.fa \
            --report $W/step3_$X/groups.tsv
    done ;;
annotate)
    for X in ltr l1; do
        n=$(ls $W/gnomad_$X/gn.chr*.vcf.gz.tbi | wc -l); [ "$n" -eq $N_CHR_GNOMAD ] || { echo "gnomAD subsets for $X: $n of $N_CHR_GNOMAD"; exit 1; }
        annotate $W/step3_$X/map_grouped.parquet $W/step3_$X/copies_grouped.tsv.gz $W/step1_all $W/gnomad_$X \
            $W/annot_$X $W/step3_$X/consensus_grouped.fa
        python $S/28_annot_subset.py --annot $W/annot_$X --subfamilies ALL --copies $W/step3_$X/copies_grouped.tsv.gz \
            --splice-sites $W/step1_all/splice_sites.tsv.gz --query "$AUTO" --relabel ${X^^} --out-dir $RW/annot_${X}_class
    done ;;
analyze)
    X=$([ "$SLURM_ARRAY_TASK_ID" -eq 1 ] && echo ltr || echo l1); CL=${X^^}; O=$RW/$X
    class_level $CL $RW/annot_${X}_class $O
    stamp "$CL: groups with >= $MIN_EXONIC_GROUP exonic copies"
    # exonic autosomal units per group in the annotation (copies without consensus positions are not there)
    mapfile -t GL < <(python - "$S" "$W/annot_$X" "$MIN_EXONIC_GROUP" << 'PY'
import sys
sys.path.insert(0, sys.argv[1])
from ridl_stats import load_annot
d = load_annot(sys.argv[2], columns=["unit_id", "subfamily", "compartment", "chrom"])
d = d[(d["compartment"] == "exon") & (d["chrom"] != "chrX")]
n = d.groupby("subfamily")["unit_id"].nunique()
for g in n.index[n >= int(sys.argv[3])]:
    print(g)
PY
)
    printf '  %s\n' "${GL[@]}"
    for G in "${GL[@]}"; do
        python $S/06_profiles.py --annot $W/annot_$X --subfamilies "$G" --bin-size 20 --paired-genes --query "$AUTO" $STRATA_TE \
            --n-perm $N_PERM --n-boot 200 --no-plots --out-dir $O/prof/${G}_b20_paired
        python $S/14_composition.py --annot $W/annot_$X --subfamilies "$G" --paired-genes --query "$AUTO" \
            --out $O/comp/te_${G}_paired
    done
    if [ ${#GL[@]} -gt 0 ]; then
        python $S/12_gene_trend.py --annot $W/annot_$X --subfamilies "$(IFS=,; echo "${GL[*]}")" --gene-table $W/gene_props.tsv \
            --query "$AUTO" --pool-orientations --n-boot 200 --baseline $W/trend/nonte_genes.tsv --out $O/trend/groups
        python $S/23_bounds.py --genes groups=$O/trend/groups_genes.tsv --comp-dir $O/comp --n-boot 2000 --out $O/bounds_groups
    fi ;;
map_lin)
    X=l1; O=$W/step3_${X}lin; mkdir -p $O
    python $S/03_map_to_consensus.py --copies $W/step2_$X/te_copies.tsv.gz --mode linear --out $O/map_grouped.parquet
    python - $W/step2_$X/te_copies.tsv.gz $O/copies_grouped.tsv.gz << 'PY'
import sys
import pandas as pd
c = pd.read_csv(sys.argv[1], sep="\t", dtype={"chrom": str})
c["subfamily_orig"] = c["subfamily"]
c.to_csv(sys.argv[2], sep="\t", index=False)
print(f"copies {len(c)}, subfamilies {c['subfamily'].nunique()}, exonic {int((c['compartment'] == 'exon').sum())}")
PY
    ;;
annotate_lin)
    X=l1
    n=$(ls $W/gnomad_$X/gn.chr*.vcf.gz.tbi | wc -l); [ "$n" -eq $N_CHR_GNOMAD ] || { echo "gnomAD subsets for $X: $n of $N_CHR_GNOMAD"; exit 1; }
    annotate $W/step3_${X}lin/map_grouped.parquet $W/step3_${X}lin/copies_grouped.tsv.gz $W/step1_all $W/gnomad_$X $W/annot_${X}lin
    python $S/28_annot_subset.py --annot $W/annot_${X}lin --subfamilies ALL --copies $W/step3_${X}lin/copies_grouped.tsv.gz \
        --splice-sites $W/step1_all/splice_sites.tsv.gz --query "$AUTO" --relabel ${X^^} --out-dir $RW/annot_${X}lin_class ;;
analyze_lin)
    class_level L1 $RW/annot_l1lin_class $RW/l1lin ;;
*) echo "unknown part: $1"; exit 1 ;;
esac
stamp "done $1"
