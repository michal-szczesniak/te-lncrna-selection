#!/bin/bash
#SBATCH --job-name=ridl_s16_patch6
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%A_%a.log
# Patch 6: analyses requested before submission, and the final figures and tables.
# Patch 11: render also draws Fig. 5 (example loci, make_fig_locus.py, from results/locus of patch 9);
# the former Fig. 5 and 6 are saved as Fig6 and Fig7.
# Outputs in work/rev/p6/, results/revision6/ and results/manuscript/. Usage: sbatch s16_patch6.sbatch <part>
#   --array=1-2 canon_filter  exonic compartment restricted to exons of the Ensembl canonical transcript (step 41):
#                             1 TE annotation, 2 repeat-free annotation
#   canon_summary             composition (step 14) and gene-level exon-intron differences (step 12) on the filtered
#                             annotations, as for the main set
#   canon_bounds              limits on TE-specific constraint with canonical exons (step 23)
#   bounds_well               limits in lncRNAs with independent support: genes with a GTEx estimate, GENCODE 29
#                             lincRNAs, both (step 37 --well-supported)
#   tss_ltr                   orientation and consensus position of lncRNA starts in LTR elements (step 42)
#   tss_ltr_sub               patch 16: clustering of LTR starts within single subfamilies (step 42b), from the
#                             step-42 table, to results/revision6/tss_ltr_subfamilies.tsv
#   facts                     counts quoted in the text without a table of their own (step 43)
#   report6                   tables to results/revision6/ and results/revision6_tables.tgz
#   synthetic                 synthetic-data runs for Supplemental Fig. S1 (as in stage s10)
#   render                    all main and supplemental figures, Table 1 and Supplemental Tables S1-S9 from the
#                             results, to results/manuscript/ (scripts in manuscript/figures, manuscript/tables)
source "$PIPE/config_rev.sh"
t=${SLURM_ARRAY_TASK_ID:-0}
P6=$RW/p6
mkdir -p $P6
case "$1" in
canon_filter)
    if [ $t -eq 1 ]; then A=$W/annot; O=$P6/annot_canon; else A=$W/annot_nonte; O=$P6/annot_nonte_canon; fi
    python $S/41_canonical_filter.py --annot $A --gtf $GTF --out-dir $O ;;
canon_summary)
    mkdir -p $P6/comp_canon $P6/trend_canon
    python $S/14_composition.py --annot $P6/annot_nonte_canon --paired-genes --out $P6/comp_canon/lnc_paired
    for G in "${MAIN_GROUPS[@]}"; do
        python $S/14_composition.py --annot $P6/annot_canon --subfamilies "$G" --query "$AUTO" --paired-genes \
            --out "$P6/comp_canon/te_${G//\*/}_paired"
    done
    python $S/12_gene_trend.py --annot $P6/annot_nonte_canon --subfamilies nonTE --gene-table $W/gene_props.tsv \
        --query "$AUTO" --n-boot 500 --out $P6/trend_canon/nonte
    python $S/12_gene_trend.py --annot $P6/annot_canon --subfamilies 'AluS*,MIR*,L2abc,AluJ*' \
        --gene-table $W/gene_props.tsv --query "$AUTO" --pool-orientations --n-boot 500 \
        --baseline $P6/trend_canon/nonte_genes.tsv --out $P6/trend_canon/te_vs_nonte ;;
canon_bounds)
    python $S/23_bounds.py --genes canonical=$P6/trend_canon/te_vs_nonte_genes.tsv --comp-dir $P6/comp_canon \
        --n-boot 2000 --out $P6/canonical_bounds ;;
bounds_well)
    python $S/37_bounds_subsets.py --genes $W/trend/te_vs_nonte_genes.tsv --comp-dir $W/comp --n-boot 2000 \
        --well-supported --gene-class $W/gene_class.tsv --skip-expression-subsets --out $P6/bounds_well ;;
tss_ltr)
    python $S/42_tss_ltr.py --gtf $GTF --lnc-genes $W/step1_all/lnc_genes.tsv.gz --rmsk-out $D/hg38.fa.out.gz \
        --n-perm 1000 --out $P6/tss_ltr ;;
tss_ltr_sub)
    T=$P6/tss_ltr_tss.tsv; [ -s "$T" ] || T=results/revision6/tss_ltr_tss.tsv
    mkdir -p results/revision6
    python $S/42b_tss_ltr_subfamilies.py --tss $T --n-perm 1000 --out results/revision6/tss_ltr_subfamilies.tsv ;;
facts)
    python $S/43_facts.py --copies $W/step2_v2/te_copies.tsv.gz --copies-grouped $W/step3/copies_grouped.tsv.gz \
        --step1 $W/step1_all --attrition $W/attrition.tsv --out $P6/facts ;;
report6)
    O=results/revision6; mkdir -p $O
    for f in $P6/canonical_bounds_*.tsv $P6/annot_canon/canonical_summary.tsv $P6/bounds_well_*.tsv \
             $P6/tss_ltr_summary.tsv $P6/tss_ltr_hist.tsv $P6/tss_ltr_tss.tsv $P6/facts.tsv; do
        [ -s "$f" ] && cp "$f" $O/
    done
    (cd $O && md5sum * > MD5SUMS && ls -l)
    tar czf results/revision6_tables.tgz -C results revision6
    ls -l results/revision6_tables.tgz ;;
synthetic)
    W=$(cd $W && pwd)
    # the synthetic-data runs of stage s10 (Fig. S1), redone when work/synthetic is missing or incomplete
    rm -rf $W/synthetic $W/synthetic_shift
    bash $S/tests/run_synthetic_test.sh $W/synthetic
    TMP=$(mktemp -d)
    cp -r $S $TMP/scripts
    cp manuscript/figures/make_synthetic_data_shift.py $TMP/scripts/tests/make_synthetic_data.py
    # the assertions at the end of the test expect no deficit outside the planted segments, so they fail here by design
    SHIFT=0.75 bash $TMP/scripts/tests/run_synthetic_test.sh $W/synthetic_shift || stamp "shifted run: assertions failed as expected"
    test -s $W/synthetic/work/prof_main/profiles.tsv
    test -s $W/synthetic_shift/work/prof_main/profiles.tsv
    rm -rf $TMP ;;
render)
    M=manuscript; O=results/manuscript; FAIL=""
    mkdir -p $O/figures $O/tables
    export MPLCONFIGDIR=$PWD/work/.mpl
    run () { local k=$1; shift; "$@" || { stamp "FAILED: $k"; FAIL="$FAIL $k"; }; }
    SYN=$W/synthetic; SHIFT=$W/synthetic_shift
    [ -s $SYN/work/prof_main/profiles.tsv ] || { stamp "no $SYN/work/prof_main: Fig. S1 needs the synthetic runs of stage s10"; SYN=none; SHIFT=none; }
    rm -f $O/figures/Fig5.* $O/figures/Fig6.* $O/figures/Fig7.*   # patch 11: Fig. 5 is the locus figure, old 5 and 6 are 6 and 7
    run main_figures env RIDL_REV5=results/revision5 python $M/figures/make_figures.py . results/revision $O/figures
    run supp_figures python $M/figures/make_supp_figures.py . results/revision $SYN $SHIFT $O/figures
    run figS3b python $M/figures/make_figS3b.py results/revision5b $O/figures/FigS3B
    [ -s results/revision6/tss_ltr_summary.tsv ] && run figS8 python $M/figures/make_figS8.py results/revision6 $O/figures/FigS8
    [ -s results/locus/bundle.json ] && run fig5_locus python $M/figures/make_fig_locus.py results/locus/bundle.json $O/figures \
        "chr6:42943970-42947607+a" "chr11:59563906-59565014-a"
    run table1 python $M/tables/make_table1.py --project . --out $O/tables/table1.md
    run supp_tables python $M/tables/make_supp_tables.py --project . --out $O/tables
    (cd $O && find . -type f ! -name MD5SUMS -exec md5sum {} + > MD5SUMS && ls -lR)
    tar czf results/manuscript_outputs.tgz -C results manuscript
    ls -l results/manuscript_outputs.tgz
    [ -z "$FAIL" ] || { stamp "render finished with failures:$FAIL"; exit 1; } ;;
*) echo "unknown part: $1"; exit 1 ;;
esac
stamp "done $1 ${t}"
