#!/bin/bash
#SBATCH --job-name=ridl_s06_summaries
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%j.log
# Usage: sbatch s06_summaries.sbatch a|b
#   a: baselines (ordinary sequence, protein-coding control), GC composition (Fig. 2), gene-level excess and
#      RNA-property trends (Fig. 3-4)
#   b: conserved stretches per copy, transcription-start check (Fig. S2), conserved fragments and localization
#      (Fig. 4a-b), TE contributions to lncRNA architecture (Fig. 5)
source "$PIPE/config.sh"
case "$1" in
a)
    C="--annot $W/annot_nonte --subfamilies nonTE --bin-size 100 --pop-strata tri,cpg --n-perm 1000 --n-boot 200 --no-plots"
    python $S/06_profiles.py $C --out-dir $W/prof_nonte/exon_vs_intron
    python $S/06_profiles.py $C --paired-genes --out-dir $W/prof_nonte/exon_vs_intron_paired
    python $S/06_profiles.py $C --fg intron --bg intergenic --out-dir $W/prof_nonte/intron_vs_intergenic
    C="--annot $W/annot_pc --subfamilies nonTE --bin-size 100 --pop-strata tri,cpg --n-perm 300 --n-boot 100 --no-plots"
    python $S/06_profiles.py $C --out-dir $W/prof_pc/exon_vs_intron
    python $S/06_profiles.py $C --paired-genes --out-dir $W/prof_pc/exon_vs_intron_paired

    stamp "step 14: GC composition"
    mkdir -p $W/comp
    python $S/14_composition.py --annot $W/annot_nonte --out $W/comp/lnc_all
    python $S/14_composition.py --annot $W/annot_nonte --paired-genes --out $W/comp/lnc_paired
    python $S/14_composition.py --annot $W/annot_pc --out $W/comp/pc_all
    python $S/14_composition.py --annot $W/annot_pc --paired-genes --out $W/comp/pc_paired
    for G in "${MAIN_GROUPS[@]}"; do
        python $S/14_composition.py --annot $W/annot --subfamilies "$G" --query "$AUTO" --paired-genes \
            --out "$W/comp/te_${G//\*/}_paired"
    done

    stamp "step 12: gene-level excess and RNA-property trends"
    mkdir -p $W/trend
    python $S/12_gene_trend.py --annot $W/annot_nonte --subfamilies nonTE --gene-table $W/gene_props.tsv \
        --query "$AUTO" --n-boot 500 --out $W/trend/nonte
    python $S/12_gene_trend.py --annot $W/annot --subfamilies 'AluS*,MIR*,L2abc,AluJ*' --gene-table $W/gene_props.tsv \
        --query "$AUTO" --pool-orientations --n-boot 500 --baseline $W/trend/nonte_genes.tsv --out $W/trend/te_vs_nonte ;;
b)
    stamp "step 10: conserved stretches per copy"
    mkdir -p $W/copycons
    python $S/10_copy_conservation.py --annot $W/annot --subfamilies 'MIR*,L2abc,AluS*,AluJ*' --query "$AUTO" \
        --paired-genes --threshold 1.0 --out $W/copycons/paired
    python $S/10_copy_conservation.py --annot $W/annot --subfamilies 'MIR*,L2abc,AluS*,AluJ*' --query "$AUTO" \
        --paired-genes --out $W/copycons/paired_win

    stamp "step 13: conservation by distance from transcription start (Fig. S2)"
    python $S/13_exon_intron_by_tss.py --ref-tx $W/step1_all/ref_tx.tsv.gz --rmsk-out $D/hg38.fa.out.gz $PHYLOP \
        --max-genes 4000 --out $W/tss_check.tsv
    python $S/13_exon_intron_by_tss.py --ref-tx $W/step1/ref_tx.tsv.gz --rmsk-out $D/hg38.fa.out.gz $PHYLOP \
        --max-genes 6000 --out $W/tss_check_strict.tsv

    stamp "steps 15-16: conserved fragments and localization"
    mkdir -p $W/side
    python $S/15_conserved_fragments.py --annot $W/annot --out $W/side/conserved
    L="--ref-tx $W/step1_all/ref_tx.tsv.gz --rmsk-out $D/hg38.fa.out.gz --genome-fasta $D/hg38.fa --units $W/side/conserved_units.tsv.gz"
    python $S/16_te_localization.py $L --gene-table $W/gene_props.tsv --expr-col synth --out $W/side/localization
    for cell in HepG2 K562; do
        python $S/16_te_localization.py $L --gene-table $RNA_PROPS --nuc-col nuc_cyto_log2_$cell \
            --expr-col synthesis_log2rpkm_$cell --out $W/side/loc_$cell
    done
    python $S/16_te_localization.py $L --gene-table $W/gene_props.tsv --min-exons 2 --out $W/side/loc_multiexon
    python $S/16_te_localization.py $L --gene-table $W/gene_props.tsv --cons-thr 1.0 --out $W/side/loc_cons1

    stamp "step 17: TE contributions to lncRNA architecture"
    python $S/17_te_architecture.py --gtf $GTF --gene-set lnc=$W/step1_all/lnc_genes.tsv.gz \
        --gene-set pc=$W/step1_pc/lnc_genes.tsv.gz --rmsk-out $D/hg38.fa.out.gz --out $W/side/architecture ;;
*) echo "unknown part: $1"; exit 1 ;;
esac
stamp "done"
