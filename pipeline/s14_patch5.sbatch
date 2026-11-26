#!/bin/bash
#SBATCH --job-name=ridl_s14_patch5
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%A_%a.log
# Patch 5 analyses (outputs in work/rev/p5/ and results/revision5/). Usage: sbatch s14_patch5.sbatch <part>
#   --array=1-2 loose_prep      compartments without our exclusions (1 = tss0: no margins around transcript starts
#                               and ends; 2 = open: also no splice-site margins, no overlap filter, sequence shared
#                               with other genes and segmental duplications kept), TE copies, mapping, repeat-free windows
#   --array=1-92 loose_gnomad   gnomAD subsets for the two variants x (TE, repeat-free) x 23 chromosomes
#   --array=1-2 loose_annot     base annotation of the two variants (TE and repeat-free)
#   --array=1-3 loose_analyze   exon-intron conservation, all genes and same genes, unmatched and GC-matched
#                               (step 32), and conserved stretches by compartment (step 15): 1 tss0, 2 open, 3 main
#   --array=1-8 scan            spike-in power scanned along the element with genes resampled (step 33)
#   usage_sum                   reads behind the usage tertiles; junctions absent from the GTEx table (step 34)
#   boundary_null               chance expectation for exonic copies that contain an exon boundary (step 35)
#   age_ca                      splice-site ages with presence requiring the copy itself (step 36)
#   bounds_sub                  gene-level limits in highly and broadly expressed lncRNAs (step 37)
#   feature_bg                  TEs at lncRNA features vs sequence next to the features (step 38)
#   report5                     tables to results/revision5/ and results/revision5_tables.tgz
source "$PIPE/config_rev.sh"
tag () { echo "${1//\*/}"; }
t=${SLURM_ARRAY_TASK_ID:-0}
P5=$RW/p5
mkdir -p $P5
VARIANTS=(tss0 open)
AUTOS=$(printf "chr%s," {1..22} | sed 's/,$//')
case "$1" in
loose_prep)
    V=${VARIANTS[$((t - 1))]}
    if [ "$V" = tss0 ]; then
        python $S/01_build_compartments.py $STEP1_COMMON --keep-unmatched-gtex --exon-tss-pad 0 --intron-tss-pad 0 \
            --tes-pad 0 --out-dir $W/step1_$V
    else
        python $S/01_build_compartments.py --gtf $GTF --chrom-sizes $D/hg38.chrom.sizes --gtex-gct $GCT \
            --keep-unmatched-gtex --exon-tss-pad 0 --intron-tss-pad 0 --tes-pad 0 --exon-splice-pad 0 \
            --intron-splice-pad 0 --skip-overlap-filter --keep-other-gene-spans --out-dir $W/step1_$V
    fi
    python $S/02_select_te_copies.py $TE_SELECT --compartments $W/step1_$V/compartments.tsv.gz --out-dir $W/step2_$V
    map_and_group $W/step2_$V $W/step3_$V
    python $S/09_nonte_windows.py --compartments $W/step1_$V/compartments.tsv.gz $NONTE_WINDOWS --chroms $AUTOS \
        --out-dir $W/nonte_$V ;;
loose_gnomad)
    i=$((t - 1)); V=${VARIANTS[$((i / 46))]}; k=$(((i % 46) / 23)); c=$((i % 23 + 1))
    if [ $k -eq 0 ]; then BED=$W/step2_$V/gnomad_regions.bed; OUT=$W/gnomad_${V}_te
    else BED=$W/nonte_$V/gnomad_regions.bed; OUT=$W/gnomad_${V}_nonte; fi
    if [ $c -eq 23 ]; then CHR=chrX; else CHR=chr$c; fi
    URL=$GNOMAD_URL/gnomad.genomes.v4.1.sites.$CHR.vcf.bgz
    mkdir -p $OUT
    awk -v c=$CHR '$1==c' $BED > $OUT/regions_$CHR.bed
    stamp "$V $(basename $OUT) $CHR: $(wc -l < $OUT/regions_$CHR.bed) regions"
    if [ ! -s $OUT/regions_$CHR.bed ]; then
        bcftools view -h $URL | bcftools annotate -x ^INFO/AC,INFO/AN -Oz -o $OUT/gn.$CHR.vcf.gz
    else
        for attempt in 1 2 3; do
            if bcftools view -R $OUT/regions_$CHR.bed -f PASS -v snps $URL \
                | bcftools annotate -x ^INFO/AC,INFO/AN -Oz -o $OUT/gn.$CHR.vcf.gz; then
                break
            fi
            stamp "attempt $attempt failed, retrying in 60 s"
            sleep 60
            [ "$attempt" -eq 3 ] && exit 1
        done
    fi
    bcftools index -f -t $OUT/gn.$CHR.vcf.gz ;;
loose_annot)
    V=${VARIANTS[$((t - 1))]}
    for x in te nonte; do
        n=$(ls $W/gnomad_${V}_$x/gn.chr*.vcf.gz.tbi 2>/dev/null | wc -l)
        [ "$n" -eq 23 ] || { echo "gnomAD subsets for $V $x: $n of 23"; exit 1; }
    done
    annotate $W/step3_$V/map_grouped.parquet $W/step3_$V/copies_grouped.tsv.gz $W/step1_$V $W/gnomad_${V}_te \
        $W/annot_$V $W/step3_$V/consensus_grouped.fa
    python $S/03_map_to_consensus.py --copies $W/nonte_$V/te_copies.tsv.gz --mode linear --out $W/nonte_$V/map.parquet
    annotate $W/nonte_$V/map.parquet $W/nonte_$V/te_copies.tsv.gz $W/step1_$V $W/gnomad_${V}_nonte $W/annot_nonte_$V ;;
loose_analyze)
    if [ $t -le 2 ]; then V=${VARIANTS[$((t - 1))]}; AT=$W/annot_$V; AN=$W/annot_nonte_$V
    else V=main; AT=$W/annot; AN=$W/annot_nonte; fi
    mkdir -p $P5/loose
    for G in "${G4[@]}" "${GX[@]}" nonTE; do
        T=$(tag "$G"); A=$AT; [ "$G" = nonTE ] && A=$AN
        python $S/32_flank_gc.py --annot $A --subfamilies "$G" --query "$AUTO" --genome-fasta $D/hg38.fa \
            --n-boot 500 --out $P5/loose/${V}_${T}_all
        [ "$V" = main ] && continue   # same-gene results of the main set exist already (patch 4)
        python $S/32_flank_gc.py --annot $A --subfamilies "$G" --query "$AUTO" --genome-fasta $D/hg38.fa \
            --paired-genes --n-boot 500 --out $P5/loose/${V}_${T}_paired
    done
    [ "$V" = main ] || python $S/15_conserved_fragments.py --annot $AT --out $P5/loose/${V}_conserved ;;
scan)
    G=${G4[$(((t - 1) % 4))]}; ORI=$([ $t -le 4 ] && echo sense || echo antisense); T=$(tag "$G")
    mkdir -p $P5/power_scan
    python $S/33_spikein_scan.py --annot $W/annot --subfamilies "$G" --orientation $ORI --scan 20 \
        --reductions 0,0.1,0.25,0.5 --fractions 1,0.5 --reps 10 --n-perm 300 --boot-genes --query "$AUTO" \
        --out $P5/power_scan/${T}_$ORI ;;
usage_sum)
    python $S/34_usage_summary.py --positions $W/splice/res_positions.tsv.gz --usage $RW/junction_usage.tsv.gz \
        --out $P5/usage_summary ;;
boundary_null)
    A=()
    for G in "${G4[@]}" "${GX[@]}"; do T=$(tag "$G"); [ -d $RW/annot_sub/$T ] && A+=(--annot-sub "$T=$RW/annot_sub/$T"); done
    python $S/35_boundary_null.py "${A[@]}" --copies $W/step3/copies_grouped.tsv.gz \
        --splice-sites $W/step1_all/splice_sites.tsv.gz --gtf $GTF --query "$AUTO" --out $P5/boundary_null ;;
age_ca)
    python $S/36_age_copyaware.py --sites $RW/splice_age_sites.tsv.gz --species "$(IFS=,; echo "${SPECIES[*]}")" \
        --out $P5/splice_age_ca ;;
bounds_sub)
    python $S/37_bounds_subsets.py --genes $W/trend/te_vs_nonte_genes.tsv --comp-dir $W/comp --n-boot 2000 \
        --out $P5/bounds_sub ;;
feature_bg)
    python $S/38_feature_local_bg.py --gtf $GTF --gene-set lnc=$W/step1_all/lnc_genes.tsv.gz \
        --gene-set pc=$W/step1_pc/lnc_genes.tsv.gz --rmsk-out $D/hg38.fa.out.gz --out $P5/feature_local ;;
report5)
    O=results/revision5; mkdir -p $O
    for f in $P5/loose/*_summary.tsv $P5/loose/*_bins.tsv $P5/power_scan/*_power.tsv $P5/usage_summary_*.tsv \
             $P5/boundary_null_summary.tsv $P5/splice_age_ca_summary.tsv $P5/splice_age_ca_reclassified.tsv \
             $P5/bounds_sub_bounds.tsv $P5/bounds_sub_orientation.tsv $P5/bounds_sub_functional.tsv \
             $P5/bounds_sub_subsets.tsv $P5/feature_local_local.tsv; do
        [ -s "$f" ] && cp "$f" $O/
    done
    for V in tss0 open; do
        for s in compartments_summary.tsv lnc_genes.tsv.gz; do
            [ -s $W/step1_$V/$s ] && cp $W/step1_$V/$s $O/step1_${V}_$s
        done
    done
    (cd $O && md5sum * > MD5SUMS && ls -l)
    tar czf results/revision5_tables.tgz -C results revision5
    ls -l results/revision5_tables.tgz ;;
*) echo "unknown part: $1"; exit 1 ;;
esac
stamp "done $1 ${t}"
