# Supplementary material

Supplementary Figures S1–S9, Supplementary Tables S1–S9 and Supplementary Text S1–S5 for “Transposable elements supply evolutionarily labile parts to human lncRNA genes with little evidence of purifying selection”. Supplementary Code, intermediate tables and the tests of the pipeline: [Zenodo DOI].

## Supplementary text

### Supplementary Text S1. Why earlier studies saw higher conservation of exonic TE sequence

Earlier genome-wide studies, which reported higher conservation of TE sequence in lncRNA exons (Kapusta et al. 2013; Carlevaro-Fita et al. 2019), compared exonic with intronic TE fragments pooled across genes, without a repeat-free baseline and without the exclusions used here. In our filtered data, such pooled comparisons showed no higher conservation of exonic copies than of intronic copies. To see where the earlier signal could come from, we rebuilt the exonic and intronic compartments without any of our exclusions. Copies containing a stretch highly conserved across mammals then became more frequent in exons than in introns for MIR, L2 and two Alu groups, much as reported earlier, but repeat-free sequence of the same compartments showed a similar excess for every group except AluY, whose excess remained larger, and within the same genes the excess of TE copies was no longer significant (Supplementary Table S3H,I). Such a signal can therefore come from the sequence we excluded, for example sequence shared with other genes or segmental duplications, which may be more conserved in exons whether or not it derives from TEs.

### Supplementary Text S2. The main test in detail

**TE-specific excess and limits.** Within each gene, the TE-specific excess is the exon–intron difference of TE copies minus the exon–intron difference of repeat-free sequence of the same gene (main text, Fig. 3B). It was computed separately for variant frequency, at non-CpG sites, and for conservation, as mean phyloP. For variant frequency we report a limit: the lower end of the 95% interval of the TE-specific reduction, expressed relative to the variant frequency of intronic copies of the group. AluY, MIR3 and L2d had too few exonic copies for useful limits.

**Robustness of the limits.** For AluS, MIR and L2a/b/c the limits stayed at or below 7% in all robustness checks: with sense and antisense copies analyzed separately; in an extended gene set, which removes only the exons of overlapping genes rather than their whole spans; and with exons restricted to the Ensembl canonical transcript of each gene, which removes rarely used alternative exons (Supplementary Table S3A). Of the twelve estimates for single orientations, two had intervals that excluded zero, without correction for multiple testing: antisense AluJ and antisense LTR copies (main text, Fig. 3B).

**Conservation.** The analyses above concern human variants; the following concern conservation across species. At equal GC content, the exon–intron difference in conservation of every TE group, including the LTR and L1 classes, was smaller than that of repeat-free lncRNA sequence across mammals and similar to it across primates (main text, Fig. 3A; Supplementary Table S3C). We also asked the simpler question whether exonic copies are more conserved than intronic copies at all, without subtracting the background of repeat-free sequence. This was the case only for L1 elements and, across primates, for LTR elements. For the TE-specific excess, no TE group was significantly more conserved than repeat-free sequence of the same genes, although for AluS and AluJ these limits are wide compared with the conservation of repeat-free exons themselves (main text, Table 1).

**Matching on flanking GC content.** The GC content of a copy partly reflects its own conservation, because copies that accumulate few substitutions keep more of their original GC. We therefore also matched copies on the GC content of 1 kb of flanking sequence on each side, which changes in the copy itself do not affect. Protein-coding exons remained strongly conserved, no TE group was more conserved in exons than in introns, and repeat-free lncRNA exons were no longer significantly more conserved than their introns (Supplementary Table S3C). Whether lncRNA exons as a whole are conserved therefore depends on how GC content is matched, but the absence of extra conservation in TE copies does not.

### Supplementary Text S3. Sensitivity of the test

**Positive control.** To ask whether human variation registers selection on TE sequence at all, we needed TE sequence that is probably under selection. We took stretches of TE copies that are conserved across species (20-nt windows of single copies with high phyloP; Methods). We chose them by conservation alone, without looking at human variants, so that the choice could not create the result we were testing for. If human variation reflects selection, these stretches should carry fewer variants than the same positions in other copies of the same family. They did: 5–13% fewer, depending on the group (main text, Fig. 3C and Table 1), with a larger share of rare variants in Alu elements; in the smaller MIR3 and L2d groups the difference was not significant. This control pools copies from exons, introns and intergenic sequence; in exonic copies alone the differences were similar in size but mostly not significant (Supplementary Table S4). Part of the difference may reflect lower local mutation rates rather than selection.

**From limits to shares of copies.** The positive control also indicates how strong a realistic protected domain might be: conserved stretches of the main groups carried 5–7% fewer variants. This lets us translate the limits of Supplementary Text S2 into the share of exonic copies that could carry such a domain without being noticed. A domain that covers a fraction *d* of the element and removes a fraction *r* of its variants, present in a fraction φ of exonic copies, lowers the variant frequency of the whole element by φ·*d*·*r* (Methods). Selection as strong as in conserved stretches, acting along the whole element, would be detected only if it were present in more than about half of exonic AluS copies, a third of MIR copies or two thirds of L2a/b/c copies; if part of the difference in conserved stretches is mutational rather than selective, these shares are minimum estimates. A domain of the same strength confined to a fifth of the element would go unnoticed even if every copy carried it. The data set a useful limit only for stronger domains: a domain covering a fifth of the element and removing half of its variants could be present in at most a quarter to a third of exonic AluS, MIR or L2a/b/c copies (main text, Table 1).

**Planted domains.** To test directly how often the analysis finds a domain of known strength, we planted one into the real data: we removed a share of the observed variants from exonic copies and reran the tests (main text, Fig. 3D; Supplementary Fig. S3A). Because detection depends on where the domain lies, we also moved a planted 20-nt domain along the element and resampled genes in every replicate (Supplementary Fig. S3B). A strong domain, with half of its variants removed in half of the exonic copies, was detected in at least 80% of replicates at most positions of AluS and MIR, but at few positions of L2a/b/c unless it was longer or present in all copies (main text, Table 1). A weak domain, with 10% of its variants removed, somewhat more than in the conserved stretches of most groups (5–7%), was detected in about half of the replicates in AluS even when every copy carried it, and in a quarter or fewer elsewhere. The analysis therefore rules out strongly constrained domains shared by many copies, but not weak or rare ones.

### Supplementary Text S4. Searches for selection that the whole-element test could miss

The whole-element test could miss selection confined to part of an element, to a subset of copies or to lncRNAs with particular properties, so we looked for each of these.

**Parts of the element.** No 20-nt window of any group stood out with fewer variants than the rest of its element (Supplementary Fig. S4). Two segments were more conserved in exonic than in intronic copies, although they did not carry fewer variants: a central segment of antisense L2a/b/c copies and the A-rich 3′ end of antisense AluS copies, where alignments to the consensus are least reliable (Supplementary Figs. S4, S5). The Alu segment coincides with the main cluster of splice acceptors in antisense Alu copies (main text, Fig. 4C), but it remained more conserved at bases at least 50 nt from splice sites, so splicing alone does not explain it. The only whole-element deficit of variants, in antisense AluJ, persisted away from splice sites, but the part of it specific to the TE was marginal and disappeared in the extended gene set (Supplementary Table S3A,B).

**Subsets of copies.** Copies containing a highly conserved stretch were no more frequent in exons than in introns (Supplementary Fig. S6).

**Particular lncRNAs.** If TE domains acted only in lncRNAs with particular properties, such as retention in the nucleus (Lubelsky and Ulitsky 2018), the TE-specific excess should change with those properties; it did not follow the stability, nuclear retention, splicing efficiency or expression of the host lncRNA (Supplementary Fig. S7). The most highly or most broadly expressed lncRNAs showed no deficit either, but with only 23–77 genes per group the limits rose to 3–14% (Supplementary Table S3G). In lncRNAs with independent support, either genes with a GTEx expression estimate or intergenic lncRNAs already annotated in GENCODE release 29, the limits were 3–8% and 3.5–7% (Supplementary Table S3G).

These searches thus found two segments with higher conservation and one marginal deficit, but no selection that the whole-element test had missed.

### Supplementary Text S5. Splice sites by junction usage and by origin

**Junction usage.** Selection might act only on splice junctions that are actually used; we therefore split sites into thirds by the number of RNA-seq reads that support their junction in GTEx. In protein-coding genes both measures of constraint rose steeply with usage, so the split can reveal it (main text, Fig. 7B). In lncRNAs, sites outside repeats became more conserved with usage but did not carry fewer variants, and TE-derived sites showed no consistent trend. Two caveats apply. About three quarters of lncRNA junctions are absent from the GTEx junction table and were analyzed separately. And even the most used third of lncRNA junctions had fewer reads than most protein-coding junctions (Supplementary Table S7B), so at comparable usage the variant data alone could not distinguish TE-derived from protein-coding sites.

**Created and preformed sites.** Created sites carried about 14% fewer variants than the matched positions in unused copies, whereas preformed sites did not (main text, Fig. 7C). Two observations argue for caution: the reduction came without a shift toward rare variants and with lower, not higher, conservation, and it weakened when CpG status, which strongly affects mutation rate, was also matched (Supplementary Table S7). One observation supports it: among the few created sites in copies shared with mouse, the dinucleotide was present in mouse more often than in the controls (Supplementary Table S6B). Hence weak selection on a subset of old created sites cannot be excluded, but most TE-derived splice sites show no sign of being maintained by selection.

## Supplementary figures

[image: Four line plots of simulated exon-intron variant profiles along the Alu consensus, with a dip only at the planted segment in sense copies]

**Supplementary Figure S1.** The pipeline recovers a planted signal, and nothing else, in simulated data. An artificial genome in the formats of the real input files carried a planted signal at consensus positions 101–130 of exonic AluS copies in sense orientation (half the variant rate, more singletons and higher phyloP; only the variant measure is shown), including a related subfamily with insertions and deletions that is placed on the shared consensus. (A) Sense copies: exon versus intron log odds ratio of a site being variable along the consensus; points, windows with q < 0.1. (B) Antisense copies. (C) Pseudo-exons placed at random in introns, compared with introns. (D) With an additional uniform 25% deficit of variants along the whole of the exonic sense copies, the raw profile (gray) shifts along its whole length, whereas the centered profile (vermilion) isolates only the planted segment. Shading, planted segment; bands, 95% intervals. The planted effects are large; sensitivity to weaker and rarer effects in the real data is shown in Fig. 3D and Supplementary Fig. S3.

[image: Two line plots of mean phyloP by distance from the transcription start for exon and intron types]

**Supplementary Figure S2.** Including the regions around transcription starts does not make lncRNA exons more conserved than introns. Mean phyloP across 447 mammals of repeat-free sequence of reference transcripts (Ensembl canonical, otherwise the longest transcript of the gene), by region and by distance from the transcription start (exons measured along the spliced transcript, introns along the genome), for 4,000 randomly chosen analyzed lncRNAs (A) and 6,000 randomly chosen lncRNAs expressed at 1 TPM or more in at least one GTEx tissue, genes absent from GTEx excluded (B). Shading, the first 500 nt, removed from exons in the main analysis (1 kb was removed from introns). Bins with fewer than 20,000 bases are omitted.

[image: Grid of heatmaps of detection rates of planted constraint by segment length, variant reduction and share of copies for four TE groups]

[image: Supplementary Figure S3B: detection of a planted 20-nt segment along the consensus, with genes resampled, for four TE groups and both orientations]

**Supplementary Figure S3.** Detection of planted constraint in the real data. (A) All settings, with the segment at one position. Variants were removed by binomial thinning in a 20-, 40- or 60-nt consensus segment or along the whole element, in 25%, 50% or 100% of the exonic copies of one orientation, and the window and whole-element tests were rerun on the same-gene design (20 replicates, 1,000 permutations each). Replicates differ only in which copies and variants were thinned, not in the data themselves, so detection is close to 0 or 1 and can depend on the position of the segment. Segments start at consensus position 101 (AluS, AluJ), 81 (MIR) or 201 (L2a/b/c). Cells give the share of replicates in which the window test (a window of the segment with q < 0.1 and fewer variants in exons) or, for the whole element, the whole-element test (P < 0.05, fewer variants) detected the planted constraint. The first column of each block, with no variants removed, shows the result on the unmodified data, which is the same in every replicate and is not a false-positive rate; the whole-element test flags antisense AluJ there because its deficit of variants is present in the real data. (B) A planted 20-nt segment moved along the element in 20-nt steps, with genes resampled with replacement in each of 20 replicates (1,000 permutations each): share of replicates in which the window test detected it, by the share of variants removed (shade) and of exonic copies carrying it (solid, all; dashed, half); dotted, nothing planted (detections inside the segment on resampled data); gray line, 80%. With half of the variants removed in half of the exonic copies, at least 80% of replicates detected the segment at 87–94% of positions in AluS and antisense AluJ, 57–71% in MIR, 14–19% in L2a/b/c and none in sense AluJ; with 10% removed in all exonic copies, median detection was 40–48% in AluS and 0–25% in the other groups.

[image: Grid of exon-intron profiles along the consensus for six TE groups, variants and phyloP, sense and antisense]

**Supplementary Figure S4.** Exon–intron profiles along the consensus for the Alu, MIR and L2 groups (AluS, AluJ, AluY, MIR, MIR3 and L2a/b/c), within the same genes. Exon minus intron differences in 20-nt windows for variants (Mantel–Haenszel log odds ratio of a site being variable, stratified by trinucleotide context, CpG status and original subfamily; below 0, fewer variants in exons) and for conservation across 447 mammals (stratified phyloP difference, not matched for GC content; above 0, more conserved in exons), in sense and antisense copies; bands, 95% intervals from 200 resamples of genes. Points, windows differing from zero (q < 0.1). Gray, segments of at least two adjacent windows differing from the rest of the element (centered profiles, q < 0.1): positions 261–300 of antisense AluS and 221–260 of antisense L2a/b/c, both more conserved in exons. No window differed from the rest of the element in variant frequency. Across all genes (not shown), the segments flagged in conservation were 1–40 of sense AluJ, 181–220 of sense MIR, 341–400 and 441–480 of antisense and 421–460 of sense L2a/b/c, all less conserved in exons, and again 221–260 of antisense L2a/b/c; none of the more conserved segments recurred as a segment of two windows in the lincRNA gene set (one window of the L2a/b/c segment did).

[image: Two dot-and-whisker panels of whole-element variant and consensus-match differences in subsets of copies]

**Supplementary Figure S5.** Whole-element estimates in subsets of copies and bases, within the same genes. (A) Exon versus intron log odds ratio of a site being variable for all copies, for exonic copies embedded in an exon without containing an exon boundary, for bases at least 50 nt from splice sites, and for complete Alu dimers (with their intronic counterparts). (B) Exon minus intron difference in the share of non-CpG bases identical to the reference consensus (consensus-match measure), which is higher in copies with fewer substitutions since insertion. Bars, 95% intervals from 200 resamples of genes. The deficit of variants in antisense AluJ copies persisted away from splice sites. In centered profiles of these subsets (not shown), no window differed from the rest of the element in variant frequency (q < 0.1); in conservation, the excess at the 3′ end of antisense Alu copies reappeared away from splice sites (AluS, windows 221–260 and 281–300; AluJ, 221–300, with 30–37 exonic units per window beyond position 240) and in AluS dimers (one window, 261–280). Antisense AluS and sense L2a/b/c copies in exons matched their consensus slightly more often than intronic copies (by 0.010 and 0.016; permutation q = 0.15 and 0.07, and q = 0.03 for antisense AluS away from splice sites); sense AluS copies tended the other way (−0.012, interval including zero).

[image: Two dot panels of the share of exonic and intronic MIR and L2 units with conserved stretches, all genes and same genes]

**Supplementary Figure S6.** Conserved stretches are not more frequent in exonic than in intronic copies. Share of exonic (circles) and intronic (squares) units of each MIR and L2 subfamily that contain a 20-nt stretch with mean phyloP across 447 mammals above 1.5, across all genes (A) and within genes that contain both exonic and intronic units of the family (B). Numbers give conserved units and all units, exonic and intronic, and the Mantel–Haenszel odds ratio of an exonic unit being conserved, stratified by unit length (quartiles) and GC content (terciles), with 95% intervals from 500 resamples of genes; subfamilies without conserved exonic units cannot be tested. Among Alu copies, no exonic and one intronic unit contained a 20-nt stretch with mean primate phyloP above 1.0.

[image: Two dot-and-whisker panels of Spearman correlations of the TE-specific excess with host-gene RNA properties]

**Supplementary Figure S7.** The TE-specific excess does not follow RNA properties of the host genes. Spearman correlation across genes between the gene-level TE-specific excess (the TE exon–intron difference minus that of repeat-free sequence of the same gene; both orientations pooled) in variant frequency (A) or mammalian phyloP (B) and properties of the host lncRNA: stability, nuclear retention, splicing efficiency and synthesis rate from ENCODE HepG2 and K562 data (Supplementary Table S8), and expression and tissue specificity from GTEx v10. Bars, 95% intervals from 500 resamples of genes. Of the 48 correlations, 2 had P < 0.05 and none had a false discovery rate below 0.1. Genes per correlation: 34–65 for nuclear retention, 51–111 for stability and synthesis, 87–182 for splicing efficiency and 121–305 for GTEx expression and tissue specificity. Panel B uses conservation across mammals for all groups; with conservation across primates, the informative measure for Alu, no correlation of AluS or AluJ had P < 0.05 either.

[image: Supplementary Figure S8: orientation of LTRs at lncRNA transcription starts against expectation, and positions of starts along the consensus of six LTR groups]

**Supplementary Figure S8.** lncRNA transcription starts in LTR elements. (A) Share of starts inside RepeatMasker LTR records that lie in records oriented with the transcript (number of starts in parentheses; lines, 95% binomial intervals), against the share expected from the orientation of LTR sequence in the exons of the same genes (vertical ticks). (B) Starts in sense LTRs along the consensus (10-nt bins) for each named LTR group (records named THE1, MST, MLT1, MLT2 or LTR12 by prefix; LTR7: LTR7, LTR7B, LTR7C and LTR7Y, the LTRs of HERVH, unlike the broader LTR7 name pattern of Supplementary Table S2); shaded, the most used 50-nt window, with its share of the group's sense starts and the share expected for the most used window if each start were placed at random in its own copy (1,000 permutations; P < 0.001 for every group). Axes end at the 98th percentile of positions. “Other LTR” pools subfamilies with different consensus sequences and is shown in A only. Within the main subfamily of each group (at least 20 sense starts, internal-region records excluded), the most used window held 71% of sense starts in THE1B, 63% in MSTA and 57% in LTR7, against 18%, 17% and 16% expected (P < 0.001 each).

[image: Mirrored histograms of lncRNA splice acceptors and donors along the MIR and L2 consensus]

**Supplementary Figure S9.** lncRNAs use recurrent positions of MIR and L2 elements as splice sites. Numbers of lncRNA splice acceptors (up) and donors (down) along the MIR consensus (A) and the 3′ part of the L2 consensus (B), in copies antisense and sense to the transcript. Consensus positions come from linear interpolation of the RepeatMasker coordinates of each subfamily, so peaks of related subfamilies may be offset by a few positions.

## Supplementary tables

### Supplementary Table S1. Selection of genes and exonic sequence

Autosomal lncRNA genes of GENCODE release 50 and the exonic sequence retained after each filter used to build the compartments. Removed alone: sequence removed by the filter applied on its own to the starting set. Left: sequence, and genes with any exonic sequence, remaining after the filters are applied in the order shown.

Genes: of 33,839 autosomal lncRNA genes, 10,626 were removed because an exon overlapped an exon of a protein-coding gene or pseudogene and 4,603 because they were below 1 TPM in all GTEx tissues; 18,610 were kept, including genes absent from GTEx.

{{TABLE:S1}}

### Supplementary Table S2. LTR and L1 groups and their mapping

Subfamilies of LTR retrotransposons (THE1, MST, MLT1, MLT2, LTR12, and LTR7, whose name pattern also matched the unrelated LTR71–LTR79 subfamilies) and of old L1 elements (L1M) selected by name pattern. Whole-class tests (LTR, L1) pool all listed subfamilies, with the original subfamily as a stratum. For profiles along the consensus, subfamilies at least 75% identical to the reference consensus of their family pattern were placed on that consensus by global alignment; the others kept their own consensus. Copies: copies selected in the compartments (intronic and intergenic copies subsampled per subfamily). Identity: identity of the subfamily consensus to the reference consensus; mapped: share of the subfamily consensus aligned to the reference. L1 copies are mapped with the 3′ RepeatMasker model; most L1M subfamilies have only 5′ and ORF2 models, so the whole-class L1 tests also used all copies with linearly interpolated consensus coordinates (Methods).

{{TABLE:S2}}

### Supplementary Table S3. Sensitivity analyses, gene composition and GC content

**A. Limits on TE-specific constraint.** Largest TE-specific reduction in the frequency of variable callable non-CpG sites consistent with the data (lower 95% bound of the gene-level excess over repeat-free sequence of the same genes, relative to the intronic frequency of such sites; 2,000 resamples of genes), with the number of genes in parentheses. Main: both orientations, unweighted mean over genes (Table 1). Weighted: genes weighted by n_exon·n_intron/(n_exon + n_intron) units. Embedded: exonic copies that lie inside an exon without containing an exon boundary. Extended: only the exons of other genes removed instead of their whole spans. 100-nt margin: 100 nt instead of 500 nt removed around transcription starts in exons. Of the twelve orientation-specific estimates, two had 95% intervals excluding zero: antisense AluJ (mean −5.3%, −10.5% to −0.1%) and antisense LTR (−3.9%, −7.1% to −0.7%). In the extended gene set the antisense AluJ estimate was +0.6% (−4.3% to +5.7%). Canonical-transcript exons: exonic compartments restricted to the exons of the Ensembl canonical transcript of each gene, so that rarely used alternative exons are removed (intronic compartments unchanged).

{{TABLE:S3A}}

**B. Whole-element exon–intron tests in five analyses.** Main: genes that contain both exonic and intronic copies, with the main exclusions. Within gene: the gene added as a stratum, so exonic copies are compared only with intronic copies of the same gene, with labels permuted within genes. 100-nt margin and extended: as in part A. Variants: Mantel–Haenszel log odds ratio of a site being variable (below 0, fewer variants in exons). phyloP: stratified difference in mean phyloP across 447 mammals (above 0, more conserved in exons); GC-matched: both orientations, intronic means weighted to the exonic GC distribution. Centered peaks: segments of at least two adjacent 20-nt windows differing from the rest of the element (q < 0.1), all measures. 95% intervals from 200 resamples of genes. Apart from the deficit of variants in antisense AluJ copies (main, 100-nt and extended analyses), the only interval excluding zero in the direction of constraint is phyloP of antisense AluS in the extended within-gene analysis, one of about 80 comparisons.

{{TABLE:S3B}}

**C. GC matching on the copy or on its flanks.** Exon minus intron mean phyloP within the same genes: unmatched; matched on the GC content of the unit itself (quintiles; as in Fig. 3A); matched on the GC content of 1 kb of flanking sequence on each side, which substitutions inside the unit do not change (quintiles); and matched on both (terciles of each). 95% intervals from 200 resamples of genes. Protein-coding exons stay strongly constrained under every matching; no TE group is significantly more conserved in exons than in introns when matched on flank GC; the exonic excess of repeat-free lncRNA sequence depends on matching on the GC content of the sequence itself.

{{TABLE:S3C}}

**D. Exonic versus intronic units across all genes.** Whole-element exon versus intron comparison without restricting to genes with both kinds of units (Fig. 2A, open circles): Mantel–Haenszel log odds ratio of a site being variable (stratified by trinucleotide context, CpG status and original subfamily; below 0, fewer variants in exons) and difference in the share of singletons among variable sites (above 0, rarer variants in exons, as expected under purifying selection). 95% intervals from 200 resamples of genes. Repeat-free: 100-nt windows of the same compartments (20,000 sampled per compartment), analyzed in the same way.

{{TABLE:S3D}}

**E. Genes that drop out of the same-gene comparison.** Genes with exonic units of a group and orientation, split by whether they also contain intronic units of that group and orientation (both) or not (exonic only). Medians over genes; splicing efficiency and synthesis rate from the ENCODE-based gene properties (Supplementary Table S8; log scale, higher = more efficiently spliced or faster synthesized); maximum GTEx v10 expression in TPM.

{{TABLE:S3E}}

**F. GC content of exonic and intronic units within the same genes.** Share of units in the two GC-richest quintiles (quintiles over exonic and intronic units of the group together, as used for matching in part C and Fig. 3A).

{{TABLE:S3F}}

**G. Limits in highly expressed, broadly expressed and independently supported lncRNAs.** Gene-level TE-specific excess (as in part A, both orientations) in subsets of the genes of the gene-level test: the quarter and the tenth with the highest maximum GTEx v10 expression, and the quarter with the broadest expression (lowest τ), with quantiles over the genes that have the respective value (GTEx expression or a defined τ). Mean change in variant frequency relative to the intronic frequency (%), 95% interval from 2,000 resamples of genes, number of genes, and the largest reduction consistent with the data. Subsets with fewer than 15 genes per group were not tested (–). Thresholds: maximum expression at least 6.7 TPM (top quarter, 215 genes) or 14.8 TPM (top tenth, 86 genes); τ at most 0.83 (broadest quarter, 183 genes); of 1,755 genes with TE and repeat-free units. Independent support: genes with a GTEx v10 expression estimate (858 genes) and genes already classified as lincRNAs in GENCODE release 29 (595 genes, all with a GTEx estimate).

{{TABLE:S3G}}

**H. Exon–intron conservation without our exclusions.** The comparison of earlier studies (exonic versus intronic units, mean phyloP across 447 mammals) in three compartment sets: main (the exclusions of Supplementary Table S1); no start/end margins (no margins around transcription starts and transcript ends); no exclusions (also no margins around exon boundaries, no overlap filter, and sequence shared with other genes and segmental duplications kept). All genes: units pooled across genes; same genes: genes with exonic and intronic units of the group. Exon minus intron difference, unmatched and matched jointly on the GC content of the unit and of its 1-kb flanks (as in part C), with 95% intervals from 500 resamples of genes (200 for the main same-gene set, from part C). Exonic / intronic units are given for each set. Repeat-free: 100-nt windows of the same compartments (20,000 sampled per compartment). –, not computed.

{{TABLE:S3H}}

**I. Conserved stretches in exonic and intronic units without our exclusions.** Units (copies or their parts in one compartment; 100-nt windows for repeat-free sequence) with a 20-nt stretch of mean phyloP above 1.5 across 447 mammals, in the three compartment sets of part H. Exonic / intronic units with a conserved stretch (of all units), and the exon versus intron Mantel–Haenszel odds ratio stratified by original subfamily and unit-length quartile (longer units contain more windows); ratio to repeat-free: the odds ratio of the TE group divided by that of repeat-free sequence of the same compartments (TE-specific excess). 95% intervals from 1,000 resamples of genes, drawn jointly for TE and repeat-free units. –, no conserved exonic unit or too few to estimate; intervals are unstable where few units are conserved.

All genes.

{{TABLE:S3I_all}}

Same genes.

{{TABLE:S3I_same}}

### Supplementary Table S4. Positive control by group, threshold and compartment

Human variants in stretches of TE copies conserved across species. Conserved stretches are 20-nt windows of single copies selected by phyloP alone: across primates for Alu, across 447 mammals for the older groups, with a fixed threshold on the window mean or the top 1% or 5% of windows of the group. Units: copies or their parts in the compartment; share: units with at least one conserved stretch. Fewer variants: 100·(1 − OR) for the Mantel–Haenszel odds ratio of a site being variable, stratified by subfamily, consensus bin, context and CpG status, (A) against the same consensus positions in copies without a stretch and (B) against the other bases of the same copies; singleton difference: share of singletons among variable sites, stretch minus comparison A. 95% intervals from 200 resamples of genes. The main text and Fig. 3C use the top 1% of windows for Alu and a mean above 1.0 for the older groups, with all compartments pooled. Rows with fewer than 10 units carrying a stretch are omitted (for AluS with a primate threshold of 1.0, almost no copy reaches it).

{{TABLE:S4}}

### Supplementary Table S5. TEs at transcription starts, splice sites and exon boundaries

The 25 subfamilies of LTR retrotransposons (RepeatMasker class LTR, including internal regions) containing the most lncRNA transcription starts. Starts: transcription starts of lncRNA transcripts in a copy of the subfamily (number and share of all lncRNA starts); exonic sequence: share of lncRNA exonic sequence in the subfamily, as the expectation; observed/expected: share of starts divided by share of exonic sequence; protein-coding: share of protein-coding transcription starts in the subfamily; lncRNA/protein-coding: ratio of the two shares of starts. LTR elements together contain 17% of lncRNA starts and about 1% of protein-coding starts (Fig. 4A).

{{TABLE:S5A}}

**B. TE shares at gene features and in the sequence next to them.** Share of features lying in any TE, in Alu and in LTR elements, compared with positions 20–100 and 100–500 nt away on the exonic side of the feature (downstream of starts and acceptors, upstream of donors and ends, along the transcript) and 20–100 nt away on the other side (upstream of starts, on the intronic side of splice sites, downstream of ends), two positions per feature, side and distance bin. Exonic-side positions falling outside an exon were not used. Ratio = feature / exonic side 100–500 nt; positions 20–100 nt away often lie in the same TE copy as the feature. Alu and LTR columns: feature / exonic side 20–100 nt. All transcripts of the analyzed lncRNA genes and of protein-coding genes; RepeatMasker classes as in Fig. 4A.

{{TABLE:S5B}}

**C. Exonic copies that contain an exon boundary, observed and expected by chance.** Observed: share of exonic copies (units of the exon-boundary analysis, Fig. 4B) whose span contains a splice site of the host gene. Expected: the same share if each copy, with its length, were placed uniformly at the positions where it overlaps the same exon (exons merged per gene). 95% intervals from 1,000 resamples of genes.

{{TABLE:S5C}}

### Supplementary Table S6. Splice-site age and presence relative to controls

Canonical GT–AG splice sites of the analyzed genes. Preformed: share of TE-derived sites whose dinucleotide was already present in the consensus of the original subfamily. Presence: share of sites whose orthologous dinucleotide is GT or AG, among sites aligned to each species by the UCSC chains (%). After the macaque split: share of all sites present at most in apes (human-specific or shared with chimpanzee, gorilla or orangutan only; the age of at most 0.4% of the sites of a subset is unresolved). 95% intervals of the presence shares, from resamples of genes, are within ±3 points for the full sets (Fig. 6).

{{TABLE:S6A}}

**B. Presence of splice sites and of their controls at the orthologous position.** Sites with at least one aligned control; controls are the same dinucleotide at the same consensus positions of unused copies (TE-derived sites) or GT or AG inside the same intron (other sites). Values: site minus control, in percentage points, with 95% intervals from 1,000 resamples of genes (positive = sites kept more often than their controls), and the share of controls present (in parentheses). For TE-derived sites, 'copy shared' restricts sites and controls to copies that are themselves present in that species, which removes non-orthologous alignments (Alu elements are absent from mouse).

{{TABLE:S6B}}

**C. Ages of TE-derived splice sites when the copy itself must be present.** Site rule (as in the first part of this table and Fig. 6C): a site is present in another species when the orthologous dinucleotide is GT or AG. Copy rule: in addition, both ends of the TE copy (20 nt inside it) must lift to the same chromosome and strand of that species at a distance of 0.5–2 times the copy length; sites in copies that were tested and failed are not counted as present there, whereas copies too short to test (under 60 nt) keep the site call. The copy rule removes non-orthologous alignments, such as Alu-derived sites 'present' in mouse (Alu elements are primate-specific), but also copies broken by indels or alignment gaps, so it makes old sites younger. Other lncRNA and protein-coding sites are unchanged (8% and 2% after the macaque split).

{{TABLE:S6C}}

In mouse, the copy rule set 26 of 35 Alu-derived, 112 of 353 MIR-derived and 72 of 180 L2-derived sites with an aligned GT or AG to not present; in macaque, 217 of 4,048 sites.

### Supplementary Table S7. Splice-site constraint by usage, age and origin

Differences between the two bases of canonical splice-site dinucleotides and their controls: protein-coding sites and lncRNA sites outside repeats against intronic bases 60 nt away; TE-derived sites against the same consensus positions and bases in copies hosting no splice site, from any compartment or from introns only, with the strata of context, CpG status and subfamily, or additionally split by CpG status. Usage: tertiles of GTEx v10 junction reads across samples within each set of sites; preformed and created: dinucleotide already present in the consensus or created by mutation after insertion. Bases: site bases in the comparison. Variants: Mantel–Haenszel log odds ratio of a base being variable (below 0, fewer variants at sites); singletons: difference in the share of singletons among variable bases; phyloP: difference in mean phyloP across 447 mammals. 95% intervals from 1,000 resamples of genes, except the rows “all (Fig. 7A)”, the pooled estimates of Fig. 7A for sites outside TEs from the main splice-site analysis (200 resamples of genes).

{{TABLE:S7A}}

**B. What the usage strata mean.** Sites whose junction is absent from the GTEx v10 junction table (or present without reads) form the stratum without usage estimate; the others were split into tertiles of the reads of their junction summed over all GTEx samples. Reads: minimum–maximum (median) per tertile; samples: median number of samples with at least one read.

{{TABLE:S7B}}

### Supplementary Table S8. ENCODE experiments used for RNA properties

ENCODE experiments (untreated or control HepG2 and K562 cells) from which the per-gene RNA properties in Supplementary Fig. S7 were computed with GENCODE v29 gene models; the per-gene values are deposited with the data (Data Availability). Values were averaged over the two cell lines. Bru-seq and BruChase-seq follow Paulsen et al. (2014). Other ENCODE experiments of these cell lines (2-h chase, chromatin and nucleoplasm fractions, total-RNA libraries, other fractions, eCLIP) are omitted. Expression and tissue specificity (tau) in Supplementary Fig. S7 come from GTEx v10, not from ENCODE (Methods).

{{TABLE:S8}}
<!-- keep:S8 -->
| Property | Assay | Cell line | ENCODE experiments |
| --- | --- | --- | --- |
| Synthesis | Bru-seq | HepG2 | ENCSR974AQD |
| Synthesis | Bru-seq | K562 | ENCSR729WFH |
| Stability (6 h chase relative to Bru-seq) | BruChase-seq, 6 h | HepG2 | ENCSR135DZR |
| Stability (6 h chase relative to Bru-seq) | BruChase-seq, 6 h | K562 | ENCSR762OPQ |
| Nucleus/cytoplasm ratio | poly(A) RNA-seq, nuclear fraction | HepG2 | ENCSR000CPC, ENCSR058OSL |
| Nucleus/cytoplasm ratio | poly(A) RNA-seq, cytosolic fraction | HepG2 | ENCSR000CPF, ENCSR931WGT |
| Nucleus/cytoplasm ratio | poly(A) RNA-seq, nuclear fraction | K562 | ENCSR000CPS, ENCSR530NHO |
| Nucleus/cytoplasm ratio | poly(A) RNA-seq, cytosolic fraction | K562 | ENCSR000BZV, ENCSR000COK, ENCSR384ZXD |
| Splicing efficiency | RNA-seq of control samples (signal tracks) | HepG2 | ENCSR225PRV, ENCSR237YZT, ENCSR452TFV, ENCSR704VPL |
| Splicing efficiency | RNA-seq of control samples (signal tracks) | K562 | ENCSR032YMP, ENCSR164MUK, ENCSR174OYC, ENCSR260BQC, ENCSR271ZAB, ENCSR292PXV, ENCSR344XID, ENCSR854OFV, ENCSR913CAE |
<!-- /keep -->

### Supplementary Table S9. Data sources and versions

All files used by the pipeline, with the download dates; md5 checksums of the phyloP files, liftOver chains and 2bit genomes were checked against the md5sum.txt lists of the same UCSC directories. gnomAD genome VCFs were not downloaded in full; the pipeline streamed the regions of the analyzed copies, windows and splice sites.

{{TABLE:S9}}
<!-- keep:S9 -->
| Data | Release | Files | Source (https://) | Downloaded |
| --- | --- | --- | --- | --- |
| Gene annotation | GENCODE release 50 (Ensembl 116), comprehensive annotation | gencode.v50.annotation.gtf.gz | ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/ | 2026-10-04 |
| Gene expression | GTEx v10, median gene TPM per tissue (RNA-SeQC 2.4.2) | GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz | storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/ | 2026-10-04 |
| Junction usage | GTEx v10, junction read counts (STAR 2.7.10a) | GTEx_Analysis_v10_STARv2.7.10a_junctions.gct.gz | storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/ | 2026-10-25 |
| Sample attributes | GTEx v10 sample attributes | GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt | storage.googleapis.com/adult-gtex/annotations/v10/metadata-files/ | 2026-10-25 |
| Genome | GRCh38 (hg38), UCSC | hg38.fa.gz, hg38.chrom.sizes | hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/ | 2026-10-04 |
| Segmental duplications | UCSC genomicSuperDups | genomicSuperDups.txt.gz | hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ | 2026-10-04 |
| Repeats | RepeatMasker open-4.0.6 with Dfam 2.0, annotation and alignments of the same run | hg38.fa.out.gz, hg38.fa.align.gz | www.repeatmasker.org/genomes/hg38/rmsk4.0.6_dfam2.0/ | 2026-10-04 |
| Conservation | phyloP on the 447-species Cactus alignment and its primate subset | hg38.phyloP447way.bw, hg38.phyloP447wayPrimates.bw | hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP447way/ | 2026-10-04 |
| Population variation | gnomAD v4.1 genomes (76,215 individuals), PASS single-nucleotide variants | gnomad.genomes.v4.1.sites.chr*.vcf.bgz (regions streamed) | gnomad-public-us-east-1.s3.amazonaws.com/release/4.1/vcf/genomes/ | 2026-10-04 |
| Allele numbers | gnomAD v4.1 genomes, allele numbers at all sites | gnomad.genomes.v4.1.allele_number_all_sites.tsv.bgz | gnomad-public-us-east-1.s3.amazonaws.com/release/4.1/tsv/genomes/ | 2026-10-04 |
| Orthologous positions | UCSC liftOver chains hg38 → panTro6, gorGor6, ponAbe3, rheMac10, calJac4, mm39 | hg38To{Species}.over.chain.gz | hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/ | 2026-10-25 |
| Other genomes | panTro6, gorGor6, ponAbe3, rheMac10, calJac4, mm39 (UCSC) | {species}.2bit | hgdownload.soe.ucsc.edu/goldenPath/{species}/bigZips/ | 2026-10-25 |
| RNA properties | ENCODE HepG2 and K562 experiments (Supplementary Table S8), per-gene values computed with GENCODE v29 gene models | gene_properties_controls_v29.tsv.gz | deposited with the data of this study (Zenodo) | – |
<!-- /keep -->

## Supplementary references

Carlevaro-Fita J, Polidori T, Das M, Navarro C, Zoller TI, Johnson R. 2019. Ancient exapted transposable elements promote nuclear enrichment of human long noncoding RNAs. Genome Res. 29:208–222. https://doi.org/10.1101/gr.229922.117

Kapusta A, Kronenberg Z, Lynch VJ, Zhuo X, Ramsay L, Bourque G, Yandell M, Feschotte C. 2013. Transposable elements are major contributors to the origin, diversification, and regulation of vertebrate long noncoding RNAs. PLoS Genet. 9:e1003470. https://doi.org/10.1371/journal.pgen.1003470

Lubelsky Y, Ulitsky I. 2018. Sequences enriched in Alu repeats drive nuclear localization of long RNAs in human cells. Nature. 555:107–111. https://doi.org/10.1038/nature25757

Paulsen MT, Veloso A, Prasad J, Bedi K, Ljungman EA, Magnuson B, Wilson TE, Ljungman M. 2014. Use of Bru-Seq and BruChase-Seq for genome-wide assessment of the synthesis and stability of RNA. Methods. 67:45–54. https://doi.org/10.1016/j.ymeth.2013.08.015
