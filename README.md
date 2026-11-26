# te-lncrna-selection

Code for the analyses in:

> Szcześniak MW, Bocheńska M, Wanowska E. Transposable elements supply evolutionarily labile parts to human lncRNA genes with little evidence of purifying selection. Manuscript submitted.

The pipeline compares human population variation (gnomAD v4.1) and evolutionary conservation (phyloP, 447 mammals and the primate subset) between exonic and intronic copies of transposable elements (TEs) in human lncRNA genes, within the same genes and against repeat-free sequence of those genes. It also covers positive controls, planted-constraint (spike-in) tests, TE-derived splice sites and their age, and LTR elements at lncRNA transcription starts.

Intermediate tables and results are deposited at Zenodo: https://doi.org/10.5281/zenodo.24560064 (all versions; version 1.0, with the code at commit 5b30018: https://doi.org/10.5281/zenodo.24560065).

## Contents

| Directory | Contents |
| --- | --- |
| `pipeline/` | Slurm stage scripts (`s*.sbatch`), the scripts that submit them as job chains (`submit_*.sh`) and the shared settings (`config.sh`, `config_rev.sh`) |
| `scripts/` | Analysis steps (Python), numbered in the order they were added; `ridl_utils.py` and `ridl_stats.py` hold shared functions |
| `tests/` | Tests of the revision steps on simulated data with planted signals |
| `figures/` | Figure scripts of the first full run (stage `s10`); kept so that stage still runs |
| `manuscript/` | Scripts for the figures and tables of the article and its supplement (stage `s16`), the supplement template and the fonts used in the figures |

## Requirements

- Linux with the Slurm workload manager
- Python 3.12 with NumPy 2.5, pandas 3.0, SciPy 1.18, pyarrow 25, pysam 0.24, pyBigWig 0.3.26, py2bit 1.0.1, pyliftover 0.4.1 and matplotlib 3.11 (versions used for the article: Python 3.12.14, NumPy 2.5.3, pandas 3.0.6, SciPy 1.18.1, pyarrow 25.0.0, pysam 0.24.1, matplotlib 3.11.2)
- bcftools and samtools 1.24

## Input data

All inputs are public (Supplementary Table S9 of the article lists files, versions and download dates): GENCODE release 50, GTEx v10 gene expression and junction read counts, GRCh38 (UCSC), UCSC segmental duplications, RepeatMasker open-4.0.6 with Dfam 2.0 (annotation and alignments of the same run), phyloP on the 447-species Cactus alignment, gnomAD v4.1 genomes (variants and allele numbers), UCSC liftOver chains and genomes of six species. Stage `s00` downloads them into `data/`; gnomAD variants are streamed for the needed regions only. Per-gene RNA properties computed from ENCODE HepG2 and K562 data (Supplementary Table S8) are part of the Zenodo deposit; place `gene_properties_controls_v29.tsv.gz` in `data/` or set `RNA_PROPS`.

## Setup

1. Clone the repository; the clone is the project directory (`PROJECT`).
2. Create a conda environment with the packages above.
3. Copy `pipeline/local.sh.example` to `pipeline/local.sh` and set the conda paths there (leave `CONDA_SH` empty if the environment is already active). `local.sh` is not tracked by git.
4. Create `data/`, or link it to an existing download.

## Running

Each `submit_*.sh` script submits a chain of Slurm jobs and accepts `DRY=1` (print the `sbatch` commands only) and `SKIP="<keys>"` (skip finished stages); see the header of each script. The stages were added in successive rounds of analysis ("revision" and "patch" stages) and are kept as they were run. Run them in this order:

| Order | Command | Stages | What it does |
| --- | --- | --- | --- |
| 1 | `DOWNLOAD=1 bash pipeline/submit_all.sh` | `s00`–`s10` | Main analysis: compartments and TE copies, mapping to consensus, gnomAD and phyloP annotation, exon–intron profiles, summaries, splice sites, report |
| 2 | `bash pipeline/submit_revision.sh` | `s00r`, `s03r`, `s11`, `s12` | Further downloads, subsets, splice-site age and junction usage, LTR and L1 elements |
| 3 | `export PIPE=$PWD/pipeline; sbatch --array=1-8 pipeline/s13_flank_gc.sbatch` | `s13` | GC matching on the flanks of each unit |
| 4 | `bash pipeline/submit_patch5.sh` | `s14` | Compartments without exclusions, bounds and further controls |
| 5 | `bash pipeline/submit_patch5b.sh` | `s15` | Planted constraint along the element with genes resampled |
| 6 | `bash pipeline/submit_patch9.sh` | `s17` | Candidate loci for the example figure |
| 7 | `bash pipeline/submit_patch6.sh` | `s16` | Last analyses, then all figures and tables of the article and supplement |

`ONLY_RENDER=1 bash pipeline/submit_patch6.sh` redraws the figures and tables from existing results.

Outputs are written to `work/` (intermediate files), `results/` (result tables) and `logs/`, none of which are tracked.

## License

The code is released under the MIT License (see `LICENSE`). The Liberation Sans fonts in `manuscript/fonts/` are distributed under the SIL Open Font License 1.1 (see `manuscript/fonts/LICENSE_OFL-1.1.txt`).
