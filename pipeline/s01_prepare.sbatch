#!/bin/bash
#SBATCH --job-name=ridl_s01_prepare
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x_%j.log
# Compartments (main set, strict expression set, protein-coding control), TE copies, ordinary-sequence
# windows, exonic attrition (Table S7) and per-gene RNA properties (Table S8, tau from GTEx v10).
source "$PIPE/config.sh"

stamp "step 1: compartments"
python $S/01_build_compartments.py $STEP1_COMMON --keep-unmatched-gtex --out-dir $W/step1_all   # main set
python $S/01_build_compartments.py $STEP1_COMMON --out-dir $W/step1          # GTEx-matched genes only (Fig. S2)
python $S/01_build_compartments.py $STEP1_COMMON --lnc-types protein_coding --skip-overlap-filter \
    --out-dir $W/step1_pc                                                     # protein-coding positive control
python $S/20_attrition.py --gtf $GTF --lnc-genes $W/step1_all/lnc_genes.tsv.gz --mask-bed $D/sd.bed \
    --out $W/attrition                                                        # Table S7

stamp "step 2: TE copies"
python $S/02_select_te_copies.py $TE_SELECT --compartments $W/step1_all/compartments.tsv.gz --out-dir $W/step2_v2

stamp "step 9: ordinary-sequence windows (lncRNA and protein-coding)"
python $S/09_nonte_windows.py --compartments $W/step1_all/compartments.tsv.gz $NONTE_WINDOWS \
    --chroms $(printf "chr%s," {1..22} | sed 's/,$//') --out-dir $W/nonte                        # autosomes only
python $S/09_nonte_windows.py --compartments $W/step1_pc/compartments.tsv.gz $NONTE_WINDOWS \
    --chroms $(printf "chr%s," {1..22} | sed 's/,$//') --out-dir $W/pc_windows

stamp "steps 11 and 21: RNA properties (ENCODE export) with tau recomputed from GTEx v10"
python $S/11_gene_properties.py --export $RNA_PROPS --out $W/gene_props_export.tsv
python $S/21_tau_gtex.py --gct $GCT --gene-table $W/gene_props_export.tsv --out $W/gene_props.tsv
python - "$RNA_PROPS" "$W/gene_class.tsv" << 'PY'
import sys
import pandas as pd
e = pd.read_csv(sys.argv[1], sep="\t", usecols=["gene_id", "project_class"])
e.to_csv(sys.argv[2], sep="\t", index=False)
print(e.project_class.value_counts().head(4).to_string())
PY
stamp "done"
