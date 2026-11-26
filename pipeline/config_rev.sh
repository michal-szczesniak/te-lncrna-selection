# Settings for the revision stages (s00r, s03r, s11, s12). Sourced by every revision stage (it sources config.sh).
# The *_STR variables can be overridden from the environment (used by the synthetic test).
source "$PIPE/config.sh"
R=$D/rev                           # downloads of the revision (chains, genomes, GTEx junctions)
RW=$W/rev                          # outputs of the revision analyses
mkdir -p $R $RW

# TE groups of the main analysis and extra groups (comma lists; read without filename expansion)
IFS=, read -r -a G4 <<< "${REV_GROUPS_STR:-AluS*,MIR*,AluJ*,L2abc}"
IFS=, read -r -a GX <<< "${REV_EXTRA_STR:-AluY*,MIR3}"
IFS=, read -r -a GALU <<< "${REV_ALU_STR:-AluS*,AluJ*}"
OLD_LIST=${REV_OLD_LIST:-MIR*,L2abc,L2d*,MIR3}           # groups old enough for mammalian phyloP
ALU_LIST=${REV_ALU_LIST:-AluS*,AluJ*,AluY*}              # primate-specific groups (primate phyloP)
N_CHR_GNOMAD=${N_CHR_GNOMAD:-23}

# other species for the age of splice sites, closest first (UCSC liftOver chains from hg38, 2bit genomes)
IFS=, read -r -a SPECIES <<< "${SPECIES_STR:-panTro6,gorGor6,ponAbe3,rheMac10,calJac4,mm39}"
species_args () {
    for sp in "${SPECIES[@]}"; do
        local cap="$(tr '[:lower:]' '[:upper:]' <<< ${sp:0:1})${sp:1}"
        local g=$R/$sp.2bit
        [ -s "$g" ] || g=$R/$sp.fa
        echo -n "--species $sp=$R/hg38To$cap.over.chain.gz=$g "
    done
}
GTEX_JUNC=$R/GTEx_Analysis_v10_STARv2.7.10a_junctions.gct.gz
GTEX_ATTR=$R/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt

# consensus segments for the spike-in power runs (bin-aligned, 20-nt bins), per group
declare -A SPIKE_SEGS=(["AluS*"]="101-120,101-140,101-160,all" ["AluJ*"]="101-120,101-140,101-160,all"
                       ["MIR*"]="81-100,81-120,81-140,all" ["L2abc"]="201-220,201-240,201-260,all")
SPIKE_DEFAULT=${SPIKE_DEFAULT:-101-120,101-140,all}
SPIKE_REPS=${SPIKE_REPS:-20}
N_PERM=${N_PERM:-1000}

# LTR and L1 extension: copies, consensus grouping by family patterns (31_make_groups.py; ';'-separated lists)
TE_SELECT_LTR=${TE_SELECT_LTR:-"--rmsk $D/hg38.fa.out.gz --rmsk-format rmout --subfamily-regex ^(THE1|MSTA|MSTB|MSTC|MSTD|MLT1|LTR12|LTR7|MLT2) --assign perbase --max-intronic 500 --max-intergenic 200"}
TE_SELECT_L1=${TE_SELECT_L1:-"--rmsk $D/hg38.fa.out.gz --rmsk-format rmout --subfamily-regex ^L1M --assign perbase --max-intronic 300 --max-intergenic 100"}
FAM_LTR_STR=${FAM_LTR_STR:-"THE1_g:^THE1[A-D]$;MST_g:^MST[A-D]1?$;MLT1A_g:^MLT1A[0-9]?$;MLT1BF_g:^MLT1[B-F][0-9]?$;MLT1GL_g:^MLT1[G-L][0-9]?$;MLT1MO_g:^MLT1[M-O][0-9]?$;LTR12_g:^LTR12;LTR7_g:^LTR7;MLT2_g:^MLT2"}
FAM_L1_STR=${FAM_L1_STR:-"L1MA_g:^L1MA[0-9]+$;L1MB_g:^L1MB[0-9]+$;L1MC_g:^L1MC[0-9]*[a-z]?$;L1MD_g:^L1MD[0-9]*$;L1ME_g:^L1ME[0-9]+[a-z]?$;L1M_g:^L1M[0-9]+$"}
MIN_EXONIC_GROUP=${MIN_EXONIC_GROUP:-200}
fam_args () {  # ';'-separated NAME:REGEX list -> one '--family' and one 'NAME:REGEX' token per line
    local fams f
    IFS=';' read -r -a fams <<< "$1"
    for f in "${fams[@]}"; do printf -- '--family\n%s\n' "$f"; done
}
