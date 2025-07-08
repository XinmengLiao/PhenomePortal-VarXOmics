#!/bin/bash

# Script: exomiser_analysis.sh
# Description: Run Exomiser analysis with HPO terms for genetic variant analysis
# Usage: sh exomiser_analysis.sh --input-sample P00110_8 --hpo-terms HP:0003002,HP:0006625 --inheritance AD --user-directory /path/to/user/dir [--vcf file.vcf]

# Default values
SCRIPTS_DIR='/mnt/storage_pool/Genomics/VarXOmics/Scripts'

# Function to display usage
usage() {
    cat << EOF
Usage: $0 [OPTIONS]

OPTIONS:
    -i, --input-sample SAMPLE       Input sample name (required)
    -o, --output-directory DIR      Output directory path (required)
    -v, --vcf FILE                  Path for input VCF file (required)
    -g, --gender GENDER             Gender (optional)
                                    Options: Male, Female
    -h, --help                      Display this help message
    --af-clinvar                    User defined allele frequency threshold for ClinVar variants, default is 1
    --af-precition                  User defined allele frequency threshold for predicted variants, default is 0.05
    --ada                           User defined ada threshold, default is 0.6
    --rf                            User defined rf threshold, default is 0.6
    --revel                         User defined revel threshold, default is 0.75
    --spliceai-al                   User defined SpliceAI acceptor loss threshold, default is 0.5
    --spliceai-ag                   User defined SpliceAI acceptor gain threshold, default is 0.5
    --spliceai-dl                   User defined SpliceAI donor loss threshold, default is 0.5
    --spliceai-dg                   User defined SpliceAI donor gain threshold, default is 0.5
    --bayesdel-addaf                User defined bayesdel (with AF) threshold, default is 0.0692655
    --bayesdel-noaf                 User defined bayesdel (no AF) threshold, default is 0.0570105
    --am-classification             User defined revel threshold, default is likely_pathogenic, ambigous
    --am-pathogenicity              User defined revel threshold, default is 0.564
    --clinvar                       User defined revel threshold, 
                                    default is Pathogenic, Likely_pathogenic, Uncertain_significance, Conflicting_classifications_of_pathogenicity
    --acmg-classification           User defined revel threshold, 
                                    default is Pathogenic, Likely_pathogenic, Uncertain_significance, Benign, Likely_benign


REQUIRED ARGUMENTS:
    - input-sample: Sample ID
    - output-directory: Output directory path
    - vcf: User uploaded vcf file
EOF
}

# Initialize variables
INPUT_SAMPLE=""
OUTPUT_DIR=""
VCF_FILE=""
GENDER=""

# Default thresholds
# AF_CLINVAR="1"
# AF_PRECITION="0.05"
# ADA="0.6"
# RF="0.6"
# REVEL="0.75"
# SPLICEAI_AL="0.5"
# SPLICEAI_AG="0.5"
# SPLICEAI_DL="0.5"
# SPLICEAI_DG="0.5"
# BAYESDEL_ADDAF="0.0692655"
# BAYESDEL_NOAF="-0.0570105"
# AM_CLASSIFICATION="likely_pathogenic,ambigous"
# AM_PATHOGENICITY="0.564"
# CLINVAR="Pathogenic,Likely_pathogenic,Uncertain_significance,Conflicting_classifications_of_pathogenicity"
# ACMG_CLASSIFICATION="Pathogenic,Likely_pathogenic,Uncertain_significance,Benign,Likely_benign"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input-sample)
            INPUT_SAMPLE="$2"
            shift 2
            ;;
        -o|--output-directory)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -g|--gender)
            GENDER="$2"
            shift 2
            ;;
        -v|--vcf)
            VCF_FILE="$2"
            shift 2
            ;;
        -s|--scripts-dir)
            SCRIPTS_DIR="$2"
            shift 2
            ;;
        -e|--exomiser-dir)
            EXOMISER_DIR="$2"
            shift 2
            ;;
        --af-clinvar)
            AF_CLINVAR="$2"
            shift 2
            ;;
        --af-precition)
            AF_PRECITION="$2"
            shift 2
            ;;
        --ada)
            ADA="$2"
            shift 2
            ;;
        --rf)
            RF="$2"
            shift 2
            ;;
        --revel)
            REVEL="$2"
            shift 2
            ;;
        --spliceai-al)
            SPLICEAI_AL="$2"
            shift 2
            ;;
        --spliceai-ag)
            SPLICEAI_AG="$2"
            shift 2
            ;;
        --spliceai-dl)
            SPLICEAI_DL="$2"
            shift 2
            ;;
        --spliceai-dg)
            SPLICEAI_DG="$2"
            shift 2
            ;;
        --bayesdel-addaf)
            BAYESDEL_ADDAF="$2"
            shift 2
            ;;
        --bayesdel-noaf)
            BAYESDEL_NOAF="$2"
            shift 2
            ;;
        --am-classification)
            AM_CLASSIFICATION="$2"
            shift 2
            ;;
        --am-pathogenicity)
            AM_PATHOGENICITY="$2"
            shift 2
            ;;
        --clinvar)
            CLINVAR="$2"
            shift 2
            ;;
        --acmg-classification)
            ACMG_CLASSIFICATION="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Error: Unknown option $1"
            usage
            exit 1
            ;;
    esac
done

# Validate required arguments
if [[ -z "$INPUT_SAMPLE" ]]; then
    echo "Error: --input-sample is required"
    usage
    exit 1
fi

if [[ -z "$VCF_FILE" ]]; then
    echo "Error: --vcf is required"
    usage
    exit 1
fi

if [[ -z "$OUTPUT_DIR" ]]; then
    echo "Error: --output-directory is required"
    usage
    exit 1
fi

# Validate directories exist
if [[ ! -d "$OUTPUT_DIR" ]]; then
    echo "Error: Output directory does not exist: $OUTPUT_DIR"
    exit 1
fi

# Validate VCF file exists
if [[ ! -f "$VCF_FILE" ]]; then
    echo "Error: VCF file does not exist: $VCF_FILE"
    exit 1
fi


# Display configuration
echo "=== Sample Infomation ==="
echo "Input Sample: $INPUT_SAMPLE"
echo "Gender: $GENDER"
echo "Output Directory: $OUTPUT_DIR"
echo "VCF File: $VCF_FILE"
echo "Scripts Directory: $SCRIPTS_DIR"
echo "Python Path: $PYTHON_PATH"
echo "=== Threshold Parameters ==="
echo "AF ClinVar: $AF_CLINVAR"
echo "AF Precition: $AF_PRECITION"
echo "ADA: $ADA"
echo "RF: $RF"
echo "REVEL: $REVEL"
echo "SpliceAI AL: $SPLICEAI_AL"
echo "SpliceAI AG: $SPLICEAI_AG"
echo "SpliceAI DL: $SPLICEAI_DL"
echo "SpliceAI DG: $SPLICEAI_DG"
echo "BayesDel AddAF: $BAYESDEL_ADDAF"
echo "BayesDel NoAF: $BAYESDEL_NOAF"
echo "AM Classification: $AM_CLASSIFICATION"
echo "AM Pathogenicity: $AM_PATHOGENICITY"
echo "ClinVar: $CLINVAR"
echo "ACMG Classification: $ACMG_CLASSIFICATION"
echo "======================================"


# Step 1: VEP Annotation

# Step 2: Python VEP result management
echo "$(date): Running Python for VEP result management..."
if ! python  "$SCRIPTS_DIR/vep_result_management_MultiVariant0701.py" \
    "$VCF_FILE" \
    "$OUTPUT_DIR/${INPUT_SAMPLE}.nodup.txt" \
    "$GENDER" \
    "$AF_CLINVAR" \
    "$AF_PRECITION" \
    "$ADA" "$RF" \
    "$REVEL" \
    "$SPLICEAI_AL" \
    "$SPLICEAI_DG" \
    "$SPLICEAI_DL" \
    "$SPLICEAI_AG" \
    "$BAYESDEL_ADDAF" \
    "$BAYESDEL_NOAF" \
    "$AM_CLASSIFICATION" \
    "$AM_PATHOGENICITY" \
    "$CLINVAR" \
    "$ACMG_CLASSIFICATION"; then
    echo "Error: Python VEP result management failed"
    exit 1
fi

echo "$(date): Python analysis completed successfully!"
echo "Results are available in: $OUTPUT_DIR"
