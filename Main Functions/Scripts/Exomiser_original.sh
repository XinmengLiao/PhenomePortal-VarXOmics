#!/bin/bash

# Script: Exomiser for REMM score (non-coding variant prediction) and HIPHIVE_PRIORITY
# Usage: sh exomiser_analysis.sh --input-sample P00110_8 --hpo-terms HP:0003002,HP:0006625 --inheritance AD --user-directory /path/to/user/dir [--vcf file.vcf]

# Default values
# for local
SCRIPTS_DIR='/Users/xinmengliao/Documents/Project/20250516_Webserver/Scripts/'
EXOMISER_DIR='/Users/xinmengliao/Documents/Project/20250516_Webserver/exomiser-cli-14.0.0'
# for server
SCRIPTS_DIR='/mnt/storage_pool/Genomics/VarXOmics/Scripts'
EXOMISER_DIR='/mnt/storage_pool/Genomics/VarXOmics/tools/exomiser-cli-14.0.0'

# Function to display usage
usage() {
    cat << EOF
Usage: $0 [OPTIONS]

OPTIONS:
    -i, --input-sample SAMPLE       Input sample name (required)
    -o, --output-directory DIR      Output directory path (required)
    -v, --vcf FILE                  Path for input VCF file (required)
    -hpo, --hpo-terms TERMS         HPO terms separated by commas (optional)
                                    If HPO terms are not provided, only the multi-omics aggregation scores will be calculated.
                                    e.g. HP:0003002,HP:0006625
    -s, --scripts-dir DIR           Scripts directory path (default: $SCRIPTS_DIR)
    -e, --exomiser-dir DIR          Exomiser directory path (default: $EXOMISER_DIR)
    -h, --help                      Display this help message

EXAMPLES:
    With HPO terms: 
    $0 --input-sample P001 --hpo-terms HP:0003002,HP:0006625 --inheritance AD --output-directory /path/to/results

    Without HPO terms:
    $0 -i P001 -o /path/to/results -v sample.vcf

REQUIRED ARGUMENTS:
    - input-sample: Sample ID
    - output-directory: Output directory path
    - vcf: User uploaded vcf file
EOF
}

# Initialize variables
INPUT_SAMPLE=""
HPO_TERMS=""
OUTPUT_DIR=""
VCF_FILE=""
GENDER=""

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input-sample)
            INPUT_SAMPLE="$2"
            shift 2
            ;;
        -hpo|--hpo-terms)
            HPO_TERMS="$2"
            shift 2
            ;;
        -o|--output-directory)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -v|--vcf)
            VCF_FILE="$2"
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
    echo "Error: --input-sample is required"
    usage
    exit 1
fi

if [[ -z "$OUTPUT_DIR" ]]; then
    echo "Error: --output-directory is required"
    usage
    exit 1
fi

# Validate directories exist
if [[ ! -d "$SCRIPTS_DIR" ]]; then
    echo "Error: Scripts directory does not exist: $SCRIPTS_DIR"
    exit 1
fi

if [[ ! -d "$EXOMISER_DIR" ]]; then
    echo "Error: Exomiser directory does not exist: $EXOMISER_DIR"
    exit 1
fi

if [[ ! -d "$OUTPUT_DIR" ]]; then
    echo "Error: Output directory does not exist: $OUTPUT_DIR"
    exit 1
fi

# Set VCF file if not provided
if [[ -z "$VCF_FILE" ]]; then
    VCF_FILE="$OUTPUT_DIR/$INPUT_SAMPLE.vcf.gz"
fi

# Display configuration
echo "=== Exomiser Analysis Configuration ==="
echo "Input Sample: $INPUT_SAMPLE"
echo "HPO Terms: $HPO_TERMS"
echo "Output Directory: $OUTPUT_DIR"
echo "VCF File: $VCF_FILE"
echo "Scripts Directory: $SCRIPTS_DIR"
echo "Exomiser Directory: $EXOMISER_DIR"
echo "======================================"

# Step 3.1: Prepare Exomiser configuration
echo "$(date): Preparing Exomiser configuration..."
EXOMISER_CONFIG="$OUTPUT_DIR/$INPUT_SAMPLE-Exomiser-genome-PHENIX_PRIORITY.yml"

if [[ ! -f "$SCRIPTS_DIR/sample-Exomiser-genome-PHENIX_PRIORITY.yml" ]]; then
    echo "Error: Template Exomiser config file not found: $SCRIPTS_DIR/sample-Exomiser-genome-PHENIX_PRIORITY.yml"
    exit 1
fi

cp "$SCRIPTS_DIR/sample-Exomiser-genome-PHENIX_PRIORITY.yml" "$EXOMISER_CONFIG"

echo "[DEBUG] Editing config file: $EXOMISER_CONFIG"

# Step 3.3: Configure HPO terms
echo "$(date): Configuring HPO terms..."
if [[ -n "$HPO_TERMS" && "$HPO_TERMS" != "" ]]; then
    FORMATTED_HPO=$(echo "$HPO_TERMS" | sed "s/[^,]*/'&'/g")
    sed -i "s|hpoIds:.*|hpoIds: [$HPO_TERMS]|" "$EXOMISER_CONFIG"
    echo "HPO terms configured: $HPO_TERMS"
else
    echo "No HPO terms provided, priortising scores will be the Multi-omics aggregation scores"
fi

# Step 3.4: Configure VCF file
echo "$(date): Configuring VCF file..."
sed -i "s|vcf:.*|vcf: $VCF_FILE|" "$EXOMISER_CONFIG"

# Step 3.5: Run Exomiser
echo "$(date): Running Exomiser analysis..."
if ! java -jar "$EXOMISER_DIR/exomiser-cli-14.0.0.jar" \
    --analysis "$EXOMISER_CONFIG" \
    --assembly hg38  \
    --output-directory "$OUTPUT_DIR" \
    --output-filename "${INPUT_SAMPLE}_Exomiser_PHENIX_PRIORITY" \
    --exomiser.data-directory="$EXOMISER_DIR/Data" \
    --spring.config.location="$EXOMISER_DIR/application.properties"; then
    echo "Error: Exomiser analysis failed"
    exit 1
fi

echo "$(date): Analysis completed successfully!"
echo "Results are available in: $OUTPUT_DIR"
