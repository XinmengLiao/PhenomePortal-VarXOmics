#!/bin/bash

# Message:
# Single variant query requires firstly transfer the input variant into vcf.gz format

# sh Scripts/SingleVariant.sh --id chr13_32338103_G_GA \
# 	-i /Users/xinmengliao/Documents/Project/20250516_Webserver/examples/chr13_32338103_G_GA/chr13_32338103_G_GA.vcf.gz \
# 	-o /Users/xinmengliao/Documents/Project/20250516_Webserver/examples/chr13_32338103_G_GA \
# 	--genome GRCH38 --only-pass yes --gender Female

usage() {
    cat << EOF
Usage: $0 [OPTIONS]

OPTIONS:
     --input-sample SAMPLE       Input sample compressed vcf file (required)
    -o, --output-directory DIR      Output directory path (required)
    --id							Queried variant information, i.e. chr13_32338103_G_GA (required). 
    								This would be used as the output files prefix 
    -g, --gender GENDER             Gender (optional)
                                    Options: Male, Female
    -h, --help                      Display this help message
    --only-pass                    	Only keep the PASS variant in vcf file. Default is yes. 
    								Options: yes, no. 
    --genome						Reference genome, default is GRCH38
    								Options: GRCH37, GRCH38


REQUIRED ARGUMENTS:
    - input-sample: Sample ID
    - output-directory: Output directory path
EOF
}

# scripts for steamline the whole pipeline
ID=""
INPUT_SAMPLE=""
OUTPUT_DIR=""
GENOME="GRCH38"
ONLY_PASS="yes"
GENDER=""

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
        --id)
            ID="$2"
            shift 2
            ;;
        --only-pass)
            ONLY_PASS="$2"
            shift 2
            ;;
        --genome)
            GENOME="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# check if the input vcf and output directory are provided
if [[ -z "$ID" || -z "$INPUT_SAMPLE" || -z "$OUTPUT_DIR" ]]; then
    echo "ERROR: --id, --input-sample (-i), and --output-directory (-o) are required."
    usage
    exit 1
fi


# Display configuration
echo "=== Single Variant Analysis Configuration ==="
echo "Input Variant Name: $ID"
echo "Output Directory: $OUTPUT_DIR"
echo "VCF File: $INPUT_SAMPLE"
echo "Reference genome version is: $GENOME"
echo "Filtered PASS: $ONLY_PASS"
echo "Gender: $GENDER"
echo "============================================="

SCRIPTS='/mnt/storage_pool/Genomics/VarXOmics/Scripts'
mkdir -p $OUTPUT_DIR

# Step 1: VEP
conda run -n vep114 bash $SCRIPTS/vep.sh \
	-v $INPUT_SAMPLE \
 	-o $OUTPUT_DIR \
 	-i $ID \
 	-g $GENOME \
	--only_pass $ONLY_PASS

# Step 2: Python (Except GeneBe can not be run due to the too permerssive)
mv ~/.netrc ~/.netrc.backup  # Remove /.netrc to remove the restriction from GeneBE
conda run -n vep114 python $SCRIPTS/vep_result_management_SingleVariant.py $OUTPUT_DIR/${ID}_vep_annotated.vcf.gz $OUTPUT_DIR/$ID.txt $GENDER 
mv ~/.netrc.backup ~/.netrc

# Step 3: R priortization
conda run -n varxomics Rscript $SCRIPTS/SingleVariant.R $OUTPUT_DIR/$ID.txt $OUTPUT_DIR

