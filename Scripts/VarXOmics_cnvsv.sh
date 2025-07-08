#!/bin/bash

# Script: CNV and SV annotation and vcf2circos plotting
# Usage: 
# sh VarXOmics_cnvsv.sh \
#     --input-sample P00110_8 \
#     -hpo HP:0003002,HP:0006625 \
#     -o /path/to/user/dir
#     -v /path/to/user/vcf

# Default values
# for local 

ANNOTSV='/mnt/storage_pool/Genomics/VarXOmics/tools/AnnotSV/bin/AnnotSV'
SCRIPTS_DIR='/mnt/storage_pool/Genomics/VarXOmics/Scripts'

# Function to display usage
usage() {
    cat << EOF
Usage: $0 [OPTIONS]

OPTIONS:
    -i, --input-sample SAMPLE       Input sample name (required)
    -o, --output-directory DIR      Output directory path (required)
    -v, --vcf VCF_FILE              Input VCF file (required)
    -g, --genome GENOME             Reference Genome, default is GRCH38. 
                                    Options: GRCH37, GRCH38
    --data-type 		    Input data type, options: CNVs or SVs (required)
    -hpo HPO Terms                  Run Exomiser priortization for HPO terms (optional)
    -h, --help                      Display this help message

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
GENOME='GRCH38'
HPO_TERMS=''
DATA_TYPE=''

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input-sample)
            INPUT_SAMPLE="$2"
            shift 2
            ;;
        -v|--vcf)
            VCF_FILE="$2"
            shift 2
            ;;
        -g|--genome)
            GENOME="$2"
            shift 2
            ;;
        -o|--output-directory)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -hpo)
            HPO_TERMS="$2"
            shift 2
            ;;
	--data-type)
	    DATA_TYPE="$2"
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

if [[ -z "$OUTPUT_DIR" ]]; then
    echo "Error: --output-directory is required"
    usage
    exit 1
fi

if [[ -z "$DATA_TYPE" ]]; then
    echo "Error: --data-type is required"
    usage
    exit 1
fi

if [[ -z "$VCF_FILE" ]]; then
    echo "Error: -v is required"
    usage
    exit 1
fi


# Display configuration
echo "=== Exomiser Analysis Configuration ==="
echo "Input Sample: $INPUT_SAMPLE"
echo "HPO Terms: $VCF_FILE"
echo "Output Directory: $OUTPUT_DIR"
echo "Reference is $GENOME"
echo "Data type is $DATA_TYPE"
echo "======================================"

# Step 1: SV annotation
# Genome in AnnotSV is GRCh38 or GRCh37, not capital H

if [[ $GENOME == "GRCH38" ]]; then 
    AnnotSV_GENOME="GRCh38"
elif [[ $GENOME == "GRCH37" ]]; then
    AnnotSV_GENOME="GRCh37"
fi

/mnt/storage_pool/Genomics/VarXOmics/tools/AnnotSV/bin/AnnotSV \
  -SVinputFile "$VCF_FILE" \
  -outputDir "$OUTPUT_DIR" \
  -overwrite 1 \
  -genomeBuild $AnnotSV_GENOME \
  -overwrite 1 -genomeBuild "$AnnotSV_GENOME" -snvIndelPASS 1 -vcf 1 -variantconvertMode full

rm -rf $OUTPUT_DIR/*.unannotated.tsv
rm -rf $OUTPUT_DIR/*AnnotSV_inputSV*.tsv
rm -rf $OUTPUT_DIR/*AnnotSV_inputSVfile.formatted.sorted*

# Step 2: vcf2circos
# make sure the extracted file does not exist
rm -rf $OUTPUT_DIR/${INPUT_SAMPLE}.short.vcf2circos.vcf
rm -rf $OUTPUT_DIR/${INPUT_SAMPLE}.short.vcf
rm -rf $OUTPUT_DIR/${INPUT_SAMPLE}.short.modified.vcf 

# extract SVtypes, svlen, and End

wholename=$(find $OUTPUT_DIR -type f -name "*.annotated.vcf")
filename=$(basename $wholename)

if [[ "$DATA_TYPE" == "cnv" ]]; then
    bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\tEND=%INFO/END;SVTYPE=%INFO/SVTYPE;SVLEN=%INFO/SVLEN\t%FORMAT[\t%SAMPLE]\n' \
        "$OUTPUT_DIR/$filename" > "$OUTPUT_DIR/${INPUT_SAMPLE}.short.vcf"
elif [[ "$DATA_TYPE" == "sv" ]]; then
    bcftools view -e 'INFO/SVLEN="."' "$OUTPUT_DIR/$filename" | \
    bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\tEND=%INFO/END;SVTYPE=%INFO/SVTYPE;SVLEN=%INFO/SVLEN\t%FORMAT[\t%SAMPLE]\n' \
        > "$OUTPUT_DIR/${INPUT_SAMPLE}.short.vcf"
fi

# some warnings will appear here, but no problem 

# change ; into |
sed 's/;/|/g' $OUTPUT_DIR/${INPUT_SAMPLE}.short.vcf > $OUTPUT_DIR/${INPUT_SAMPLE}.short.modified.vcf

# Add INFO atthe beginning
awk 'BEGIN{OFS="\t"} 
     /^#/ {print; next} 
     { $8 = "INFO=" $8; print }' $OUTPUT_DIR/${INPUT_SAMPLE}.short.modified.vcf > "$OUTPUT_DIR/${INPUT_SAMPLE}.short.vcf2circos.vcf"
	     
VCF2CIRCOS='/mnt/storage_pool/Genomics/miniconda3/envs/vep114/bin/vcf2circos'
#VCF2CIRCOS='/mnt/storage_pool/Genomics/miniconda3/envs/vep114/bin/vcf2circos'
OPTION='/mnt/storage_pool/Genomics/VarXOmics/tools/vcf2circos/config_vcf2circos_21082023/Static/options.json'

if [[ GENOME == "GRCH38" ]]; then
    ASSEMBLY="hg38"
elif [[ GENOME == "GRCH37" ]]; then
    ASSEMBLY="hg37"
fi

# html output
$VCF2CIRCOS -i "$OUTPUT_DIR/${INPUT_SAMPLE}.short.vcf2circos.vcf" \
    -o $OUTPUT_DIR/${INPUT_SAMPLE}_vcf2circos.html \
    -p $OPTION \
    -a "$ASSEMBLY"

echo "$(date): Analysis completed successfully!"
echo "Results are available in: $OUTPUT_DIR"

# Step 3: Exomiser (if required)
if [[ -n $HPO_TERMS && $HPO_TERMS != "" ]]; then
    # create config file
    EXOMISER_CONFIG="$OUTPUT_DIR/$INPUT_SAMPLE-Exomiser-genome-PHENIX_PRIORITY.yml"
    cp "$SCRIPTS_DIR/sample-Exomiser-genome-PHENIX_PRIORITY.yml" "$EXOMISER_CONFIG"

    # add HPO inside
    FORMATTED_HPO=$(echo "$HPO_TERMS" | sed "s/[^,]*/'&'/g")
    sed -i "s|hpoIds:|hpoIds: [$HPO_TERMS]|" "$EXOMISER_CONFIG"
    echo "HPO terms configured: $HPO_TERMS"

    # Add VCF file 
    echo "$(date): Configuring VCF file..."
    sed -i "s|vcf:|vcf: $VCF_FILE|" "$EXOMISER_CONFIG"

    # for server
    SCRIPTS_DIR='/mnt/storage_pool/Genomics/VarXOmics/Scripts'
    EXOMISER_DIR='/mnt/storage_pool/Genomics/VarXOmics/tools/exomiser-cli-14.0.0'

    if [[ $GENOME == "GRCH38" ]]; then
        ASSEMBLY="hg38"
    elif [[ $GENOME == "GRCH37" ]]; then 
        ASSEMBLY="hg37"
    fi

    echo "$(date): Running Exomiser analysis..."
    if ! java -jar "$EXOMISER_DIR/exomiser-cli-14.0.0.jar" \
        --analysis "$EXOMISER_CONFIG" \
        --assembly $ASSEMBLY  \
        --output-directory "$OUTPUT_DIR" \
        --output-filename "${INPUT_SAMPLE}_Exomiser_PHENIX_PRIORITY" \
        --exomiser.data-directory="$EXOMISER_DIR/Data" \
        --spring.config.location="$EXOMISER_DIR/application.properties"; then
        echo "Error: Exomiser analysis failed"
        exit 1
    fi

    echo "$(date): Analysis completed successfully!"
    echo "Results are available in: $OUTPUT_DIR"
fi

