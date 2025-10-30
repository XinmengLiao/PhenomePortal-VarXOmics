#!/bin/bash

# Script: CNV and SV annotation and vcf2circos plotting
# Usage: 
# sh VarXOmics_cnvsv.sh \
#     --input-sample P00110_8 \
#     -o /path/to/user/dir \
#     -v /path/to/user/vcf --data-type cnv/sv

# Default values
ANNOTSV='/mnt/storage_pool/Genomics/VarXOmics/tools/AnnotSV/bin/AnnotSV'
SCRIPTS_DIR='/mnt/storage_pool/Genomics/VarXOmics/Scripts'

# Function to display usage
usage() {
    cat << EOF
Usage: $0 [OPTIONS]

OPTIONS:
    -i, --input-file FILE           Input txt file for queried CNVs and SVs (required)
    --id							Queried variant information, i.e. chr13_32338103_G_GAATTTCCG (required). 
    								This would be used as the output files prefix, normally is the queried variant. 
    -o, --output-directory DIR      Output directory path (required)
    -g, --genome GENOME             Reference Genome, default is GRCH38. 
                                    Options: GRCH37, GRCH38
    -h, --help                      Display this help message

REQUIRED ARGUMENTS:
    - input-file: Sample ID
    - output-directory: Output directory path
    - vcf: User uploaded vcf file
EOF
}

# Initialize variables
INPUT_FILE=""
ID=""
OUTPUT_DIR=""
GENOME='GRCH38'

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input-file)
            INPUT_FILE="$2"
            shift 2
            ;;
        -g|--genome)
            GENOME="$2"
            shift 2
            ;;
        --id)
            ID="$2"
            shift 2
            ;;
        -o|--output-directory)
            OUTPUT_DIR="$2"
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
if [[ -z "$INPUT_FILE" ]]; then
    echo "Error: --input-file is required"
    usage
    exit 1
fi

if [[ -z "$OUTPUT_DIR" ]]; then
    echo "Error: --output-directory is required"
    usage
    exit 1
fi

if [[ -z "$ID" ]]; then
    echo "Error: --id is required"
    usage
    exit 1
fi


# Display configuration
echo "=== Exomiser Analysis Configuration ==="
echo "Sample ID: $ID"
echo "Input Sample: $INPUT_FILE"
echo "Output Directory: $OUTPUT_DIR"
echo "Reference is $GENOME"
echo "======================================"

mkdir -p $OUTPUT_DIR

# Step 1: Convert txt into vcf file 
# Genome in AnnotSV is GRCh38 or GRCh37, not capital H

if [[ $GENOME == "GRCH38" ]]; then 
    AnnotSV_GENOME="GRCh38"
elif [[ $GENOME == "GRCH37" ]]; then
    AnnotSV_GENOME="GRCh37"
fi

if grep -q '^chr' $INPUT_FILE; then 
    sed -i 's/^chr//' $INPUT_FILE
fi

echo -e "#CHROM\tStart\tEnd\tSV_id\tSV_type" > "$OUTPUT_DIR/$ID.bed"
awk -F"_" '{print $1"\t"$2"\t"$3"\t"$1":"$2"_"$3"\t"$4}' $INPUT_FILE >> $OUTPUT_DIR/$ID.bed

# Step 2: SV annotation
conda run -n vep114 /mnt/nas/Genomics/VarXOmics/tools/AnnotSV/bin/AnnotSV \
  -SVinputFile "$OUTPUT_DIR/$ID.bed" \
  -outputDir "$OUTPUT_DIR" \
  -overwrite 1 \
  -genomeBuild $AnnotSV_GENOME \
  -overwrite 1 -genomeBuild "$AnnotSV_GENOME" -snvIndelPASS 1 -vcf 1 -variantconvertMode full

rm -rf $OUTPUT_DIR/*.unannotated.tsv
rm -rf $OUTPUT_DIR/*AnnotSV_inputSV*.tsv
rm -rf $OUTPUT_DIR/*AnnotSV_inputSVfile.formatted.sorted*

# # Step 2: vcf2circos
# # make sure the extracted file does not exist
# rm -rf $OUTPUT_DIR/${ID}.short.vcf2circos.vcf
# rm -rf $OUTPUT_DIR/${ID}.short.vcf2circos.tmp.vcf
# rm -rf $OUTPUT_DIR/${ID}.short.vcf
# rm -rf $OUTPUT_DIR/${ID}.short.modified.vcf 

# # # extract SVtypes, svlen, and End

# wholename=$(find $OUTPUT_DIR -type f -name "*.annotated.vcf")
# filename=$(basename $wholename)

# if [[ "$DATA_TYPE" == "cnv" ]]; then
#     conda run -n vep114 bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\tEND=%INFO/END;SVTYPE=%INFO/SVTYPE;SVLEN=%INFO/SVLEN\t%FORMAT[\t%SAMPLE]\n' \
#         "$OUTPUT_DIR/$filename" >> "$OUTPUT_DIR/${ID}.short.vcf"
# elif [[ "$DATA_TYPE" == "sv" ]]; then
#     conda run -n vep114 bcftools view -e 'INFO/SVLEN="."' "$OUTPUT_DIR/$filename" -o "$OUTPUT_DIR/${ID}.filtered.vcf"
#     conda run -n vep114 bcftools query \
#         -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\tEND=%INFO/END;SVTYPE=%INFO/SVTYPE;SVLEN=%INFO/SVLEN\t%FORMAT[\t%SAMPLE]\n' \
#         "$OUTPUT_DIR/${ID}.filtered.vcf" >> "$OUTPUT_DIR/${ID}.short.vcf"
#     #conda run -n vep114 bcftools view -e 'INFO/SVLEN="."' "$OUTPUT_DIR/$filename" | \
#     #conda run -n vep114 bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\tEND=%INFO/END;SVTYPE=%INFO/SVTYPE;SVLEN=%INFO/SVLEN\t%FORMAT[\t%SAMPLE]\n' \
#     #    >> "$OUTPUT_DIR/${ID}.short.vcf"
# fi

# # # some warnings will appear here, but no problem 

# # # change ; into |
# sed 's/;/|/g' $OUTPUT_DIR/${ID}.short.vcf > $OUTPUT_DIR/${ID}.short.modified.vcf

# # # Add INFO at the beginning
# awk 'BEGIN{OFS="\t"} 
#      /^#/ {print; next} 
#      { $8 = "INFO=" $8; print }' $OUTPUT_DIR/${ID}.short.modified.vcf > "$OUTPUT_DIR/${ID}.short.vcf2circos.tmp.vcf"

# sed '$d' "$OUTPUT_DIR/${ID}.short.vcf2circos.tmp.vcf" > "$OUTPUT_DIR/${ID}.short.vcf2circos.vcf"

	     
# VCF2CIRCOS='/mnt/storage_pool/Genomics/miniconda3/envs/vep114/bin/vcf2circos'
# OPTION='/mnt/storage_pool/Genomics/VarXOmics/tools/vcf2circos/config_vcf2circos_21082023/Static/options.json'

# if [[ GENOME == "GRCH38" ]]; then
#     ASSEMBLY="hg38"
# elif [[ GENOME == "GRCH37" ]]; then
#     ASSEMBLY="hg37"
# fi

# # # html output
# conda run -n vep114 $VCF2CIRCOS -i "$OUTPUT_DIR/${ID}.short.vcf2circos.vcf" \
#     -o $OUTPUT_DIR/${ID}_vcf2circos.html \
#     -p $OPTION \
#     -a "$ASSEMBLY"
