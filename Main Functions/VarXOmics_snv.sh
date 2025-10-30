#!/bin/bash

# scripts for steamline the whole pipeline
SAMPLEID=$1
SCRIPTS='/mnt/storage_pool/Genomics/VarXOmics/Scripts'
CONFIG_FILE="/mnt/storage_pool/Genomics/VarXOmics/examples/${1}/${1}_multivariant_config.yaml"

source ~/miniconda3/etc/profile.d/conda.sh
conda activate vep114

echo "User uploaded sample file is $(yq '.input_vcf' $CONFIG_FILE)"
echo "Output directory is $(yq '.output_dir' $CONFIG_FILE)"
echo "User defined reference genome is $(yq '.genome' $CONFIG_FILE)"
echo "Sample ID is $(yq '.sampleID' $CONFIG_FILE)"

# Step 1: VEP
# sh $SCRIPTS/vep.sh \
# 	-v $(yq '.input_vcf' $CONFIG_FILE) \
# 	-o $(yq '.output_dir' $CONFIG_FILE) \
# 	-i $(yq '.sampleID' $CONFIG_FILE) \
# 	-g $(yq '.genome' $CONFIG_FILE) \
# 	--only-pass $(yq '.only_pass' $CONFIG_FILE)

# Step 2: Python (Except GeneBe can not be run due to the too permerssive)
bash $SCRIPTS/python_MultiVariant.sh \
	-v $(yq -r '.input_vcf' $CONFIG_FILE) \
	-o $(yq -r '.output_dir' $CONFIG_FILE) \
	-g $(yq -r '.gender' $CONFIG_FILE) \
	-i $(yq -r '.sampleID' $CONFIG_FILE) \
	--af-clinvar $(yq -r '.AF_CLINVAR' $CONFIG_FILE) \
	--af-precition $(yq -r '.AF_PRECITION' $CONFIG_FILE) \
	--ada $(yq -r '.ADA' $CONFIG_FILE) \
	--rf $(yq -r '.RF' $CONFIG_FILE) \
	--revel $(yq -r '.REVEL' $CONFIG_FILE) \
	--spliceai-al $(yq -r '.SPLICEAI_AL' $CONFIG_FILE) \
	--spliceai-ag $(yq -r '.SPLICEAI_AG' $CONFIG_FILE) \
	--spliceai-dl $(yq -r '.SPLICEAI_DL' $CONFIG_FILE) \
	--spliceai-dg $(yq -r '.SPLICEAI_DG' $CONFIG_FILE) \
	--bayesdel-addaf $(yq -r '.BAYESDEL_ADDAF' $CONFIG_FILE) \
	--bayesdel-noaf $(yq -r '.BAYESDEL_NOAF' $CONFIG_FILE) \
	--am-classification $(yq -r '.AM_CLASSIFICATION' $CONFIG_FILE) \
	--am-pathogenicity $(yq -r '.AM_PATHOGENICITY' $CONFIG_FILE) \
	--clinvar $(yq -r '.CLINVAR' $CONFIG_FILE) \
	--acmg-classification $(yq -r '.ACMG_CLASSIFICATION' $CONFIG_FILE)

# Step 3: Exomiser orignial 
HPO_VALUE=$(yq '.hpo' "$CONFIG_FILE")
if [ "$HPO_VALUE" != "null" ] && [ -n "$HPO_VALUE" ]; then
	bash $SCRIPTS/Exomiser_REMM.sh \
		-v $(yq -r '.input_vcf' $CONFIG_FILE) \
		-o $(yq -r '.output_dir' $CONFIG_FILE) \
		-i $(yq -r '.sampleID' $CONFIG_FILE) \
		-hpo $(yq -r '.hpo' $CONFIG_FILE)

	bash $SCRIPTS/Exomiser_original.sh \
		-v $(yq -r '.input_vcf' $CONFIG_FILE) \
		-o $(yq -r '.output_dir' $CONFIG_FILE) \
		-i $(yq -r '.sampleID' $CONFIG_FILE) \
		-hpo $(yq -r '.hpo' $CONFIG_FILE)
fi 

# Step 4: R priortization
Rscript $SCRIPTS/variant_priortizeV3.R $(yq -r '.sampleID' $CONFIG_FILE) $(yq -r '.hpo' $CONFIG_FILE) $(yq -r '.inheritance' $CONFIG_FILE) $(yq -r '.output_dir' $CONFIG_FILE) 

# Step 5: GO-KEGG results and plots
mkdir -p $(yq -r '.output_dir' $CONFIG_FILE)/Figures
Rscript $SCRIPTS/GO-KEGG-all.R $(yq -r '.sampleID' $CONFIG_FILE) $(yq -r '.prioritised_file' $CONFIG_FILE) 

# Step 6: Other plots
mkdir -p $(yq -r '.output_dir' $CONFIG_FILE)/Figures_data
Rscript $SCRIPTS/RPlotting-Multivariant.R $(yq -r '.sampleID' $CONFIG_FILE) $(yq -r '.prioritised_file' $CONFIG_FILE) $(yq -r '.output_dir' $CONFIG_FILE)


# If user selected specific variants, regenerate the GO-KEGG plots and all other plots 
#sh $SCRIPTS/GO-KEGG-all.R $(yq '.sampleID' $CONFIG_FILE) $(yq '.output_dir' $CONFIG_FILE) $(yq '.user_selected_go_kegg_geneset' $CONFIG_FILE)  
#sh $SCRIPTS/RPlotting-Multivariant.R $(yq '.sampleID' $CONFIG_FILE) $(yq '.user_selected_variant' $CONFIG_FILE) $(yq '.output_dir' $CONFIG_FILE)


