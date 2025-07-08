#!/bin/bash


#sh Scripts/SingleVariant.sh chr13_32338103_G_GA \
#	/Users/xinmengliao/Documents/Project/20250516_Webserver/examples/chr13_32338103_G_GA/chr13_32338103_G_GA.vcf.gz \
#	/Users/xinmengliao/Documents/Project/20250516_Webserver/examples/chr13_32338103_G_GA \
#	GRCH38 yes Female

# scripts for steamline the whole pipeline
ID=$1
INPUT_VCF=$2
OUTPUT_DIR=$3
GENOME=$4
ONLY_PASS=$5
GENDER=$6

SCRIPTS='/mnt/storage_pool/Genomics/VarXOmics/Scripts'

#source /mnt/storage_pool/Genomics/miniconda3/etc/profile.d/conda.sh
#conda activate vep114

# Step 1: VEP
conda run -n vep114 bash $SCRIPTS/vep.sh \
	-v $2 \
 	-o $3 \
 	-i $1 \
 	-g $4 \
	--only_pass $5

# Step 2: Python (Except GeneBe can not be run due to the too permerssive)
conda run -n vep114 python $SCRIPTS/vep_result_management_SingleVariant.py $OUTPUT_DIR/${ID}_vep_annotated.vcf.gz $3/$ID.txt $6 

# Step 3: R priortization
conda run -n varxomics Rscript $SCRIPTS/SingleVariant.R $3/$ID.txt $3

