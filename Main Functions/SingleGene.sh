#!/bin/bash

# scripts for steamline the whole pipeline
GENE=$1
SCRIPTS='/mnt/storage_pool/Genomics/VarXOmics/Scripts'
OUTPUT_DIR=$2

# Step 1: R priortization

mkdir -p $OUTPUT_DIR
conda run -n varxomics Rscript $SCRIPTS/SingleGene.R $GENE $OUTPUT_DIR

