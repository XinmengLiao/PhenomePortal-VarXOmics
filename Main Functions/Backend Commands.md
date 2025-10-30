Single Gene Query: 
```
varxomics='/mnt/storage_pool/Genomics/VarXOmics'
sh $varxomics/SingleGene.sh [Gene] [output-dir]

#example
sh $varxomics/SingleGene.sh APOB example/APOB
```

Single Variant Query:
```
varxomics='/mnt/storage_pool/Genomics/VarXOmics'
bash $varxomics/SingleVariant.sh \
	--id [FILENAME] \
	-i [INPUT_TEXT] \
	-o [OUTPUT_DIR] \
	--genome [GENOME] --only-pass [yes|no] --gender [Male|Female]

# example
bash $varxomics/SingleVariant.sh --id chr13_32338103_G_GA \
	-i /mnt/storage_pool/Genomics/VarXOmics/examples/chr13_32338103_G_GA/chr13_32338103_G_GA.txt \
	-o /mnt/storage_pool/Genomics/VarXOmics/examples/chr13_32338103_G_GA \
	--genome GRCH38 --gender Female
```
* Users input the text in **"CHR_POS_REF_ALT"** format.
* The text will be converted into VCF format in the first step. 

Single small variant analysis 
```
varxomics='/mnt/nas/Genomics/VarXOmics'
bash VarXOmics_snv.sh \
  -i test \
  -o $varxomics/examples/test \
  -v $varxomics/examples/test/test.vcf.gz \
  --genome GRCH38 --only-pass yes \
  --gender Female --inheritance AD \
  --hpo HP:0003002,HP:0006625
```

Multiple CNVs analysis 
```
bash VarXOmics_cnvsv.sh --input-sample P001_167 -hpo HP:0003002 -o examples/cnv/ -v examples/cnv_org/P001_167.cnv.vcf.gz --data-type cnv
```

Multiple SVs analysis 
```
varxomics="/mnt/nas/Genomics/VarXOmics"
bash $varxomics/SingleCNVSV.sh --id [id] -i [input_txt_file] -o [output_dir] --genome GRCH38

# example
bash $varxomics/SingleCNVSV.sh --id chr11_65736013_65736014_INS -i examples/chr11_65736013_65736014_INS/chr11_65736013_65736014_INS.txt -o examples/chr11_65736013_65736014_INS --genome GRCH38
```
