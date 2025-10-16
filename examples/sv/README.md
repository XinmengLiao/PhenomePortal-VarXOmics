These are the files would be used in constructing the webpage for CNV and SV. (Actually CNVs and SVs commands and results are same)

## Single CNV/SV query
Command:
```
varxomics="/mnt/nas/Genomics/VarXOmics"
bash $varxomics/SingleCNVSV.sh --id [id] -i [input_txt_file] -o [output_dir] --genome GRCH38

# example
bash $varxomics/SingleCNVSV.sh --id chr11_65736013_65736014_INS -i examples/chr11_65736013_65736014_INS/chr11_65736013_65736014_INS.txt -o examples/chr11_65736013_65736014_INS --genome GRCH38
```

#### Variant Summary
 - Table:
    - CNV: `chr10_17268558_17274050_DEL.annotated.tsv`
    - SV: `chr11_65736013_65736014_INS.annotated.tsv`
 - Decipher


## Multiple CNVs analysis
Command:
```bash
# for CNV
bash VarXOmics_cnvsv.sh --input-sample P001_167 -hpo HP:0003002 -o examples/cnv -v examples/cnv_org/P001_167.cnv.vcf.gz --data-type cnv
```
Input: `P001_167.cnv.vcf.gz`

### UI design

#### Variant Summary 
  - Table: `P001_167.cnv.annotated.tsv`
  - Circos figure: `P001_167_vcf2circos.html`
  - Decipher

#### Exomiser Prioritization
`P001_167_Exomiser_PHENIX_PRIORITY.json`

## Multiple SVs analysis
Command:
```bash
bash VarXOmics_cnvsv.sh -v examples/cnv/P001_167.cnv.vcf.gz -o examples/cnv -i P001_167 -g GRCH38 --data-type sv
```

Input: `P001_167.sv.vcf.gz`

### UI design
#### Variant Summary
 - Table: `P001_167.sv.annotated.tsv`
 - Circos figure: `P001_167_vcf2circos.html`
 - Decipher 

#### Exomiser Prioritization
 - `P001_167_Exomiser_PHENIX_PRIORITY.json` or other formats
