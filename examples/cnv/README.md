These are the files would be used in constructing the webpage for CNV. (Actually CNVs and SVs commands and results are same)

Command:
```bash
# for CNV
bash VarXOmics_cnvsv.sh --input-sample P001_167 -hpo HP:0003002 -o examples/cnv/ -v examples/cnv_org/P001_167.cnv.vcf.gz --data-type cnv

# for SV
bash VarXOmics_cnvsv.sh --input-sample P001_167 -hpo HP:0003002 -o sv -v examples/sv/P001_167.sv.vcf.gz --data-type sv
```

Input: `P001_167.cnv.vcf.gz`

### UI design

#### Variant Summary 
  - Table: `P001_167.cnv.annotated.tsv`
  - Circos figure: `P001_167_vcf2circos.html`
  - Decipher

#### Exomiser Prioritization
`P001_167_Exomiser_PHENIX_PRIORITY.json`
