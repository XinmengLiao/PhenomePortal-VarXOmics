This is the showcase for single variant query. Only the files will used for webserver are listed here. Missing files means no data generated.

Command for query: 
```
sh SingleVariant.sh chr13_32338103_G_GA \
	/mnt/storage_pool/Genomics/VarXOmics/examples/chr13_32338103_G_GA/chr13_32338103_G_GA.vcf.gz \
	/mnt/storage_pool/Genomics/VarXOmics/examples/chr13_32338103_G_GA \
	GRCH38 yes Female
```

Input: `chr13_32338103_G_GA.vcf.gz`

#### Variant summary
  - Table: `chr13_32338103_G_GA.txt`
  - Lollipop: `Lollipop.html`
  - AF world map: `AFMap.png`
  - PubMed

#### xQTL, GWAS, MR, and PGx
  - GWAS Catalog: `chr13_32338103_G_GA.gwas.txt`
  - eqtl Catalog: `chr13_32338103_G_GA.eqtl_catalog.txt`
  - eQTL GTEx: `chr13_32338103_G_GA.eqtl_gtex.txt`
  - pqtl EpigraphDB: `chr13_32338103_G_GA-pqtl_epigraphdb.txt`
  - PGx: `chr13_32338103_G_GA.pgx.txt`
  - SNP MR EpigraphDB: `chr13_32338103_G_GA.snpmr.txt`


### UI design
It's good to split the page into two side bars: variant summary and xQTL, GWAS, MR, and PGx. \
In the Variant summary page, I hope more annotations can be splitted into different sections for more straightforward display. 

#### Basic Information
columns: HGVSc, HGVSp, VARIANT_CLASS, Consequence, IMPACT, EXON, ClinVar_CLNSIG, acmg_classification, Existing_variation, MAX_AF, MAX_AF_POPS
mutated protein figure: `Lollipop.html`

#### Clinical Pathogenicities
 - ClinVar columns: ClinVar, ClinVar_CLNSIG, ClinVar_CLNDN,  ClinVar_CLNHGVS ClinVar_CLNSIGINCL, ClinVar_CLNVC,  ClinVar_GENEINFO, ClinVar_CLNDISDB, ClinVar_CLNSIGCONF, ClinVar_CLNREVSTAT, ClinVar_CLNDNINCL
 - ACMG columns: acmg_score, acmg_classification, acmg_criteria

#### Prediction Score
 - Missense variant columns: ada_score, rf_score, REVEL, REVEL_rankscore, REVEL_score, AlphaMissense_class,  AlphaMissense_genome,  AlphaMissense_pathogenicity, AlphaMissense_protein_variant, AlphaMissense_transcript_id, AlphaMissense_uniprot_id
 - Splicing variant columns: SpliceAI_pred_DP_AG, SpliceAI_pred_DP_AL, SpliceAI_pred_DP_DG, SpliceAI_pred_DP_DL, SpliceAI_pred_DS_AG, SpliceAI_pred_DS_AL, SpliceAI_pred_DS_DG, SpliceAI_pred_DS_DL

#### Allele Frequency
 - 1000 Genomes Project: AF
 - gnomAD exonome columns: gnomADe_AF, gnomADe_AFR_AF, gnomADe_AMR_AF, gnomADe_ASJ_AF, gnomADe_EAS_AF, gnomADe_FIN_AF, gnomADe_MID_AF, gnomADe_NFE_AF, gnomADe_REMAINING_AF, gnomADe_SAS_AF
 - gnomAD genome columns: gnomADg_AF, gnomADg_AFR_AF, gnomADg_AMI_AF, gnomADg_AMR_AF, gnomADg_ASJ_AF, gnomADg_EAS_AF, gnomADg_FIN_AF, gnomADg_MID_AF, gnomADg_NFE_AF, gnomADg_REMAINING_AF, gnomADg_SAS_AF
 - Map figure: `AFMap.png`

examples: \
![ALT TEXT](https://github.com/XinmengLiao/APMI-VarXOmics/blob/main/images/SingleVariant-VariantSummary.png)

#### xQTL, GWAS, MR, and PGx
Just keep the current design of this page with different bar of **eqtl Catalog, eQTL GTEx, pqtl, GWAS, MR, PGx**.
