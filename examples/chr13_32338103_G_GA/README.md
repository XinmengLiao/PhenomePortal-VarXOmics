This is the showcase for single variant query.

Command for query: 
```
varxomics='/mnt/storage_pool/Genomics/VarXOmics'
bash $varxomics/SingleVariant.sh \
	--id [FILENAME] \
	-i [INPUT_VCF] \
	-o [OUTPUT_DIR] \
	--genome [GENOME] --only-pass [yes|no] --gender [Male|Female]

# example
bash $varxomics/SingleVariant.sh --id chr13_32338103_G_GA \
	-i /mnt/storage_pool/Genomics/VarXOmics/examples/chr13_32338103_G_GA/chr13_32338103_G_GA.vcf.gz \
	-o /mnt/storage_pool/Genomics/VarXOmics/examples/chr13_32338103_G_GA \
	--genome GRCH38 --only-pass yes --gender Female
```

Input: `/mnt/storage_pool/Genomics/VarXOmics/examples/chr13_32338103_G_GA/chr13_32338103_G_GA.vcf.gz` \
*Input must be in compressed **vcf.gz** format

The following outputs could be found in: `/mnt/storage_pool/Genomics/VarXOmics/examples/chr13_32338103_G_GA`

#### Variant summary
  - Table: `chr13_32338103_G_GA.txt`
  - Lollipop: JSON file: `Lollipop.json` and HTML file: `Lollipop.html`
  - AF world map: `AFMap.png`
  - PubMed

#### xQTL, GWAS, MR, and PGx
  - GWAS Catalog: `chr13_32338103_G_GA.gwas.txt`
  - SNP MR EpigraphDB: `chr13_32338103_G_GA.snpmr.txt`
  - eqtl Catalog: `chr13_32338103_G_GA.eqtl_catalog.txt`
  - eQTL GTEx: `chr13_32338103_G_GA.eqtl_gtex.txt`
  - pqtl EpigraphDB: `chr13_32338103_G_GA.pqtl-epigraphdb.txt`
  - pqtl Open Targets Genomics: `chr13_32338103_G_GA.pqtl_otg.txt`
  - PGx: `chr13_32338103_G_GA.pgx.txt`


### UI design
It's good to split the page into two panels: **variant summary** and **xQTL, GWAS, MR, and PGx**. \
In the Variant summary page, I hope annotations can be splitted into different sections for more straightforward display as following: 

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


### UI Examples: 
#### Variant Summary
![ALT TEXT](https://github.com/XinmengLiao/APMI-VarXOmics/blob/main/images/SingleVariant-VariantSummary.png)

#### xQTL, GWAS, MR, and PGx
Just keep the current design of this page with different side bars and show the tables of **eqtl Catalog, eQTL GTEx, pqtl, GWAS, MR, PGx**.
