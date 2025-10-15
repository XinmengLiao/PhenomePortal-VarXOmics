This is the showcase for multiple variant analysis, with breast cancer patient WGS data.

Input files:
vcf: `P00110_11.vcf.gz` \

Command for analysis
```bash
varxomics='/mnt/nas/Genomics/VarXOmics'
bash VarXOmics_snv.sh \
  -i test \
  -o $varxomics/examples/test \
  -v $varxomics/examples/test/test.vcf.gz \
  --genome GRCH38 --only-pass yes \
  --gender Female --inheritance AD \
  --hpo HP:0003002,HP:0006625
```

#### Variant Summary section:
  - Basic information: `P00110_11.priorisation.txt`
  - Variant accumulated Genes: `Figures_data/user_selected_GeneAccumultate-all.txt` and figure `Figures/user_selected_GeneAccumultate-AllorTop10.pdf`
  - Variant accumulated ClinVar disease: `Figures_data/user_selected_ClinVarPathogenicity.txt` and figure `Figures/user_selected_ClinVarPathogenicity.pdf`
  - Variants type: `Figures_data/user_selected_VariantClass.txt` and figure `Figures/user_selected_VariantClass.pdf`
  - Variants consequence: `Figures_data/user_selected_VariantConsequence.txt` and figure `Figures/user_selected_VariantConsequence.pdf`
  - Variant SNP density: `Figures_data/user_selected_SNP_Density.txt` and figure `Figures/user_selected_SNP_density.jpg`
  - ClinVar Pathogenicity Distribution: `Figures_data/user_selected_ClinVarPathogenicity.txt` and figure `Figures/user_selected_ClinVarPathogenicity.pdf`
  - ClinVar variants on all chromosomes: `Figures/ideogram.png`

#### GO and KEGG: (skip all figures for the current version)
  - Gene-GO terms (need to split by BP, MF, CC types): `GO-All.txt`
  - GO terms count (need to split by BP, MF, CC types): `GO-All-count.txt`
  - GO Enrichment results (need to split by BP, MF, CC types): `GO-All-Enrich.txt`
  - GO BP enrichment figures: `Figures/BPEnrich-AllorTop10.png` and `Figures/BPEnrichCircus.pdf`
  - GO MF enrichment figures: `Figures/MFEnrich-AllorTop10.png` and `Figures/MFEnrichCircus.pdf`
  - GO CC enrichment figures: `Figures/CCEnrich-AllorTop10.png` and `Figures/CCEnrichCircus.pdf`
  - Gene-KEGG terms: `KEGG-All.txt`
  - KEGG terms count: `KEGG-All-count.txt`
  - KEGG Enrichment results: `KEGG-All-Enrich.txt`
  - KEGG enrichment figures: `Figures/KEGGEnrich-AllorTop10.png` and `Figures/KEGGEnrichCircus.pdf`

#### Network:
  `P00110_11.network.txt`

#### GWAS, xQTL, and PGx
  - GWAS: `P00110_11.gwas_filtered.txt` and `P00110_11.gwas.txt`
  - eQTL Catalog: `P00110_11.eqtl_catalog_filtered.txt` and `P00110_11.eqtl_catalog.txt` 
  - eQTL GTEx: `P00110_11.eqtl_gtex_filtered.txt` and `P00110_11.eqtl_gtex.txt`
  - pQTL Open Target Genetics: `P00110_11.pqtl_otg_filtered.txt` and `P00110_11.pqtl_otg.txt`
  - pQTL EpigraphDB: `P00110_11.pqtl_epigraphdb_filtered.txt`
  - PGx: `P00110_11.pgx_filtered.txt` and `P00110_11.pgx.txt`

#### Original Exomiser priortization 
If user requires to run Exomiser original priortization, then also show another side bar of: \
**Exomiser PHENIX_PRIORITY** `P00110_11_Exomiser_PHENIX_PRIORITY.json`


 #### Interactive function
 Users can selected their variants of interest, and all the sections only shows results of their selections. 
 1. Variant Summary section: the subset of users selection.
 2. Rerun script `Scripts/GO-KEGG-subset.R` and figures regenerated in:
  - Gene-GO terms (need to split by BP, MF, CC types): `user_selected_GO-All.txt`
  - GO terms count (need to split by BP, MF, CC types): `user_selected_GO-All-count.txt`
  - GO Enrichment results (need to split by BP, MF, CC types): `user_selected_GO-All-Enrich.txt`
  - GO BP enrichment figures: `Figures/user_selected_BPEnrich.png` and `Figures/user_selected_BPEnrichCircus.pdf`
  - GO MF enrichment figures: `Figures/user_selected_MFEnrich.png` and `Figures/user_selected_MFEnrichCircus.pdf`
  - GO CC enrichment figures: `Figures/user_selected_CCEnrich.png` and `Figures/user_selected_CCEnrichCircus.pdf`
  - Gene-KEGG terms: `user_selected_KEGG-All.txt`
  - KEGG terms count: `user_selected_KEGG-All-count.txt`
  - KEGG Enrichment results: `user_selected_KEGG-All-Enrich.txt`
  - KEGG enrichment figures: `Figures/user_selected_KEGGEnrich.png` and `Figures/user_selected_KEGGEnrichCircus.pdf`
 3. Rerun script `Scripts/GO-KEGG-subset.R` and figures are regenerated in:
  - Variant accumulated Genes: `Figures_data/user_selected_GeneAccumultate-all.txt` and figure `Figures/user_selected_GeneAccumultate-AllorTop10.pdf`
  - Variant accumulated ClinVar disease: `Figures_data/user_selected_ClinVarPathogenicity.txt` and figure `Figures/user_selected_ClinVarPathogenicity.pdf`
  - Variants type: `Figures_data/user_selected_VariantClass.txt` and figure `Figures/user_selected_VariantClass.pdf`
  - Variants consequence: `Figures_data/user_selected_VariantConsequence.txt` and figure `Figures/user_selected_VariantConsequence.pdf`
  - Variant SNP density: `Figures_data/user_selected_SNP_Density.txt` and figure `Figures/user_selected_SNP_density.jpg`
  - ClinVar Pathogenicity Distribution: `Figures_data/user_selected_ClinVarPathogenicity.txt` and figure `Figures/user_selected_ClinVarPathogenicity.pdf`
  - ClinVar variants on all chromosomes: `Figures/ideogram.png`
 4. Network: the subset of users selection.
 5. GWAS, xQTL, and PGx: the subset of users selection.
