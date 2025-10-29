**VarXOmics** is a versitile web server for genomic data querying and variants analysis and priortisation. The data sources are integrated from diverse open-source databases, including multi-omics informations (eQTL, pQTL, GWAS, Traits), pharmacogenomics infomration, protein-protein interactions, gene-drug interactions. The detailed data could be accessed as described followed:


### Multi-Omics insights

#### 1) eQTL

GTEx eQTL: [GTEx v10-eQTL](https://www.gtexportal.org/home/downloads/adult-gtex/qtl).

eQTL Catalog: [eQTL/sumstats](https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/).

Pre-calculated scores: `Data Access/pre-calculated eQTL scores.txt` Scores are calculated only with significant eQTL results.

#### 2) pQTL and MR

[EpiGraphDB](https://github.com/MRCIEU/epigraphdb-r) Data are retrieved by `epigraphdb-r`.

Pre-calculated scores: `Data Access/pre-calculated pQTL scores.txt` Scores are calculated by considering significant pQTL results and high-confidence protein-protein coexpression interactions.

#### 3) GWAS

GWAS Catalog: [GWAS Catalog](https://www.ebi.ac.uk/gwas/docs/file-downloads) (Version: 1.0-associations_e113_r2025-04-28).

#### 4) Tratis

Curated trait list: `Data Access/Trait lists.txt` The list contains summary statistic for Athleticism, Behavior, Cardiovascular, Hormones, Immune System, Longevity, Metabolism, Nutrition and Diet, Physical Appearance, Sensory Perception, Wellness, Neurogenic and Cognitive functions, and Personality traits.

### Pharmacogenomics

PharmGKB: [PharmGKB](https://www.pharmgkb.org/downloads) (Version: 2025-May-05).

PharmGKB-pediatric drug information: [PharmGKB](https://www.pharmgkb.org/downloads) (Version: 2024-May-15).

PharmVar: [PharmVar](https://www.pharmvar.org/download) (Version: 6.1.2).

### Gene-Drug Interactions

Gene-Drug interaction database: [DGIdb](https://dgidb.org/downloads) (Version: 2025-May).

### Protein-Protein Interactions

STRING database: [STRING-homo spaiens](https://string-db.org/cgi/download?sessionId=bMqCTjlbF9J3&species_text=Homo+sapiens&settings_expanded=0&min_download_score=0&filter_redundant_pairs=0&delimiter_type=txt) (Version: 12.0).

### Gene expression levels in different tissues

Human Protein Atlas: [HPA](https://www.proteinatlas.org/), which includes `RNA expression (consensus) (50 tissues)` and `RNA expression (HPA) (40 tissues)`.

GTEx v10: [GTEx v10-RNA bulk_tissue_expressionl](https://www.gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression).

### Gene-Disease list

Curated gene-disease list: `Data Access/Gene-Disease panel.xlsx` It consists of lists from ACMG secondary findings v3.2, lists for carrier screening and newborn screening.

### Software requirements:
Small varinat querying and analysis is conducted using [Ensembl Variant Effect Predictor (VEP)](https://www.ensembl.org/info/docs/tools/vep/index.html) (e114), while structural variants are analysed by [AnnotSV](https://github.com/lgmgeo/AnnotSV) (v3.4). The pipeline for data management and scoring sytem is built on `R(v4.3.3)` and and `python (v3.12)`.
