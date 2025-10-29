## 🔹 Multi-Variant Analysis

VarXOmics supports **multi-variant analysis** for **small variants (SNPs & INDELs)**, **copy number variants (CNVs)**, and **structural variants (SVs)**.  
This module allows users to upload VCF files, configure advanced filtering parameters, and perform variant interpretation in a structured and reproducible manner.

### 1. Input Data

- From the **Home** page, choose **“Multi Variant”**.  
- Drop and upload a **VCF file**.  
- Click **Run** to proceed to the variant filtering configuration window.


### 2. Variant Filtering Settings

After uploading the VCF, a configuration panel titled **“Multi-Variant Analysis Settings”** will appear.  
This panel allows fine-tuning of variant filtration and prioritization parameters.


### 🧩 Sample Filters

| Setting | Description |
|----------|-------------|
| **Reference Genome** | Select the genome build used for alignment — **GRCh37** or **GRCh38**. |
| **Gender** | Choose the biological sex of the sample (**Male**, **Female**, or **Unknown**) for filtering X/Y-linked variants. |
| **Only PASS Variants** | When set to **Yes**, filters out non-PASS variants from the VCF. |
| **Inheritance Pattern** | Specify the expected mode of inheritance (**Dominant**, **Recessive**, **Compound Heterozygous**, or **Any pattern**). |
| **HPO Terms (CSV)** | Provide one or multiple **HPO IDs** (comma-separated) to integrate phenotype-based filtering, e.g. `HP:0001250,HP:0000707`. |

### 🌍 Allele Frequency Filters

| Setting | Description |
|----------|-------------|
| **AF ClinVar** | Sets a frequency cap for ClinVar-derived evidence. Variants above this threshold are excluded. |
| **AF (Predicted variants)** | Defines the **allele frequency cutoff** for predicted variants (e.g., 0.05). Lower thresholds retain rarer variants. |

### 🤖 Predictive Scores

| Score | Description | Typical Cutoff |
|--------|--------------|----------------|
| **ADA** | ADA score for splicing impact prediction. | 0.6 |
| **RF** | Random Forest model score for splicing variant prediction. | 0.6 |
| **REVEL** | REVEL score for missense variant pathogenicity prediction. | 0.75 |

### 🧬 SpliceAI Thresholds

| Parameter | Description | Default |
|------------|--------------|----------|
| **Acceptor Loss (AL)** | Minimum threshold for acceptor site loss prediction. | 0.5 |
| **Acceptor Gain (AG)** | Minimum threshold for acceptor site gain prediction. | 0.5 |
| **Donor Loss (DL)** | Minimum threshold for donor site loss prediction. | 0.5 |
| **Donor Gain (DG)** | Minimum threshold for donor site gain prediction. | 0.5 |

### 🧮 BayesDel Scores

| Parameter | Description | Default |
|------------|--------------|----------|
| **BayesDel (with AF)** | Bayesian deleteriousness score integrating allele frequency. | 0.069 |
| **BayesDel (no AF)** | Bayesian deleteriousness score excluding frequency adjustment. | -0.057 |

> These thresholds align with backend default values for consistent filtering.


### 🧾 Classification Filters

#### **AlphaMissense (AM) Classification**
Select which **ACMG/AMP pathogenicity tiers** to retain in the analysis:
- Pathogenic  
- Likely Pathogenic  
- Ambiguous  
- Benign  

#### **AlphaMissense (AM) Pathogenicity**
- Numerical cutoff representing aggregated pathogenicity evidence (e.g., 0.564).

#### **ClinVar Status**
Select ClinVar review statuses to include:
- Pathogenic  
- Likely Pathogenic  
- Uncertain Significance  
- Conflicting Classifications  
- Benign  
- Likely Benign  

#### **ACMG Classification**
Specify ACMG evidence tiers for filtering:
- Pathogenic  
- Likely Pathogenic  
- Uncertain Significance  
- Benign  
- Likely Benign  

### 3. Run Analysis

Once all filters are configured, click **Start Analysis** to execute the variant prioritization pipeline.  


### 4. Result Overview
VarXOmics processes and prioritizes all variants based on its integrated filtering and scoring strategy. The final results are divided into six major sections, each providing an in-depth view of the analysis outcome.

---

#### 🧬 Variant Overview

This section displays all filtered variants prioritized by the VarXOmics strategy. Each variant is annotated comprehensively, including genomic, transcript, clinical, and functional evidence.

**Key Features:**
- **Detailed Annotation Table:**  
  Displays complete variant-level information.  
  Users can customize visible columns using the **“Columns”** button.
- **Quick Summary Subpanels:**  
  Provides concise visual summaries of:
  - **Variant Consequence Distribution**  
  - **Variant Types**  
  - **Clinical Pathogenicity**  
  - **Affected Genes**

These visual overviews allow users to rapidly assess the biological and clinical relevance of the prioritized variants.

---

#### 📊 Variant Summary

This section provides an integrated summary of all filtered variants through seven visual plots, allowing users to explore the overall variant landscape.

**Displayed Figures:**
  - **SNP Density Plot:**  
   Visualizes the distribution of SNPs across all chromosomes.
  - **ClinVar Ideogram:**  
   Shows the chromosomal positions of ClinVar-annotated variants and their pathogenicity status.
  - **ClinVar Pathogenicity:**  
   Summarizes the proportion of variants categorized by ClinVar pathogenicity.
  - **Variant Class:**  
   Displays the proportion of SNPs, INDELs, and other variant types.
  - **Top 10 Disease Accumulators:**  
   Highlights the diseases most frequently associated with the filtered variants.
  - **Top 10 Gene Accumulators:**  
   Shows the genes where variants are most enriched.
  - **Variant Consequence Bar Plot:**  
   Displays the frequency of different variant consequences (e.g., missense, nonsense, synonymous).

Together, these plots offer an intuitive view of the biological and clinical impact of the variant set.

---

#### 🔬 Enrichment Analysis

This section performs functional enrichment for all genes linked to the filtered variants, helping users interpret the biological processes and pathways involved.

**Components:**
- **GO Term Enrichment Table:**  
  Lists significantly enriched **Gene Ontology (GO)** terms (Biological Process, Molecular Function, Cellular Component).
- **KEGG Pathway Table:**  
  Shows **KEGG pathway enrichment** results for the filtered gene set.
- **GO Enrichment Visualization:**  
  - Top enriched GO terms ranked by gene count.
  - GO Term Counts.
  - Statistics of GO Biological Processes, Molecular Functions, Cellular Components. 
- **KEGG Enrichment Visualization:**  
  Presents top pathways ranked by enrichment significance.

The enrichment section provides direct insight into the biological themes underlying the filtered variant set.

---

#### 📚 GWAS & QTL Associations

This section presents variant- and gene-level associations collected from multiple public resources, allowing users to explore known evidence linking the queried variants or genes to complex traits, diseases, and regulatory effects.  
Contents cover: GWAS, eQTL, pQTL, Pharmacogenomics (PGx).  
Users can customize visible columns with the **“Columns”** button (top-right corner).  
Together, these resources allow users to assess the translational and regulatory significance of their filtered variants.

---

#### 🕸️  Network

This section visualizes all biological interactions derived from the filtered variant and gene set.  
This integrative network provides a global view of how variants, genes, proteins, drugs, and phenotypes interconnect through multiple evidence layers.  
This network provides an integrated, evidence-weighted visualization of how variants influence genes, pathways, and drug responses — bridging molecular biology and clinical relevance.


**Interactive Features:**
- **Edge Weight Filter:**  
  Adjust the slider to dynamically filter edges weights. For each interactioins, the weights are set as:
  
| Interaction Type | Weight / Significance Cutoff |
|------------------|------------------------------|
| **PPIs** | Interaction confidence **> 0.7** |
| **GWAS** | **p < 5 × 10⁻⁸** |
| **eQTL / pQTL** | **p_beta < 0.05** |
| **MR (Mendelian Randomization)** | **p < 0.05** |

- **Search Box:**  
  Locate nodes by label (e.g., gene symbol, variant ID, drug name).
- **Node Summary Panel:**  
  Displays total counts of nodes and edges, categorized by node and interaction type.
- **Table View:**  
  Click the **table icon** in the bottom-right corner to open a structured view of the network, allowing users to inspect and filter interaction data directly.

---

✅ **Tip:**  
Users can combine interaction filters and weight thresholds to create **focused sub-networks**, such as *gene–drug pharmacogenomic networks* or *variant–phenotype associations* for specific pathways.


#### 🧠 Exomiser

VarXOmics integrates **Exomiser**(https://exomiser.readthedocs.io/en/stable/#) analysis to support phenotype-driven variant prioritization. By linking variant annotations with **HPO terms**, Exomiser scores potential candidate variants and genes according to their phenotypic relevance.

This section presents:  
- Ranked list of candidate genes and variants associated with input HPO phenotypes.  
- Exomiser scores indicating phenotype similarity and variant pathogenicity.  
- Cross-references with ClinVar and other annotations for interpretation consistency.

---

### 5. Export results  
Users can export all results, figures, and enrichment tables as a single compressed archive by clicking **“Export All”** in the top-right corner of the interface.
