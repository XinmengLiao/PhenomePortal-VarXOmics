# VarXOmics Variant Query Guide

VarXOmics supports querying **small variants (SNPs & INDELs)**, **copy number variants (CNVs)**, and **structural variants (SVs)**. Both **GRCh37** and **GRCh38** reference genomes are supported.

## 🔹 Small Variant Query

### 1. Input Variant  
- From the **Home** page, choose **“Single Variant”**.  
- Enter the variant location using the format:  `chromosome_position_reference allele_alternative allele` \
    Example: `chr1_123456_A_T`  
- Click **Run** to proceed.

### 2. Select Reference Genome  
- Choose the reference genome: **GRCh37** or **GRCh38**.  
- Click **Start Analysis** to begin.  
> **Note:** Variants aligned with **GRCh37** will be **lifted over to GRCh38** before VEP (v114) annotation.

### 3. View Results  
Results will appear automatically once the data is loaded. Annotations are organized into four sections:

#### 🧬 Variant Annotation
- Contents includes:
  - Transcript and protein alterations  
  - **ClinVar** information  
  - **Global allele frequencies** (gnomAD, 1000 Genomes Project)  
  - **In silico** prediction scores  
- Users could customize displayed columns using the **“Columns”** button.

#### 🍭 ClinVar Lollipop
- Visualizes **protein alterations** as a lollipop plot.  
- Supports **zooming** and **export** in **PNG** or **PDF** format.

#### 📊 GWAS & QTLs
- Displays variant-associated data in structured tables:
  - **GWAS**
  - **eQTLs**
  - **pQTLs**
  - **Pharmacogenomics (PGx)**

#### 🌍 Global Allele Frequency
- Shows the **maximum global minor allele frequency (MAF)** on a **world map**, illustrating population-specific frequency distributions.

### 4. Export Data  
- Download **all results** (tables and figures) as a **compressed file** via the **“Export All”** button (top-right corner).


## 🔹 CNV & SV Query

### 1. Input Variant  
- From the **Home** page, choose **“Single CNV/SV”**.  
- Enter the variant in the format: `chromosome_starting position_ending position_mutation type`

  Example: `chr3_100000_200000_dup`  
- Click **Run** to proceed.

### 2. Select Reference Genome and Type  
- Reference genome options: **GRCh37**, **GRCh38**  
- Mutation type options: **CNV**, **SV**  
- Click **Start Analysis** to begin.  
> **Note:** CNVs and SVs are annotated using **AnnotSV**.

### 3. View Results  
- Annotated results appear once loaded.  
- Customize displayed columns using the **“Columns”** button.

### 4. Export Data  
- Export all annotation tables as a **compressed file** via the **“Export All”** button.
