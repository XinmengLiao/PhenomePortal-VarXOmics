## 🔹 Multi SV/CNV Analysis

VarXOmics supports **multi-sample CNV (Copy Number Variant)** and **SV (Structural Variant)** analysis.  

### 1. Input Data

- From the **Home** page, choose **“Multi SV/CNV”**.  
- Drop to upload **CNV/SV VCF file** (e.g., `sampleID.cnv.vcf.gz`).  
- Click "Run" and the analysis settings will pop-up:
  - Reference Genome (**GRCh37** or **GRCh38**)  
  - Variant Type (**CNV** or **SV**)  
  - HPO ID, separated by comma. (e.g., HP:0031218,HP:0002071) 
- After configuration, click **Start Analysis** to begin.

<img width="60%" src="https://github.com/user-attachments/assets/a3a31be2-ef6c-43f0-9f6d-05c6d8871259" />

<br>

### 2. Results Overview

#### 🧾 Annotated TSV Table

Displays all **annotated CNVs/SVs** identified across the uploaded samples. CNVs and SVs are annotated by AnnotSV. Users can choose which columns to display via the **“Columns”** button (top-right corner).  

#### 🧬 Circos Visualization
This section presents a **genome-wide Circos plot** of all identified CNVs and SVs.

### 3. Export
Users can export all results, figures, and enrichment tables as a single compressed archive by clicking **“Export All”** in the top-right corner of the interface.
