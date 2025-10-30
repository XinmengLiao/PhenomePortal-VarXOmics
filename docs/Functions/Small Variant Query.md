## 🔹 Small Variant Query

### 1. Input Variant  
- From the **Home** page, choose **“Single Variant”**.  
- Enter the variant location using the format:  `chromosome_position_reference allele_alternative allele` \
    Example: `chr1_123456_A_T`  
- Click **Run** to proceed.  

<img width="40%" src="https://github.com/user-attachments/assets/c8833545-a799-4584-b172-02afab10f9a9" />

<br>

### 2. Select Reference Genome  
- Choose the reference genome: **GRCh37** or **GRCh38**.  
- Click **Start Analysis** to begin.  
> **Note:** Variants aligned with **GRCh37** will be **lifted over to GRCh38** before VEP (v114) annotation.

<img width="60%" src="https://github.com/user-attachments/assets/64558385-5fad-48ed-98c1-e0ba8fb5d10b" />

<br>

### 3. View Results  
Results will appear automatically once the data is loaded. Annotations are organized into four sections:

#### 🧬 Variant Annotation
- Contents includes:
  - Transcript and protein alterations  
  - ClinVar information  
  - Global allele frequencies (gnomAD, 1000 Genomes Project)  
  - In silico prediction scores  
- Users could customize displayed columns using the **“Columns”** button.

<img width="80%" src="https://github.com/user-attachments/assets/72278cf6-74e9-4bb5-b89b-4c01f94e0c7a" />


#### 🍭 ClinVar Lollipop
- Visualizes **protein alterations** as a lollipop plot.  
- Supports **zooming** and **export** in **PNG** or **PDF** format.

<img width="80%" src="https://github.com/user-attachments/assets/748bf42a-76b6-4e91-9fdf-b0a59105a5d9" />

<br>

#### 📊 GWAS & QTLs
- Displays variant-associated data in separated subpanels:
  - **GWAS**
  - **eQTLs**
  - **pQTLs**
  - **Pharmacogenomics (PGx)**

<img width="80%" src="https://github.com/user-attachments/assets/a04ad5d4-aff5-490b-a067-dbe4d25a20a3" />


#### 🌍 Global Allele Frequency
- Shows the **maximum global minor allele frequency (MAF)** on a **world map**, illustrating population-specific frequency distributions.

<img width="80%" src="https://github.com/user-attachments/assets/2cbfb444-c88b-4084-95e7-183b1af03803" />

<br>

### 4. Export Data  
- Download **all results** (tables and figures) as a **compressed file** via the **“Export All”** button (top-right corner).
