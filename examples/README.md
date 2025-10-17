Here are the examples showcasing 4 functions of VarXOmics. 

1. Single small variant searching: `chr13_32338103_G_GA`
2. Single gene query: `BRCA2`
3. Multiple variant analysis (small variants): `P00110_11/`
4. Multiple variant analysis (SVs and CNVs): `P001_167.cnv.gz` and `P001_167.sv.gz`
5. Single CNV/SV searching:
   - CNV: `10_17268559_17274050_DEL`
   - SV: `11_65736014_65736014_INS`

## 2025-October-17
### Home page
<img width="1902" height="952" alt="image" src="https://github.com/user-attachments/assets/0ca85db1-56e8-4a48-a647-2c058a2aa8ec" />

1. Descriptions change to:
Title: Features:
    - Comprehensive variant annotation: HGVSc, HGVSp, consequences on proteins, clinical evidence, global allele frequency, prediction scores.
    - Multi-omics insights: insights of eQTL, pQTL, GWAS, MR, and PGx.
    - Diverse functions: functional pathway enrichment, biological network interactions.
    - Variant prioritizations: VarXOmics and Exomiser approaches for stratifying significant and/or disease-causing variants.
    - Visualizations: Downable publication-ready figures and structured tabular data.
2. The "Help & Support" could link to https://github.com/XinmengLiao/PhenomePortal-VarXOmics/tree/main/Data%20Access

### Single gene searching
<img width="1830" height="494" alt="image" src="https://github.com/user-attachments/assets/cb864a5f-4057-4e58-a7e4-6a37c388dc04" />

1. Change the title from "Gene Evidence Summary" to "Clinical Evidence"
2. Remove "Expression & Links" section since there are own three link-outs. Could you move these link-outs to somewhere maybe under the Gene Tag?
3. Columns of ClinGen.Disease (ClinGen), Diseases (ClinVar), GenCC.disease_title (GenCC) could be used for showing the disease. Particularly, Diseases (ClinVar) contains multiple diseases in one line, so need to seperate them by '|'

### Single variant searching 
<img width="1812" height="412" alt="image" src="https://github.com/user-attachments/assets/8dc80d08-b797-4fd1-babf-09ad40b676b6" />

1. I tried to search chr13_32338103_G_GA (JobID: ID: 32180abf), but "Variant Annotation" and "GWAS & QTLs" sections are empty. There should be some results. 
2. Change "Figures" to "Global Allele Frequency", and rescale the figure. 


### Multiple small variant anlaysis
<img width="920" height="933" alt="image" src="https://github.com/user-attachments/assets/f5b2b07a-bdcb-42e6-afb8-30598673351c" />

1. Remove "Mitochondrial" and "Compound heterzygous" in the inheritance filtering options, since they are not inheritance.
2. Move GWAS and PGx results into the "GWAS & QTLs & PGx" side bar.
3. Change while to black text in the coloumn filtration. 
4. Remove the "REMM Genes" and "REMM Variants" sections. Don't need to show these two results in the interface.

### Single CNV/SV searching and Multiple CNV/SV analysis
Add DECIPHER searching function inside the result. Directly show the exact location of the searching regions, but users could also scroll around the genome. 
