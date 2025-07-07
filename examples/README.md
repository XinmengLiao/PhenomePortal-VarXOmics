Here are the examples showcasing 4 functions of VarXOmics. 

1. Single variant query: `chr13_32338103_G_GA`
2. Single gene query: `APOB`
3. Multiple variant analysis (small variants): `P00110_11/`
4. Multiple variant analysis (SVs and CNVs): `cnvsv`


### UI design for the home page
The structure looks overall really nice. I hope some functions could be added. 
1. Change the description to: 
2. vcf file could only be vcf.gz file currently.
3. More filter options could be added. These filterations will be applied only for multiple variant anlaysis. \

   **Clinical Pathogenicity (multi-Selections)** \
       a) ClinVar Pathogenicity (6 options): Pathogenic, Likely pathogenic, Uncertain Significance, Conflicting classifications, Benign, Likely benign \
       b) ACMG Predicted Classification (5 options):  Pathogenic, Likely pathogenic, Uncertain Significance, Benign, Likely benign

      **Allele Frequency (float)** \
      a) ClinVar variant: default display: 1 \
      b) Predicted variant: default display: 0.05 \
      
      **Variant Deleteriouseness** \
      a) PolyPhen Score: default display: 0.05 \
      b) BayesDel Score (addAF): \
      c) BayesDel Score (noAF): \
      
      **Missense Impact** \
      a) AlphaMissense Classification (multi-selection, 3 options): Likely pathogenic, Uncertain Significance, Likely benign \
      b) AlphaMissense Score: default display: 0.564 \
      c) REVEL Score: default display: 0.5 \
      d) SIFT Score: default display: 0.05 \
      f) Ada Score: \
      e) Rf score:  \
      
      **Splicing Impact** \
      a) Acceptor Gain :default display: 0.8 \
      b) Acceptor Loss: default display: 0.8 \
      c) Donor Gain: default display: 0.8 \
      d) Donor Loss: default display: 0.8 \
   
   All these filterations will be passed to configuration files `Scripts/multivariant_config.yaml`

5. Descriptions change to:
    Multi-database variant annotation \
    Multi-omics variant evidence \
    Functional enrichment analysis \
    Biological network analysis \
    Variant evidence aggregated scoring \
    Variant-Phenotype significance ranking \
    Pharmacogenomic associations \
    Custom selecting and data exporting \

6. Add more options:
   a) Reference genome version: GRCH37 or GRCH38
   b) 'PASS only' or not
   
