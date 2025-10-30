library(dplyr)
library(tidyverse)
library(data.table)
library(epigraphdb)
library(httr)
library(jsonlite)
for (i in c("select","filter", "mutate","rename", "left_join", "slice")){
  conflicted::conflict_prefer(i, "dplyr")
}
rm(i)
conflicted::conflicts_prefer(stats::sd)
conflicted::conflicts_prefer(httr::content)
conflicted::conflicts_prefer(plotly::layout)

## ---- Load files ----
#user input
# Rscript x.R P00110_8 HP:0003002,HP:0006625 AD
args <- commandArgs(trailingOnly = TRUE)
sampleid <- args[1]
hpo_terms <- args[2]
inheritance <- args[3]
print(paste0("Inheritance: ",inheritance))
user_dir = paste0(args[4],"/")
path = "sysmed" # local or sysmed

# Defined files
#user_dir = "/Users/xinmengliao/Documents/Project/20250516_Webserver/Case_Study/BreastCancer/result/"
#sampleid <- "P0097_182"
#sampleid <- "P00110_11"
#hpo_terms = "HP:0001627"
#hpo_terms = "HP:0003002,HP:0006625"
#inheritance = "AD"
#inheritance = ""

#outputs from the python
exomiser_file <-paste0(user_dir, sampleid, "_Exomiser_REMM.variants.tsv")
pgx_file <- paste0(user_dir, sampleid, ".pgx_filtered.txt")
gwas_file <- paste0(user_dir, sampleid, ".gwas_filtered.txt")
eqtl_catalog_file <- paste0(user_dir, sampleid, ".eqtl_catalog_filtered.txt")
eqtl_gtex_file <- paste0(user_dir, sampleid, ".eqtl_gtex_filtered.txt")
otg_file <- paste0(user_dir, sampleid, ".pqtl_otg_filtered.txt")

if(path == "local"){
  eqtl.score_file <- "/Users/xinmengliao/Documents/Project/20250516_Webserver/Datasets/xQTL/eQTL.combined.score.txt"
  pqtl.score_file <- "/Users/xinmengliao/Documents/Project/20250516_Webserver/Datasets/xQTL/eQTL-pQTL.score.txt"
  dgi_file <- "/Users/xinmengliao/Documents/Project/20250516_Webserver/Datasets/DGIdb/interactions.tsv"
  string_file <- "/Users/xinmengliao/Documents/Project/20250516_Webserver/Datasets/xQTL/String/String.all.score.txt"
}else if (path == "sysmed"){
  eqtl.score_file <- "/mnt/storage_pool/Genomics/VarXOmics/Databases/eQTL.combined.score.txt"
  pqtl.score_file <- "/mnt/storage_pool/Genomics/VarXOmics/Databases/eQTL-pQTL.score.txt"
  dgi_file <- "/mnt/storage_pool/Genomics/VarXOmics/Databases/DGIdb_interactions.tsv"
  string_file <- "/mnt/storage_pool/Genomics/VarXOmics/Databases/String.all.score.txt"
}


# load files
print(paste0(user_dir,sampleid, ".txt"))
result <- read.csv(paste0(user_dir,sampleid, ".txt"),header = T,sep = "\t") %>% unique()
colnames(result) <- gsub("am_","AlphaMissense_", colnames(result))
eqtl.score <- read.csv(eqtl.score_file, header = T,sep = "\t")
pqtl.score <- read.csv(pqtl.score_file,header = T,sep = "\t")

## ---- pQTL evidence  ----
# pQTL MR, and information, and scores
# EpiGraphDB contains the MR resutls, only for table and network 
# pQTL scores only come from Open Target Genetics

cat("\nAdding pQTL Evidence\n")

pqtl.score.fun <- function(data){
  cols <- c('expID','outID','outID_mrbase','nsnp','method','beta','se','pvalue','rsID')
  pqtl_tpm_df = as_tibble(setNames(as.list(rep(NA, length(cols))), cols))
  
  rsid.list <- unique(data$rsID) 
  
  if(length(rsid.list) == 0){
    return(pqtl_tpm_df)
  }else{
    if(length(which(rsid.list == "")) >0){
      rsid.list <- rsid.list[- which(rsid.list =="")]
    }
    print(paste0("Totally ", length(rsid.list), " rsID are using for retrieving pQTL evidence."))
    
    for (i in 1:length(rsid.list)){
      if(i %% 500 == 0){
        print(paste0("Now processing ", i, "th rsID."))
      }
      rela.pro <- pqtl_pleio(rsid =  rsid.list[i])
      if(nrow(rela.pro) > 0){
        for (j in 1:nrow(rela.pro)){
          pqtl_tpm <- pqtl(query = rela.pro$expID[j], searchflag = "proteins",rtype = "mrres",pvalue = 1) 
          if(nrow(pqtl_tpm)> 0 ){
            pqtl_tpm <- pqtl_tpm %>% mutate(rsID = rsid.list[i])
            pqtl_tpm_df <- rbind(pqtl_tpm_df, pqtl_tpm)
          }
        }
      }
    }
    
    colnames(pqtl_tpm_df) <- paste0("EpiGraphDB-",colnames(pqtl_tpm_df) )
    pqtl_tpm_df <- pqtl_tpm_df %>% na.omit(`EpiGraphDB-expID`) %>% 
      mutate(`EpiGraphDB-rsID`= as.character(`EpiGraphDB-rsID`))

    return(pqtl_tpm_df)
  }
  
}

pqtl.df <- pqtl.score.fun(result)
write.table(pqtl.df,file = paste0(user_dir, sampleid, ".pqtl_epigraphdb_filtered.txt"),quote = F,sep = "\t",row.names = F)



## ---- Scoring functions ----

# Normalize score
# normalize_minmax <- function(x) {
#   if (length(unique(x)) == 1) {
#     return(rep(0, length(x)))
#   }
#   (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
# }

normalize_minmax <- function(x, epsilon = 1e-6) {
  if (length(unique(x)) == 1) {
    return(rep(0.5, length(x)))  # 如果所有值都相同，返回0.5中间值
  }
  res <- (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  res <- res * (1 - epsilon) + epsilon
  return(res)
}


# ClinVar and ACMG scores and AF scores
patho.score <- function(data){
  min_af_nonzero <- min(data$MAX_AF[data$MAX_AF > 0], na.rm = TRUE)
  
  final.data <- data %>% 
    mutate(ClinVar.score = case_when(grepl("^Pathogenic", fixed = F, ClinVar_CLNSIG) ~ 1,
                                     grepl("^Likely_pathogenic", fixed = F, ClinVar_CLNSIG) ~ 0.8,
                                     (grepl("^Conflicting_classifications_of_pathogenicity",ClinVar_CLNSIG,fixed = F) & grepl("^Pathogenic", fixed = F, ClinVar_CLNSIGCONF)) ~ 0.6, 
                                     (grepl("^Conflicting_classifications_of_pathogenicity",ClinVar_CLNSIG,fixed = F) & !grepl("^Pathogenic", fixed = F, ClinVar_CLNSIGCONF)) ~ 0.5, 
                                     grepl("^Uncertain_significance", fixed = F, ClinVar_CLNSIG) ~ 0.5,
                                     grepl("^Benign", fixed = F, ClinVar_CLNSIG) ~ 0.2,
                                     grepl("^Likely_benign", fixed = F, ClinVar_CLNSIG) ~ 0.2, TRUE ~ 0.1), 
           ACMG.score = case_when(acmg_classification == "Pathogenic" ~ 1,
                                  acmg_classification == "Likely_pathogenic" ~ 0.8,
                                  acmg_classification == "Uncertain_significance" ~ 0.5,
                                  acmg_classification == "Benign" ~ 0.2,
                                  acmg_classification == "Likely_benign" ~ 0.2, TRUE ~ 0.1),
           Pathogenicity.score = ClinVar.score * ACMG.score,
           MAX_AF = case_when(is.na(MAX_AF) ~ min_af_nonzero, 
                              MAX_AF == "" ~  min_af_nonzero,
                              MAX_AF == 0 ~ min_af_nonzero,
                              TRUE ~ MAX_AF)) %>% 
    mutate(MAX_AF = as.numeric(MAX_AF),
           AF.score = -log10(MAX_AF),
           AF.score = normalize_minmax(AF.score))

  return(final.data)
}

# prediction scores
predict.score <- function(data){
  # Normalise scores to 0~1 
  normalize_alphamissense <- function(x) {
    out <- numeric(length(x))
    out[x <= 0.34] <- (x[x <= 0.34] / 0.34) * 0.33
    out[x > 0.34 & x <= 0.564] <- ((x[x > 0.34 & x <= 0.564] - 0.34) / (0.564 - 0.34)) * 0.33 + 0.33
    out[x > 0.564] <- ((x[x > 0.564] - 0.564) / (1 - 0.564)) * 0.34 + 0.66
    return(out)
  }
  
  normalize_spliceai <- function(x) {
    out <- numeric(length(x))
    out[x <= 0.1] <- (x[x <= 0.1] / 0.1) * 0.33
    out[x > 0.1 & x <= 0.5] <- ((x[x > 0.1 & x <= 0.5] - 0.1) / 0.4) * 0.33 + 0.33
    out[x > 0.5] <- ((x[x > 0.5] - 0.5) / 0.5) * 0.34 + 0.66
    return(out)
  }
  
  normalize_revel <- function(x){
    out <- numeric(length(x))
    out[x <= 0.29] <- (x[x <= 0.29] / 0.29) * 0.33
    out[x > 0.29 & x <= 0.643] <- ((x[x > 0.29 & x <= 0.643] - 0.29) / (0.643-0.29)) * 0.33 + 0.33
    out[x > 0.644] <- ((x[x > 0.644] - 0.644) / 0.644) * 0.34 + 0.66
    return(out)
  }
  
  clean_revel <- function(x) {
    if (is.na(x) | x == "." | x == "") {
      return(NA_real_)
    } else {
      splits <- unlist(strsplit(x, split = "&"))
      splits <- splits[splits != "." & splits != ""]
      nums <- as.numeric(splits)
      nums <- unique(nums)
      nums <- nums[!is.na(nums)]
      if (length(nums) == 0) {
        return(NA_real_)
      } else {
        return(mean(nums))
      }
    }
  }
  
  final.data <- data %>%
    mutate(spliceAI_max = pmax(SpliceAI_pred_DS_AG, SpliceAI_pred_DS_AL,SpliceAI_pred_DS_DG, SpliceAI_pred_DS_DL, na.rm = TRUE),
           spliceAI_max = if_else(is.na(spliceAI_max), 0 , spliceAI_max),
           AlphaMissense_pathogenicity = if_else(is.na(AlphaMissense_pathogenicity), 0 , AlphaMissense_pathogenicity),
           REVEL_score = sapply(REVEL_score,clean_revel),
           REVEL_score = as.numeric(if_else(is.na(REVEL_score), 0, REVEL_score))) %>% 
    mutate(Prediction.score = case_when(
      (IMPACT == "HIGH" & !grepl("^splice",Consequence) & spliceAI_max == 0) ~ 0.99,
      (grepl("^splice",Consequence) & spliceAI_max != 0) ~ normalize_spliceai(spliceAI_max),
      # since missense variant are moderate impact, the weight is only 0.66
      grepl("^missense_variant", Consequence) ~
       ((normalize_alphamissense(AlphaMissense_pathogenicity) + normalize_revel(REVEL_score)) / 2 ) * 0.66,
      TRUE ~ 0)
    )

  return(final.data)
}


# ClinVar phenotypes equal to HPO terms
check_phenotype_match <- function(clinvar_string, hpo_term) {
  clinvar_string <- as.character(clinvar_string)
  clinvar_string <- gsub("&_","_", clinvar_string)
  hpo_term <- as.character(hpo_term)
  
  if (is.na(clinvar_string) | is.na(hpo_term)) {
    return(NA)
  }
  
  hpo_terms <- toupper(unlist(strsplit(split = ",", hpo_term)))
  clinvar_terms <- unlist(strsplit(clinvar_string, split = "&"))
  
  for (hpo in hpo_terms) {
    hpo_words <- toupper(unlist(strsplit(trimws(hpo), " ")))
    
    for (cv in clinvar_terms) {
      cv_clean <- gsub("_", " ", cv) 
      cv_clean <- gsub("-"," ", cv_clean)
      cv_words <- toupper(unlist(strsplit(trimws(cv_clean), " ")))
      
      if (all(hpo_words %in% cv_words) || all(cv_words %in% hpo_words)) {
        return("matched")
      }
    }
  }
  
  return(NA)
}


# Pheno2Gene
pheno2gene.score <- function(hpo){
  hpo <- unlist(strsplit(hpo,split = ","))
  
  cols <- c('Gene','Pheno2Gene.Rank','Pheno2Gene.Score')
  pheno2gene_df = as_tibble(setNames(as.list(rep(NA, length(cols))), cols))
  hpo_tern_all <- NULL
  
  for(i in 1:length(hpo)){
    print(paste0("Now processing Pheno2Gene for ", hpo[i]))
    hpo.id <- hpo[i]
    pheno2gene.url <- paste0("https://phen2gene.wglab.org/api?HPO_list=", hpo.id)
    hpo.url <- paste0("https://ontology.jax.org/api/hp/terms/HP%3A",unlist(strsplit(split = ":", hpo.id))[2])
    
    # API
    response1 <- GET(pheno2gene.url)
    response2 <- GET(hpo.url)
    
    if (status_code(response1) == 200 & status_code(response2) == 200) {
      parsed1 <- fromJSON(content(response1, as = "text", encoding = "UTF-8"))
      parsed2 <- fromJSON(content(response2, as = "text", encoding = "UTF-8"))
      
      
      if (!is_empty(parsed1$results)){
        pheno2gene_df_tmp <- parsed1$results %>%
          select(Gene, Rank, Score) %>% 
          rename(Pheno2Gene.Rank = Rank, Pheno2Gene.Score = Score) %>%
          mutate(Pheno2Gene.Score = as.numeric(Pheno2Gene.Score))
        pheno2gene_df <- rbind(pheno2gene_df, pheno2gene_df_tmp)
        if (grepl("cancer", parsed2$name, ignore.case = T)){
          add.neoplasm <- gsub("cancer","neoplasm", parsed2$name)
        }else if (grepl("carcinoma", parsed2$name, ignore.case = T)){
          add.neoplasm <- gsub("carcinoma","neoplasm", parsed2$name)
        }else{
          add.neoplasm <- NULL
        }
        
        hpo_tern_all <- c(hpo_tern_all,as.character(parsed2$name), as.character(parsed2$synonyms), add.neoplasm)
        print(hpo_tern_all)
      }
    }
  }
  pheno2gene_df <- pheno2gene_df %>% group_by(Gene) %>% mutate(Pheno2Gene.Score = mean(Pheno2Gene.Score)) %>% ungroup()
  return(list(results = pheno2gene_df, hpo_name = hpo_tern_all))
  
}

# MR scores (EpiGraphDB, Gene - SNP - outcomes)
# if Gene and SNP and outcomes are same, then need to add the scores, and if the scores are smaller than 0.05 
mr.data <- function(hpo){
  hpo <- unlist(strsplit(hpo, split = ","))
  mr <- data.frame()
  for (j in hpo){
    #eqtl
    eqtl.col <- c('gene.ensembl_id','gene.name','gwas.id','gwas.trait','r.beta','r.se','r.p','r.rsid','Type')
    eqtl.mr <- as_tibble(setNames(as.list(rep(NA, length(eqtl.col))), eqtl.col))
    eqlt.mr.tpm <- xqtl_single_snp_mr(outcome_trait = j, qtl_type = "eQTL") %>% mutate(Type = "eQTL") 
    eqtl.mr <- rbind(eqtl.mr, eqlt.mr.tpm)
    mr <- rbind(mr, eqtl.mr)
    
    #pqtl
    pqtl.col <- c('gene.ensembl_id','gene.name','gwas.id','gwas.trait','r.beta','r.se','r.p','r.rsid','Type')
    pqtl.mr <- as_tibble(setNames(as.list(rep(NA, length(eqtl.col))), eqtl.col))
    pqtl.mr.tpm <- xqtl_single_snp_mr(outcome_trait = j, qtl_type = "pQTL") %>% mutate(Type = "pQTL")
    pqtl.mr <- rbind(pqtl.mr, pqtl.mr.tpm)
    
    mr <- rbind(mr, pqtl.mr)
  }
  return(mr)
}

## ---- Adding scores (with HPO terms) ----

cat("\nStart to calculate the prioritizing scores.\n(1) Pathogenicity and prediction scores.\n")

# Add scores for basic information
result <- patho.score(result)
result <- predict.score(result)

# Add scores for Phenotypes
if(hpo_terms!= ""){
  cat("\n(2) xQTL scores.\n")
  
  # Add scores for xQTL
  result <- result %>% 
    left_join(., eqtl.score, by = c("variant_info","SYMBOL" = "eqtl.gene")) %>% unique() %>% 
    left_join(., pqtl.score %>% filter(eqtl.gene == pqtl.gene) %>% select(-eqtl.gene), by = c("variant_info", "SYMBOL" = "pqtl.gene")) %>% unique() %>% 
    mutate(final.eqtl.normalised.score = if_else(is.na(final.eqtl.normalised.score), 0, final.eqtl.normalised.score),
           final.pqtl.normalised.score = if_else(is.na(final.pqtl.normalised.score), 0, final.pqtl.normalised.score))
  
  cat("\n(3) Phenotype scores.\n")
  
  pheno2gene.df <- pheno2gene.score(hpo_terms)
  pheno2gene.data <- pheno2gene.df$results %>% na.omit(Gene)
  pheno2gene.hponame <- paste(pheno2gene.df$hpo_name,collapse = ",")
  if(length(grep("Breast cancer", pheno2gene.df$hpo_name))> 0 ){
    pheno2gene.hponame <- paste(pheno2gene.hponame, "Malignant tumor of breast" ,sep = ",")
  }
  pheno2gene.data$Gene = as.character(pheno2gene.data$Gene)
  exomiser <- read.csv(exomiser_file,header = T,sep = "\t") %>% 
    mutate(variant_info = paste(CONTIG, START, REF, ALT,sep = "_"),
           variant_info = paste0("chr", variant_info))
  
  if(as.character(inheritance) != "" ){
    print(paste0("Inheritance matched with ", inheritance ," will be used for scoring."))
    exomiser = exomiser %>% filter(MOI == inheritance) %>% select(variant_info,GENE_SYMBOL, EXOMISER_GENE_PHENO_SCORE, EXOMISER_VARIANT_SCORE) %>% unique()
  }else{
    print("No inheritance for filtering exomiser. The maximum gene-phenotype score will be used. ")
    exomiser = exomiser %>% select(variant_info, GENE_SYMBOL, EXOMISER_GENE_PHENO_SCORE, EXOMISER_VARIANT_SCORE) %>% unique() %>% 
      group_by(variant_info,GENE_SYMBOL) %>%
      slice_max(EXOMISER_GENE_PHENO_SCORE, n = 1, with_ties = FALSE) %>%
      ungroup()
  }
  
  result <- result %>% 
    left_join(., pheno2gene.data, by = c("SYMBOL" = "Gene")) %>% unique() %>% mutate(HPO_term = pheno2gene.hponame) %>% 
    left_join(., exomiser, by = c("SYMBOL" = "GENE_SYMBOL", "variant_info")) %>% unique() %>% 
    mutate(EXOMISER_GENE_PHENO_SCORE = if_else(is.na(EXOMISER_GENE_PHENO_SCORE), 0, EXOMISER_GENE_PHENO_SCORE),
           EXOMISER_VARIANT_SCORE = if_else(is.na(EXOMISER_VARIANT_SCORE), 0, EXOMISER_VARIANT_SCORE),
           Pheno2Gene.Score = if_else(is.na(Pheno2Gene.Score), 0 , Pheno2Gene.Score)) %>% 
    rowwise() %>% mutate(clinvar.match = check_phenotype_match(ClinVar_CLNDN, HPO_term)) %>% ungroup() %>% 
    mutate(Phenotype.score = ((Pheno2Gene.Score + EXOMISER_GENE_PHENO_SCORE)/2)) 
    #mutate(Phenotype.score = ((Pheno2Gene.Score + EXOMISER_GENE_PHENO_SCORE)/2) + Clinvar.pheno.score) 
      #mutate(Phenotype.score = ((Pheno2Gene.Score + EXOMISER_GENE_PHENO_SCORE)/2) + Clinvar.pheno.score)
  
  # Add MR phenotype scores
  mr <- mr.data(hpo_terms)
  mr <- mr %>% filter(r.p < 0.05) %>% na.omit(Type) %>% unique()
  
  if(nrow(mr) > 0 ){
    print("Adding MR results")
    result <- result %>% 
      mutate(MR.score = 0) %>% 
      left_join(., mr, by = c("rsID"="r.rsid", "SYMBOL"="gene.name")) %>% unique() %>% 
      mutate(MR.score = -log10(`r.p`) * `r.beta`,
             MR.score = normalize_minmax(MR.score))
  }else{
    print("No MR results to be added. ")
    result <- result %>% 
      mutate(MR.score = 0)
  }
  
  # GWAS score
  gwas.df <- read.csv(gwas_file,header = T,sep = "\t",quote = "") 
  gwas_data <- gwas.df %>% filter(grepl(paste(pheno2gene.df$hpo_name, collapse = "|"), ignore.case = T,fixed = F,DISEASE.TRAIT)) %>% 
    filter(P.VALUE < 0.05)
  if (nrow(gwas_data) > 0){
    gwas_data = gwas_data %>% mutate(GWAS.score = -log10(`P.VALUE`) * abs(OR.or.BETA)) %>% 
      select(SNPS, GWAS.score)
    result <- result %>% left_join(., gwas_data,by = c("rsID" = "SNPS")) %>% unique() %>% 
      mutate(GWAS.score = if_else(is.na(GWAS.score), 0, normalize_minmax(GWAS.score)))
  }else{
    result <- result %>% mutate(GWAS.score = 0)
  }
  
  final.result <- result %>% 
    mutate(Priortise.score = Pathogenicity.score + AF.score + Prediction.score + final.eqtl.normalised.score + final.pqtl.normalised.score + Phenotype.score + EXOMISER_VARIANT_SCORE + GWAS.score) %>% 
    arrange(clinvar.match, desc(Priortise.score)) %>% unique() 
  loc <- final.result %>% filter(grepl("^LOC", SYMBOL))
  final.result <- final.result %>% filter(!grepl("^LOC", SYMBOL)) %>% rbind(., loc) 
  
  # rearrange by inheritance
  if(inheritance != "" & inheritance %in% c("XR","AR")){
    final.result <- final.result %>% 
      mutate(homo_flag = if_else(Zygosity == "Homozygous", 1, 0)) %>%
      arrange(clinvar.match, desc(homo_flag), desc(Priortise.score)) %>% select(-homo_flag)
  }
  
  final.result <- final.result %>% mutate(Rank = rownames(.))
  
}else{
  cat("\n(2) xQTL scores.\n")
  
  # Add scores for xQTL
  result <- result %>% 
    left_join(., eqtl.score, by = c("variant_info")) %>% unique() %>% 
    left_join(., pqtl.score, by = c("variant_info")) %>% unique() %>% 
    mutate(final.pqtl.normalised.score = if_else(is.na(final.pqtl.normalised.score), 0, final.pqtl.normalised.score),
           final.eqtl.normalised.score = if_else(is.na(final.eqtl.normalised.score), 0, final.eqtl.normalised.score))
  
  cat("\n(3) Phenotype scores.\n")
  
  final.result <- result %>% 
    mutate(Priortise.score = Pathogenicity.score + AF.score + Prediction.score + final.eqtl.normalised.score + final.pqtl.normalised.score) %>% 
    arrange(desc(Priortise.score)) %>% unique() %>% 
    mutate(Rank = rownames(.))
}

# df <- final.result %>% select(Priortise.score,SYMBOL, Rank,variant_info,IMPACT, Consequence, ClinVar_CLNSIG, ClinVar_CLNSIGCONF,acmg_classification,
#                               Prediction.score,final.eqtl.normalised.score, final.pqtl.normalised.score, Pheno2Gene.Score,EXOMISER_GENE_PHENO_SCORE,
#                               EXOMISER_VARIANT_SCORE, Pathogenicity.score,clinvar.match)

write.table(final.result, paste0(user_dir, sampleid, ".priorisation.txt"),quote = F,sep = "\t",row.names = F)


## ---- Network generation ----

# 1. variant - gene (from raw data)
# 2. variant - drug (PharmGKB)
# 3. variant - gene (eQTL, qval)
# 4. variant - gene (pQTL, pval)
# 5. gene -drug (DGIdb)
# 6. gene - gene (String)
# 7. variant - phenotype (raw data, GWAS, MR)

# network.df (source, target, interaction type, score)

cat("\nGenerating Network Files.")

# 1. variant - gene (from raw data)
vg.raw <- final.result %>% select(variant_info, rsID, SYMBOL) %>% unique() %>% 
  mutate(Type = "Variant-Gene-Original", Weight = 0,
         from.type = "Variant", to.type = "Gene",from = variant_info) %>% 
  rename(to = SYMBOL) %>% select(from, from.type, to, to.type, Type, Weight)

# 2. variant - drug (PharmGKB)
pgx.df <- read.csv(pgx_file,header = T,sep = "\t",quote = "") 
if(nrow(pgx.df) > 0){
  vd.pharm <- pgx.df %>% select(Variant.Haplotypes, pgx.rsID, Drug.s.) %>% unique() %>% 
    mutate(Type = "Variant-Drug", Weight = 0,
           from.type = "Variant", to.type = "Drug",from = Variant.Haplotypes) %>% 
    rename(to = Drug.s.) %>% select(from, from.type, to, to.type, Type, Weight) %>% 
    filter(!is.na(to))
}else{
  vd.pharm = NULL
}


# 3. variant - gene (eQTL, qval)
eqtl.catalog.df <- read.csv(eqtl_catalog_file,header = T,sep = "\t")
vg.eqtl.catalog <- eqtl.catalog.df %>% select(eQTL_Catalog_variant,eQTL_Catalog_Gene.name, eQTL_Catalog_p_beta) %>% unique() %>% 
  mutate(Type = "eQTL", from.type = "Variant", to.type = "Gene",from = eQTL_Catalog_variant) %>% 
  rename(to = eQTL_Catalog_Gene.name,Weight = eQTL_Catalog_p_beta,) %>% select(from, from.type, to, to.type, Type, Weight) %>% 
  filter(!is.na(to) & to != "") %>% 
  mutate(Weight = paste0("eQTL Catalog p-beta:", Weight))
eqtl.gtex.df <- read.csv(eqtl_gtex_file,header = T,sep = "\t")
vg.eqtl.gtex <- eqtl.gtex.df %>% select(eQTL_GTEx_variant_id,eQTL_GTEx_gene_name, eQTL_GTEx_qval) %>% unique() %>% 
  mutate(Type = "eQTL", from.type = "Variant", to.type = "Gene",from = eQTL_GTEx_variant_id) %>% 
  rename(to = eQTL_GTEx_gene_name,Weight = eQTL_GTEx_qval,) %>% select(from, from.type, to, to.type, Type, Weight) %>% 
  filter(!is.na(to) & to != "") %>% 
  mutate(Weight = paste0("GTExv10 eQTL qvalue:", Weight))
vg.eqtl <- rbind(vg.eqtl.catalog, vg.eqtl.gtex)

# 4. variant - gene (pQTL, pval)
otg.df <- read.csv(otg_file, header = T,sep = "\t")
vg.pqtl <- otg.df %>% select(variant_info, OTG.pqtl.gene, OTG.pqtl_pval) %>% unique() %>% 
  mutate(Type = "pQTL", from.type = "Variant", to.type = "Protein",
         from = variant_info,
         OTG.pqtl_pval = paste0("Open Target Genetics pQTL pvalue:",OTG.pqtl_pval)) %>% 
  rename(to = OTG.pqtl.gene, Weight = OTG.pqtl_pval) %>% 
  select(from, from.type, to, to.type, Type, Weight) %>% 
  filter(!is.na(to))

# 5. gene - drug (DGIdb)
dgi <- read.csv(dgi_file,header = T,sep = "\t") %>% select(gene_claim_name, drug_claim_name,interaction_score)
gd.dgi <- dgi %>% rename(from = gene_claim_name, to = drug_claim_name, Weight = interaction_score) %>% 
  mutate(Type = "Gene-Drug",from.type = "Gene", to.type = "Drug") %>%
  select(from, from.type, to, to.type, Type, Weight) %>% unique() %>% 
  right_join(., result %>% select(SYMBOL), by = c("from" = "SYMBOL")) %>% unique() %>% 
  na.omit(to) %>% 
  mutate(Weight = if_else(is.na(Weight) |Weight == "NULL", "0", Weight)) %>% 
  mutate(Weight = as.numeric(Weight),
         Weight = paste0("Gene-drug interaction score:", Weight))

# 6. gene - gene (String)
involved.genes <- c(vg.raw %>% pull(to),vg.eqtl %>% pull(to), vg.pqtl %>% pull(to), gd.dgi %>% pull(from)) %>% unique()
string <- fread(string_file, header = T,sep = "\t")
gg.string <- string %>% select(protein1, protein2, combined_score) %>% 
  mutate(Type = "PPI",from.type = "Protein", to.type = "Protein",
         combined_score = combined_score/1000 ) %>% 
  filter(protein1 %in% involved.genes & protein2 %in% involved.genes) %>% 
  rename(from = protein1, to = protein2, Weight = combined_score) %>% select(from, from.type, to, to.type, Type, Weight) %>% unique() %>% 
  mutate(Weight = paste0("Protein-protein interaction score:", Weight))

# 7. variant - phenotype (raw data, GWAS, MR)
gwas.df <- read.csv(gwas_file,header = T,sep = "\t",quote = "") 
vp.raw <- final.result %>% select(variant_info, rsID, ClinVar_CLNDN) %>% unique() %>% 
  mutate(Type = "Variant-Phenotype-Original", Weight = 0,
         from.type = "Variant", to.type = "Phenotype",
         from = variant_info) %>% 
  rename(to = ClinVar_CLNDN) %>% select(from, from.type, to, to.type, Type, Weight) 

vp.gwas <- gwas.df %>% select(SNPS, DISEASE.TRAIT,P.VALUE) %>% unique() %>% 
  mutate(Type = "Variant-Phenotype-GWAS", Weight = P.VALUE,
         from.type = "Variant", to.type = "Phenotype",
         from = SNPS) %>% 
  rename(to = DISEASE.TRAIT) %>% select(from, from.type, to, to.type, Type, Weight) %>% 
  mutate(Weight = paste0("GWAS Catalog pvalue:", Weight))


if(nrow(pqtl.df) > 0){
  vp.mr1 <- pqtl.df %>% select(`EpiGraphDB-rsID`,`EpiGraphDB-expID`,`EpiGraphDB-pvalue`) %>% unique() %>% 
    mutate(Type = "Variant-Phenotype-GWAS", Weight = `EpiGraphDB-pvalue`,
           from.type = "Variant", to.type = "Gene",
           from = `EpiGraphDB-rsID`) %>% 
    rename(to = `EpiGraphDB-expID`) %>% select(from, from.type, to, to.type, Type, Weight) %>% 
    mutate(Weight = paste0("EpiGraphDB mendelian randomization qvalue:", Weight))
  
  vp.mr2 <- pqtl.df %>% select(`EpiGraphDB-rsID`,`EpiGraphDB-outID`,`EpiGraphDB-pvalue`) %>% unique() %>% 
    mutate(Type = "Variant-Phenotype-GWAS", Weight = `EpiGraphDB-pvalue`,
           from.type = "Variant", to.type = "Phenotype",
           from = `EpiGraphDB-rsID`) %>% 
    rename(to = `EpiGraphDB-outID`) %>% select(from, from.type, to, to.type, Type, Weight) %>% 
    mutate(Weight = paste0("EpiGraphDB mendelian randomization qvalue:", Weight))
  vp.mr <- rbind(vp.mr1, vp.mr2) %>% filter(!is.na(to))
}else{
  vp.mr <- NULL
}

network.all <- rbind(vg.raw, vd.pharm, vg.eqtl, vg.pqtl, gd.dgi,gg.string,vp.mr) %>% na.omit(Type)

write.table(network.all,paste0(user_dir,sampleid, ".network.txt"),quote = F,sep = "\t",row.names = F)
