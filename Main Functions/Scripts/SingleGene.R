library(dplyr)
library(tidyverse)
library(data.table)
library(epigraphdb)
library(httr)
library(jsonlite)
library(g3viz)
for (i in c("select","filter", "mutate","rename", "left_join", "slice")){
  conflicted::conflict_prefer(i, "dplyr")
}
rm(i)
conflicted::conflicts_prefer(stats::sd)
conflicted::conflicts_prefer(httr::content)
conflicted::conflicts_prefer(plotly::layout)

args <- commandArgs(trailingOnly = TRUE)
gene <- args[1]
output_dir <- args[2]
cat(paste0("Queried gene: ", gene))
cat(paste0("Output directory: ", output_dir))

## ---- load files ----
path = "sysmed"
cat(paste0("\nSetting path to ", path))
if(path == "local"){
  result_file <- "/Users/xinmengliao/Documents/Project/20250516_Webserver/usr_input/user_input_gene.txt"
  clinvar.data <- fread("/Users/xinmengliao/Documents/Project/20250516_Webserver/Datasets/ClinVar/clinvar_20250504.simplify.txt", header = TRUE, sep = "\t")
  clingen.data <- read.csv("/Users/xinmengliao/Documents/Project/20250516_Webserver/Datasets/Clingen/ClingenforGene.txt", header = TRUE, sep = "\t")
  gencc.data <- read.csv("/Users/xinmengliao/Documents/Project/20250516_Webserver/Datasets/GenCC/GenCCforGene.txt", header = TRUE, sep = "\t", colClasses = "character")
  dgi <- read.csv("/Users/xinmengliao/Documents/Project/20250516_Webserver/Datasets/DGIdb/interactions.tsv",header = T,sep = "\t")
  string <- fread("/Users/xinmengliao/Documents/Project/20250516_Webserver/Datasets/xQTL/String/String.all.score.txt",header = T,sep = "\t")
  clinvar.hgsvp <- read.csv("/Users/xinmengliao/Documents/Project/20250516_Webserver/Datasets/ClinVar/hgvsp-plp.txt",header = TRUE, sep = "\t")
  output_file <- paste0(output_dir, "/", gene,"query.txt")
  letter <- unlist(strsplit(split = "",gene))[1]
  eqtl.catalog <- fread(
    paste0("/Users/xinmengliao/Documents/Project/20250516_Webserver/Datasets/xQTL/eQTL_Catalogue/forSingleGene/genes_",letter,".txt"),
    header = T,sep = "\t") %>% 
    filter(Gene.name == gene)
  eqtl.gtex <- read.csv(
    paste0("/Users/xinmengliao/Documents/Project/20250516_Webserver/Datasets/GTEx/forSingleGene/genes_",letter,".txt"),
    header = T,sep = "\t") %>% 
    filter(gene_name == gene)
  sqtl.gtex <- read.csv("/Users/xinmengliao/Documents/Project/20250516_Webserver/Datasets/xQTL/GTEx_Analysis_v10_sQTL_updated/GTEx.all.tissue.known.gene.txt",
    header = T,sep = "\t")
  
}else if (path == "sysmed"){
  clinvar.data <- fread("/mnt/storage_pool/Genomics/VarXOmics/Databases/clinvar_20250504.simplify.txt", header = TRUE, sep = "\t")
  clingen.data <- read.csv("/mnt/storage_pool/Genomics/VarXOmics/Databases/ClingenforGene.txt", header = TRUE, sep = "\t")
  gencc.data <- read.csv("/mnt/storage_pool/Genomics/VarXOmics/Databases/GenCCforGene.txt", header = TRUE, sep = "\t", colClasses = "character")
  dgi <- read.csv("/mnt/storage_pool/Genomics/VarXOmics/Databases/DGIdb_interactions.tsv",header = T,sep = "\t")
  string <- fread("/mnt/storage_pool/Genomics/VarXOmics/Databases/String.all.score.txt",header = T,sep = "\t")
  clinvar.hgsvp <- read.csv("/mnt/storage_pool/Genomics/VarXOmics/Databases/clinvar-hgvsp-plp.txt",header = TRUE, sep = "\t")
  output_file <- paste0(output_dir, "/", gene,"query.txt")
  letter <- unlist(strsplit(split = "",gene))[1]
  eqtl.catalog <- fread(
    paste0("/mnt/storage_pool/Genomics/VarXOmics/Databases/eQTL_Catalogue_forSingleGene/genes_",letter,".txt"),
    header = T,sep = "\t") %>% 
    filter(Gene.name == gene)
  eqtl.gtex <- read.csv(
    paste0("/mnt/storage_pool/Genomics/VarXOmics/Databases/eQTL_GTEx_forSingleGene/genes_",letter,".txt"),
    header = T,sep = "\t") %>% 
    filter(gene_name == gene)
  sqtl.gtex <- read.csv("/mnt/storage_pool/Genomics/VarXOmics/Databases/sqtl_GTEx.all.tissue.known.gene.txt",
                        header = T,sep = "\t")
}


## ---- Basic information: combine ClinVar, ClinGen, GenCC ---- 
# ClinVar
cat("\n 1. Searching ClinVar information. \n")
gene.df <- clinvar.data %>% 
  filter(str_detect(GENEINFO, gene))

if(nrow(gene.df) == 0){
    gene.df = data.frame("Note" = "No results or data found. " )
}else{
  gene.df <- gene.df %>% 
    separate_rows(GENEINFO, sep = "\\|") %>%
    filter(str_detect(GENEINFO, gene)) %>%
    mutate(Variant = paste(CHROM, POS, REF, ALT, sep = "_"),
           Variant = paste0("chr", Variant)) %>% 
    mutate(GENEINFO = sapply(strsplit(split = ":",GENEINFO),`[`,1)) %>% unique() %>% 
    rename(`Diseases sources` = CLNDISDB, Diseases = CLNDN, HGVS = CLNHGVS, `Review Criteria` = CLNREVSTAT, Pathogenicity = CLNSIG, 
           `Variant Type` = CLNVC, Gene = GENEINFO, rsID = RS, `Review Star` = Review.Status, variant_info = Variant) %>% 
    mutate(Pathogenicity = gsub("_", " ", Pathogenicity)) %>% unique()
}

write.table(gene.df, paste0(output_dir, "/", gene,"-clinvar.txt"), quote = F,sep = "\t",row.names = F)


# ClinGen
cat("2. Searching ClinGen information. \n")
clingen.data <- clingen.data %>%
  filter(Gene == gene) %>% unique() %>% 
  select(Gene, Disease, Inheritance,MONDO, Classification,)
if(nrow(clingen.data) > 0){
  colnames(clingen.data) <- paste0("ClinGen.",colnames(clingen.data))
  write.table(clingen.data, paste0(output_dir, "/", gene,"-clingen.txt"), quote = F,sep = "\t",row.names = F)
}else{
  write.table(data.frame("Note" = "No results or data found. " ),
              paste0(output_dir, "/", gene,"-clingen.txt"), quote = F,sep = "\t",row.names = F)
}


# GenCC
cat("3. Searching GenCC information. \n")
gencc.data <- gencc.data %>% 
  filter(Gene == gene) %>% unique()
if(nrow(gencc.data) > 0){
  colnames(gencc.data) <- paste0("GenCC.",colnames(gencc.data))
  write.table(gencc.data, paste0(output_dir, "/", gene,"-gencc.txt"), quote = F,sep = "\t",row.names = F)
}else{
  write.table(data.frame("Note" = "No results or data found. " ),
              paste0(output_dir, "/", gene,"-gencc.txt"), quote = F,sep = "\t",row.names = F)
}

## ---- xQTL: eQTL and pQTL ----
#### ---- eQTL Catalog 
cat("4. Searching eQTL information from eQTL Catalog. \n")
colnames(eqtl.catalog) <- paste0("eQTL_Catalog.",colnames(eqtl.catalog))
if(nrow(eqtl.catalog) > 0 ){
  write.table(eqtl.catalog, paste0(output_dir, "/", gene,"-eqtl_catalog.txt"),  quote = F,sep = "\t",row.names = F)
}else{
  write.table(data.frame("Note" = "No results or data found. " ),
              paste0(output_dir, "/", gene,"-eqtl_catalog.txt"), quote = F,sep = "\t",row.names = F)
}

#### ---- GTEx
cat("5. Searching eQTL information from GTEx. \n")
colnames(eqtl.gtex) <- paste0("eQTL_GTEx.",colnames(eqtl.gtex))
if(nrow(eqtl.gtex) > 0 ){
  write.table(eqtl.gtex, paste0(output_dir, "/", gene,"-eqtl_gtex.txt"),  quote = F,sep = "\t",row.names = F)
}else{
  write.table(data.frame("Note" = "No results or data found. " ),
              paste0(output_dir, "/", gene,"-eqtl_gtex.txt"), quote = F,sep = "\t",row.names = F)
}

#### ---- sQTL 
cat("6. Searching sQTL information from GTEx. \n")
colnames(sqtl.gtex) <- paste0("sQTL_GTEx.",colnames(sqtl.gtex))
sqtl.data <- sqtl.gtex %>% filter(sQTL_GTEx.gene_name == gene) %>% unique()
if(nrow(sqtl.data) > 0){
  write.table(sqtl.data, paste0(output_dir, "/", gene,"-sqtl_gtex.txt"),  quote = F,sep = "\t",row.names = F)
}else{
  write.table(data.frame("Note" = "No results or data found. " ),
              paste0(output_dir, "/", gene,"-sqtl_gtex.txt"), quote = F,sep = "\t",row.names = F)
}

#### ---- pQTL
cat("7. Searching pQTL information from EpiGraphDB. \n")
library(epigraphdb)
pqtl_tpm <- pqtl(query = gene, searchflag = "proteins",rtype = "simple",pvalue = 1)

if(nrow(pqtl_tpm) > 0){
  pqtl_tpm <- pqtl_tpm %>% select(-outID_mrbase, -nsnp, -direction) 
  colnames(pqtl_tpm) <- c('pQTL.Protein','pQTL.Disease/Phenotype','pQTL.pvalue','pQTL.rsID','pQTL.steiger_pvalue',
                          'pQTL.coloc_prob','pQTL.beta','pQTL.method','pQTL.trans_cis',
                          'pQTL.q_pvalue','pQTL.ld_check','pQTL.se')
  write.table(pqtl_tpm, paste0(output_dir, "/", gene,"-pqtl_epigraphdb.txt"),  quote = F,sep = "\t",row.names = F)
}else{
  write.table(data.frame("Note" = "No results or data found. " ),
              paste0(output_dir, "/", gene,"-pqtl_epigraphdb.txt"), quote = F,sep = "\t",row.names = F)
}

#### ---- eQTL MR (EpigrphDB)
cat("8. Searching eQTL MR information from EpiGraphDB. \n")
eqtlmr <- xqtl_single_snp_mr(exposure_gene  = gene, qtl_type = "eQTL")
if(nrow(eqtlmr) > 0){
  write.table(eqtlmr, paste0(output_dir, "/", gene,"-eqtlmr.txt"),  quote = F,sep = "\t",row.names = F)
}else{
  write.table(data.frame("Note" = "No results or data found. " ),
              paste0(output_dir, "/", gene,"-eqtlmr.txt"), quote = F,sep = "\t",row.names = F)
}

## ---- DGIdb ----
cat("9. Searching Gene-Drug interactions from DGIdb. \n")

dgi <- dgi %>% filter(gene_name == gene) %>% unique()
colnames(dgi) <- paste0("DGI.",colnames(dgi))

if(nrow(dgi) > 0){
  write.table(dgi, paste0(output_dir, "/", gene,"-dgidb.txt"),  quote = F,sep = "\t",row.names = F)
}else{
  write.table(data.frame("Note" = "No results or data found. " ),
              paste0(output_dir, "/", gene,"-dgidb.txt"), quote = F,sep = "\t",row.names = F)
}

## ---- Network generation ----
# 1. gene - variant (eQTL, qval)
# 2. gene - variant (MR, qval)
# 4. gene -drug (DGIdb)
# 5. gene - gene (String)

# network.df (source, target, interaction type, score)

cat("10. Generating Network for single gene query. \n")

#  1. gene - varian{t (eQTL, qval)
if(nrow(eqtl.catalog) > 0){
  gv.eqtl.catalog <- eqtl.catalog %>% select(eQTL_Catalog.Gene.name,eQTL_Catalog.variant,eQTL_Catalog.p_beta) %>% unique() %>%
    filter(!is.na(eQTL_Catalog.Gene.name) & !is.na(eQTL_Catalog.variant)) %>%
    mutate(Type = "eQTL", Weight = eQTL_Catalog.p_beta,
           from = eQTL_Catalog.Gene.name) %>%
    rename(to = eQTL_Catalog.variant) %>% select(from, to, Type, Weight)
}else{
  gv.eqtl.catalog <- NULL
}

if(nrow(eqtl.gtex) > 0){
  gv.eqtl.gtex <- eqtl.gtex %>% select(eQTL_GTEx.gene_name,eQTL_GTEx.variant_id,eQTL_GTEx.qval) %>% unique() %>%
    filter(!is.na(eQTL_GTEx.gene_name) & !is.na(eQTL_GTEx.variant_id)) %>%
    mutate(Type = "eQTL", Weight = eQTL_GTEx.qval,
           from = eQTL_GTEx.gene_name) %>%
    rename(to = eQTL_GTEx.variant_id) %>% select(from, to, Type, Weight)
}else{
  gv.eqtl.gtex = NULL
}

gv.eqtl <- rbind(gv.eqtl.catalog,gv.eqtl.gtex )

# 2. gene - variant (pQTL, qval)
if(nrow(pqtl_tpm) > 0){
  gv.pqtl1 <- pqtl_tpm %>% select(pQTL.Protein,`pQTL.Disease/Phenotype`,`pQTL.q_pvalue`) %>% unique() %>%
    filter(!is.na(pQTL.Protein) & !is.na(`pQTL.Disease/Phenotype`)) %>%
    mutate(Type = "pQTL", Weight = `pQTL.q_pvalue`,
           from = pQTL.Protein) %>%
    rename(to = `pQTL.Disease/Phenotype`) %>% select(from, to, Type, Weight)
}else{
  gv.pqtl1 <- NULL
}

if(nrow(pqtl_tpm) > 0){
  gv.pqtl2 <- pqtl_tpm %>% select(pQTL.Protein,`pQTL.rsID`,`pQTL.q_pvalue`) %>% unique() %>%
    filter(!is.na(pQTL.Protein) & !is.na(pQTL.rsID)) %>%
    mutate(Type = "pQTL", Weight = `pQTL.q_pvalue`,
           from = pQTL.Protein) %>%
    rename(to = pQTL.rsID) %>% select(from, to, Type, Weight)
}else{
  gv.pqtl2 <- NULL
}

gv.pqtl <- rbind(gv.pqtl1, gv.pqtl2 )

# 3. gene -drug (DGIdb)
if(nrow(dgi) > 0){
  gd <- dgi %>% select(DGI.gene_name, DGI.drug_name) %>% unique() %>% 
    filter(!is.na(DGI.drug_name) & !is.na(DGI.gene_name)) %>% 
    mutate(Type = "Gene-Drug", from = DGI.gene_name,Weight = 0) %>% 
    rename(to = DGI.drug_name) %>% select(from, to, Type, Weight) 
}else{
  gd = NULL
}

# 4. gene - gene (String)
gg.string <- string %>% select(protein1, protein2, combined_score) %>% 
  filter(protein1 ==gene) %>% 
  filter(!grepl("^ESNG",protein2)) %>% 
  mutate(protein2 = gsub("_HUMAN","",protein2), protein2 = toupper(protein2),
         Type = "PPI",combined_score = combined_score/1000 ) %>% 
  rename(from = protein1, to = protein2, Weight = combined_score) %>% select(from, to, Type, Weight) %>% unique()

all <- rbind(gv.eqtl, gv.pqtl, gd,gg.string)
if(nrow(all)>0){
  write.table(all, paste0(output_dir, "/", gene,"-network.txt"), quote = F,sep = "\t",row.names = F)
}else{
  write.table(data.frame("Note" = "No results or data found. " ),
              paste0(output_dir, "/", gene,"-network.txt"), quote = F,sep = "\t",row.names = F)
}


## ---- Lollipop plots (ClinVar for single gene search) ----
cat("11. Plotting Lollipop plot for single gene \n")

clinvar.lol.data <- clinvar.hgsvp %>% filter(GENEINFO == gene)
if(nrow(clinvar.lol.data)>0){
  plot.options <- g3Lollipop.theme(theme.name = "simple",title.text = gene,
                                   y.axis.label = paste0(gene, " Mutations"),
                                   legend.title = "Consequence")
  
  clinvar.lol.chart <- g3Lollipop(
    mutation.dat = clinvar.lol.data,gene.symbol = gene,aa.pos.col = 'AA_Position',
    protein.change.col = "Protein_Change",
    gene.symbol.col = "GENEINFO",output.filename = paste0(output_dir,"/Lollipop-", gene),
    plot.options =  plot.options)
  htmlwidgets::saveWidget(clinvar.lol.chart, paste0(output_dir,"/ClinVarPLP-Lollipop-", gene,".html"))
  jsonlite::write_json(clinvar.lol.data, path = paste0(output_dir, "/ClinVarPLP-Lollipop-", gene, ".json"), pretty = TRUE)
}

