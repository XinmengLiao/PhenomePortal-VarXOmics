library(dplyr)
library(tidyverse)
library(data.table)
library(httr)
library(jsonlite)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(GO.db)
library(KEGGREST)
for (i in c("select","filter", "mutate","rename", "left_join", "slice")){
  conflicted::conflict_prefer(i, "dplyr")
}
rm(i)
conflicted::conflicts_prefer(stats::sd)
conflicted::conflicts_prefer(httr::content)
conflicted::conflicts_prefer(plotly::layout)


## ---- Load files ----
#user input
args <- commandArgs(trailingOnly = TRUE)
sampleid <- args[1]
output_dir <- args[2]
genes_file <- args[3]
genes <- read.csv(genes_file,header = F,sep = "\t")
genes <- unlist(genes$V1) %>% unique()
print("Now doing GO and KEGG for the user selected genes,")

load("KEGG_hsa_20251021.Rdata")
load("GO_hsa_20240624.Rdata")

## ---- Functions ----

# GO 
get_GO_annotation <- function(genes) {
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(GO.db)
  library(dplyr)
  
  gene_df <- bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  
  # Retrieve the terms 
  go_annot <- AnnotationDbi::select(org.Hs.eg.db, keys=gene_df$ENTREZID, 
                                    columns=c("GO", "ONTOLOGY"), keytype="ENTREZID")
  go_annot <- left_join(go_annot, gene_df, by="ENTREZID")
  
  go_ids <- unique(go_annot$GO)
  index <- which(is.na(go_ids))
  if(length(index) > 0){
    go_ids <- go_ids[-index]
  }
  go_term_map <- data.frame(GO = go_ids, TERM = Term(GOTERM[go_ids]))
  go_annot_full <- left_join(go_annot, go_term_map, by="GO") %>%
    select(SYMBOL, GO, TERM, ONTOLOGY) %>% filter(!is.na(GO))
  
  go_count_list <- go_annot_full %>% group_by(TERM, ONTOLOGY) %>% summarise(count = n_distinct(SYMBOL), .groups = "drop") %>% 
    arrange(desc(count))
  
  # Enrichment Analysis
  ego_bp <- enricher(
    gene = gene_df$ENTREZID,
    TERM2GENE = BPterm2gene_df,
    TERM2NAME = BPterm2name_df,
    pAdjustMethod = "BH",
    pvalueCutoff = 1) %>% as.data.frame() %>% mutate(Type = "BP")
  
  ego_mf <- enricher(
    gene = gene_df$ENTREZID,
    TERM2GENE = MFterm2gene_df,
    TERM2NAME = MFterm2name_df,
    pAdjustMethod = "BH",
    pvalueCutoff = 1) %>% as.data.frame() %>% mutate(Type = "MF")
  
  ego_cc <- enricher(
    gene = gene_df$ENTREZID,
    TERM2GENE = CCterm2gene_df,
    TERM2NAME = CCterm2name_df,
    pAdjustMethod = "BH",
    pvalueCutoff = 1) %>% as.data.frame() %>% mutate(Type = "CC")
  
  enrich <- rbind(ego_bp, ego_mf, ego_cc)
  
  return(list(GO.result = go_annot_full, GO.count = go_count_list,enrich = enrich))
}


# KEGG
get_KEGG_annotation <- function(genes) {
  
  # SYMBOL to ENTREZID
  gene_df <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  
  kk <- enricher(gene = gene_df$ENTREZID,
                TERM2GENE = path2gene_df,
                TERM2NAME = path2name_df,pvalueCutoff = 1) %>% as.data.frame()
  
  return(kegg = kk)
}

## ---- Generating outputs ---- 

# GO
go <- get_GO_annotation(genes)

# KEGG
kegg <- get_KEGG_annotation(genes)

print(kegg$kegg)

# Output text
write.table(go$GO.result,  paste0(output_dir,"/user_selected_GO-all.txt"), quote = F, sep = '\t', row.names = F)
write.table(go$GO.count, paste0(output_dir,"/user_selected_GO-all-count.txt"), quote = F, sep = '\t', row.names = F)
write.table(go$enrich, paste0(output_dir,"/user_selected_GO-all-enrich.txt"), quote = F, sep = '\t', row.names = F)

write.table(kegg,paste0(output_dir,"/user_selected_KEGG-all-enrich.txt"), quote = F, sep = '\t', row.names = F)
