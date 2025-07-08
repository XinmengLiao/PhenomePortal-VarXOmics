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
library(xlsx)
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
#sampleid = 'P00110_11'
output_dir <- args[2]
#output_dir = "/Users/xinmengliao/Documents/Project/20250516_Webserver/Case_Study/BreastCancer/result"
genes_file <- args[3]
genes <- read.csv(genes_file,header = F,sep = "\t")
#genes <- read.csv("Case_Study/BreastCancer/result/filtered.genes.txt",header = F,sep = "\t")
genes <- unlist(genes$V1) %>% unique()
print("Now doing GO and KEGG for the user selected genes,")

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
  ego_bp <- enrichGO(gene = gene_df$ENTREZID,
                     OrgDb = org.Hs.eg.db,
                     keyType = "ENTREZID",
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 1) %>% as.data.frame() %>% mutate(Type = "BP")
  
  ego_mf <- enrichGO(gene = gene_df$ENTREZID,
                     OrgDb = org.Hs.eg.db,
                     keyType = "ENTREZID",
                     ont = "MF",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 1) %>% as.data.frame() %>% mutate(Type = "MF")
  
  ego_cc <- enrichGO(gene = gene_df$ENTREZID,
                     OrgDb = org.Hs.eg.db,
                     keyType = "ENTREZID",
                     ont = "CC",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 1) %>% as.data.frame() %>% mutate(Type = "CC")
  
  enrich <- rbind(ego_bp, ego_mf, ego_cc)
  
  return(list(GO.result = go_annot_full, GO.count = go_count_list,enrich = enrich))
}


# KEGG
get_KEGG_annotation <- function(genes) {
  
  # SYMBOL to ENTREZID
  gene_df <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  
  # Retrieve the terms 
  # Get KEGG Pathway IDs
  kegg_annot <- AnnotationDbi::select(org.Hs.eg.db, 
                                      keys = gene_df$ENTREZID, 
                                      columns = "PATH", 
                                      keytype = "ENTREZID") %>%
    left_join(gene_df, by = "ENTREZID") %>%
    filter(!is.na(PATH))  
  
  kegg_ids <- unique(kegg_annot$PATH)
  
  kegg_term_map <- lapply(kegg_ids, function(pid) {
    res <- tryCatch({
      info <- keggGet(paste0("hsa", pid))[[1]]
      data.frame(PATH = pid, PATHWAY_NAME = info$NAME, stringsAsFactors = FALSE)
    }, error = function(e) {
      data.frame(PATH = pid, PATHWAY_NAME = NA, stringsAsFactors = FALSE)
    })
    return(res)
  }) %>% bind_rows() %>% filter(!is.na(PATHWAY_NAME))
  
  # Merge
  kegg_annot_full <- left_join(kegg_annot, kegg_term_map, by = "PATH") %>%
    select(SYMBOL, PATH, PATHWAY_NAME) %>%
    distinct() %>% filter(!is.na(PATH))
  
  kegg_count_list <- kegg_annot_full %>% group_by(PATHWAY_NAME) %>% summarise(count = n_distinct(SYMBOL), .groups = "drop") %>% 
    arrange(desc(count ))
  
  
  # Enrichment Analysis
  kk <- enrichKEGG(gene = gene_df$ENTREZID,
                   organism = "hsa",
                   keyType = "kegg",
                   pvalueCutoff = 1) %>% as.data.frame()
  
  return(list(kegg.result = kegg_annot_full, kegg.count = kegg_count_list, kegg = kk))
}

## ---- Generating outputs ---- 

# GO
go <- get_GO_annotation(genes)

# KEGG
kegg <- get_KEGG_annotation(genes)


# Output text
write.xlsx(go$GO.result, paste0(output_dir,"/user_selected_GO-KEGG.xlsx"),sheetName = "GO-all",append = T)
write.xlsx(go$GO.count, paste0(output_dir,"/user_selected_GO-KEGG.xlsx"),sheetName = "GO-all-count",append = T)
write.xlsx(go$enrich, paste0(output_dir,"/user_selected_GO-KEGG.xlsx"),sheetName = "GO-all-enrich",append = T)

gc()

write.xlsx(kegg$kegg.result, paste0(output_dir,"/user_selected_GO-KEGG.xlsx"), sheetName = "KEGG-all",append = T)
write.xlsx(kegg$kegg.count, paste0(output_dir,"/user_selected_GO-KEGG.xlsx"), sheetName = "KEGG-all-count",append = T)
write.xlsx(kegg$kegg, paste0(output_dir,"/user_selected_GO-KEGG.xlsx"), sheetName = "KEGG-all-enrich",append = T)


# kegg bar plot
keggbarplot <- function(data){
  if(nrow(data)>0){
    bar.data <- data %>% mutate(logp = -log10(p.adjust)) %>% arrange(category,Description) 
    bar.data$Description <- factor(bar.data$Description, levels =  bar.data$Description)
    p.kegg.bar <- ggplot(bar.data, aes(y = Count, x = Description,fill = category)) +
      geom_bar(stat = 'identity') +
      scale_fill_scico_d(palette = "romaO", direction = 1, alpha = 0.9,begin = 0.4, end = 0.9) +
      geom_line(aes(y = logp, x =  Description, group = 1), color = "#ffe6a7",size = 0.7, show.legend = FALSE) +
      geom_point(aes(y = logp,x = Description), color = "#ffe6a7", size = 3, show.legend = FALSE) +
      theme_classic() + labs(y = "Gene count", fill = "", x = "")+ 
      theme(axis.title.x.top = element_text(), legend.position = 'top',
            axis.text.x = element_text(angle = 270,hjust = 0),
            axis.title.x.bottom = element_text(),
            axis.text.x.bottom = element_text(),
            axis.title.y = element_text(),
            panel.grid.major.y = element_blank())+
      scale_y_continuous(
        name = expression(-log[10](p.adjust)),
        sec.axis = sec_axis(~ ., name = "Number of genes", breaks = bar.data$Count, labels = bar.data$Count)
      )
    
    y_levels <- max(bar.data$logp) - min(bar.data$logp)
    x_levels <- length(unique(bar.data$Description))
    width <- y_levels * 6
    height <- x_levels
    max_pixels <- 50000
    dpi = 600
    if(width * dpi > max_pixels){
      width = 50000/300
      dpi = 300
    }else if (height * dpi > max_pixels){
      height = 50000/300
      dpi = 300
    }
    ggsave(p.kegg.bar, filename = paste0(output_dir, '/Figures/user_selected_KEGGEnrich.png'),
           width = width, height = height, dpi = dpi,limitsize = F)
  }
}

gobarplot <- function(data, type){
  if(nrow(data)>0){
    bar.data <- data %>% mutate(logp = -log10(p.adjust)) %>% arrange(logp)
    bar.data$Description <- factor(bar.data$Description, levels =  bar.data$Description)
    p.go.bar <- ggplot(bar.data, aes(x = logp, y = Description, color = logp, size = logp)) +
      geom_point(stat = 'identity') +
      theme_bw()+
      scale_color_gradient(low = "blue", high = "red")+
      labs(y = paste0("GO-",type), x = "-log10(p.adj)", size ="-log10(p.adj)", color = "-log10(p.adj)" )+ 
      theme(axis.title.x.top = element_text(), legend.position = 'top',
            axis.text.x = element_text(angle = 270,hjust = 0),
            axis.title.y = element_text())+
      scale_x_continuous(breaks = seq(round(min(bar.data$logp),1), round(max(bar.data$logp),1), 1))
    
    x_levels <- max(bar.data$logp) - min(bar.data$logp)
    y_levels <- length(unique(bar.data$Description))
    width <- x_levels * 4
    height <- y_levels / x_levels
    max_pixels <- 50000
    dpi = 600
    if(width * dpi > max_pixels){
      width = 50000/300
      dpi = 300
    }else if (height * dpi > max_pixels){
      height = 50000/300
      dpi = 300
    }
    ggsave(p.go.bar, filename = paste0(output_dir, '/Figures/user_selected_', type, 'Enrich.png'),
           width = width, height = height, dpi = dpi,limitsize = F)
  }
}

# Save bar plots
keggbarplot(kegg$kegg)
gobarplot(go$enrich %>% filter(Type == "BP"), type = "BP")
gobarplot(go$enrich %>% filter(Type == "MF"), type = "MF")
gobarplot(go$enrich %>% filter(Type == "CC"), type = "CC")


## ---- circus plot ----
circusplot <- function(data, type){
  if(nrow(data) > 0 ){
    circusdata <- data %>% separate_rows(geneID, sep = "/") %>% 
      left_join(gene_df, by = c("geneID" = "ENTREZID")) %>% 
      rename(Genes = SYMBOL, Pathway = Description) %>% select(Pathway, Genes) %>% unique() %>% 
      mutate(Pathway = str_wrap(Pathway, width = 15))
    pdf(paste0(output_dir, '/Figures/user_selected_', type,'EnrichCircus.pdf'), width = 20, height = 20)
    circos.clear()
    circos.par(track.margin = c(0.05, 0.05), 
               cell.padding = c(0.05, 0.05, 0.05, 0.05),
               start.degree = 90)
    
    chordDiagram(circusdata,annotationTrack = "grid")
    circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
      sector.name = get.cell.meta.data("sector.index")
      xlim = get.cell.meta.data("xlim")
      ylim = get.cell.meta.data("ylim")
      circos.text(mean(xlim), ylim[2] + 1, sector.name,  # 这里的 ylim[2] + 1 决定往外推
                  facing = "clockwise", 
                  niceFacing = TRUE, rot = 45,
                  adj = c(0, 0.5), 
                  cex = 1.5)
    }, bg.border = NA)
    dev.off()
  }
}

# Save circus plots
safe_circusplot <- function(data, type) {
  tryCatch({
    circusplot(data, type = type)
  }, error = function(e) {
    cat("Circus figure generate for", type, "failed, please check if data are too much\n")
  })
}

safe_circusplot(kegg$kegg, "KEGG")
safe_circusplot(go$enrich %>% filter(Type == "BP"), "BP")
safe_circusplot(go$enrich %>% filter(Type == "MF"), "MF")
safe_circusplot(go$enrich %>% filter(Type == "CC"), "CC")
