library(dplyr)
library(ggplot2)
library(viridis)
library(karyoploteR)
library(forcats)
library(g3viz)
library(circlize)
library(tidyr)
for (i in c("select","filter", "mutate","rename", "left_join", "slice")){
  conflicted::conflict_prefer(i, "dplyr")
}
rm(i)

args <- commandArgs(trailingOnly = TRUE)
sampleid <- args[1]
user_selected_file <- args[2]
output_dir <- args[3]
setwd(output_dir)
print(paste0("Now generating plots for variant summary section", user_selected_file))
df <- read.csv(user_selected_file,header = T,sep = "\t") %>% unique()

## ---- Variant class ----
class <- df %>% select(variant_info, VARIANT_CLASS) %>% unique() %>% 
  mutate(VARIANT_CLASS = case_when(VARIANT_CLASS == "insertion" ~ "INS", VARIANT_CLASS == "deletion" ~ "DEL", VARIANT_CLASS == "SNV" ~ "SNV", T ~ "Other")) %>% 
  group_by(VARIANT_CLASS) %>% summarise(count = n(),.groups = "drop") %>%
  arrange(desc(VARIANT_CLASS)) %>%
  mutate(
    fraction = count / sum(count),
    ymax = cumsum(fraction),
    ymin = c(0, head(ymax, n=-1)),
    labelPosition = (ymax + ymin) / 2,
    label = paste0(VARIANT_CLASS)
  )
my_colors <- c("SNV" = "#fce38a", "INS" = "#f38181", "DEL" = "#95e1d3", "Other" = "#9f9f9f")

if(nrow(class) > 0 ){
  p.class <- ggplot(class, aes(y = count, x = "", fill = VARIANT_CLASS)) +
    geom_bar(stat = 'identity') +
    #geom_text(x = 0.5, aes(y = labelPosition, label = label), size = 5) +
    scale_fill_manual(values = my_colors) +
    coord_polar(theta = "y") +
    theme_void() + 
    labs(fill = "Variant class")+
    theme(legend.position = "right")
  
  ggsave(p.class, filename = "Figures/user_selected_VariantClass.pdf",width = 6,height = 6,units = "in",dpi = 600)
  write.table(class, "Figures_data/user_selected_VariantClass.txt",quote = F,sep = "\t",row.names = F)
}

## ---- IMPACT ----

high <- c('transcript_ablation','splice_acceptor_variant','splice_donor_variant','stop_gained','frameshift_variant','stop_lost',
          'start_lost','transcript_amplification','feature_elongation','feature_truncation')
moderate <- c('inframe_insertion','inframe_deletion','missense_variant','protein_altering_variant')
low <-  c('splice_donor_5th_base_variant','splice_region_variant','splice_donor_region_variant','splice_polypyrimidine_tract_variant',
          'incomplete_terminal_codon_variant','start_retained_variant','stop_retained_varian','synonymous_variant')
modifier <- c('coding_sequence_variant','mature_miRNA_variant','5_prime_UTR_variant','3_prime_UTR_variant','non_coding_transcript_exon_variant',
              'intron_variant','NMD_transcript_varian','non_coding_transcript_variant','coding_transcript_variant','upstream_gene_variant',
              'downstream_gene_variant','TFBS_ablation','TFBS_amplification','TF_binding_site_variant','regulatory_region_ablation',
              'regulatory_region_amplification','regulatory_region_variant','intergenic_variant','sequence_variant')
conseq <- df %>%
  select(variant_info, Consequence) %>%
  unique() %>%
  mutate(Consequence = sapply(strsplit(split = "&",Consequence),`[`,1)) %>% 
  group_by(Consequence) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(desc(count)) %>% na.omit(Consequence) %>% 
  mutate(Impact = case_when(Consequence %in% high ~ "High", Consequence %in% moderate ~ "Moderate",
                            Consequence %in% low ~ "Low", Consequence %in% modifier ~ "Modifier", TRUE ~ NA))
my_colors <- c("High" = "#eb99a7", "Moderate" = "#f9dca7", "Low" = "#b9dbab", "Modifier" = "#b4d0e8")

if(nrow(conseq)> 0 ){
  conseq$Consequence <- factor(conseq$Consequence, levels = rev(c(high, moderate, low ,modifier)))
  
  p.conseq <- ggplot(conseq, aes(x = Consequence, y = count, fill = Impact)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = my_colors) +
    labs(x = "", y = "No. of Variants") +
    theme_bw() +
    theme(
      legend.position = "right",
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_line(color = "#f0f0f0"),
      panel.grid.minor.x = element_blank(),
      axis.text.y = element_text(size = 11),
      axis.text.x = element_text(size = 10,angle = 270, vjust = 0,hjust = 0),
      axis.title.x = element_text(size = 12, margin = margin(t = 10))
    ) +
    scale_y_continuous(breaks= seq(min(conseq$count), max(conseq$count), 20))
  
  width = max(conseq$count) - min(conseq$count)
  if(width > 100){width = width / 10}
  height = length(unique(conseq$Consequence))
  if(height > 10){height = height/2}
  
  ggsave(p.conseq, filename = "Figures/user_selected_VariantConsequence.pdf",width = width,height = height,units = "in",dpi = 600)
  write.table(conseq, "Figures_data/user_selected_VariantConsequence.txt",quote = F,sep = "\t",row.names = F)
}

## ---- ClinVar distribution ----
clinvar_data <- df %>% select(variant_info, ClinVar_CLNSIG) %>% 
  filter(ClinVar_CLNSIG != "") %>% 
  mutate(ClinVar_CLNSIG = sapply(strsplit(split = "&", ClinVar_CLNSIG),`[`,1)) %>% unique() %>% 
  group_by(ClinVar_CLNSIG) %>% summarise(count = n(),.groups = "drop")

if(nrow(clinvar_data)> 0 ){
  clinvar_data$ClinVar_CLNSIG <- factor(clinvar_data$ClinVar_CLNSIG, 
                                        levels = c("Pathogenic","Likely_pathogenic",
                                                   "Uncertain_significance","Conflicting_classifications_of_pathogenicity"))
  
  my_colors <- c("Pathogenic" = "#edafb8", "Likely_pathogenic" = "#f7e1d7", "Uncertain_significance" = "#ccd5ae", "Conflicting_classifications_of_pathogenicity" = "#e9edc9")
  
  p.clinvar <- ggplot(clinvar_data, aes(y = count , x = "", fill = ClinVar_CLNSIG)) + 
    geom_bar(stat = 'identity')+
    #geom_text(x = 0.5, aes(y = labelPosition, label = label), size = 5) +
    scale_fill_manual(values = my_colors) +
    coord_polar(theta = "y") +
    theme_void() + 
    labs(fill = "ClinVar Pathogenicity")+
    theme(legend.position = "right")
  
  ggsave(p.clinvar, filename = "Figures/user_selected_ClinVarPathogenicity.pdf",width = 6,height = 6,units = "in",dpi = 600)
  write.table(clinvar_data, "Figures_data/user_selected_ClinVarPathogenicity.txt",quote = F,sep = "\t",row.names = F)
}

## ---- SNP density ----
library(CMplot)
snp_data <- df %>% select(X.CHROM, POS) %>% mutate(SNP = 1:nrow(.)) %>% unique() %>%
  mutate(Chromosome = X.CHROM, Position = POS) %>% select(-X.CHROM, -POS) %>% 
  mutate(Chromosome = gsub("chr","", Chromosome),
         Chromosome = case_when(Chromosome == "X" ~ "23", 
                                Chromosome == "Y" ~ "24",
                                Chromosome == "M" ~ "25", T~ Chromosome),
         Chromosome = as.numeric(Chromosome))

if(nrow(snp_data)> 0 ){
  write.table(snp_data, "Figures_data/user_selected_SNP_Density.txt",quote = F,sep = "\t",row.names = F)
  original_wd <- getwd()
  setwd(paste0(output_dir, '/Figures'))
  CMplot(snp_data, 
         plot.type = "d",  # density plot
         bin.size = 1e6,   # 每 Mb 计算一个 bin
         chr.den.col = c("darkgreen", "yellow", "red"),  
         file = "jpg", cex = 5, lab.cex = 10,
         axis.cex = 1.2, 
         file.output = T, dpi = 600,
         verbose = TRUE, width = 16,height = 7,
         file.name ="user_selected_SNP_density")
  setwd(original_wd)
}

## ---- Gene Accumulation ----
gene.accu <- df %>% 
  select(variant_info, SYMBOL) %>% unique() %>% 
  group_by(SYMBOL) %>% 
  summarise(count = n(), .groups = "drop") %>% 
  arrange(desc(count),SYMBOL) 

gene.accu.fun <- function(gene.accu.set){
  gene.accu.set$SYMBOL <- factor(gene.accu.set$SYMBOL, levels = rev(gene.accu.set$SYMBOL))
  
  p.gene.accu <- ggplot(gene.accu.set, aes(x = count, y = SYMBOL, color = count,size = count)) +
    geom_point()+
    scale_size(range = c(3, 8)) + 
    #geom_bar(stat = "identity") +
    scale_color_viridis(option = "G", discrete = F,direction = 1, begin = 0.93, end = 0.7) +
    labs(x = "", y = "") +
    theme_bw()  +
    guides(size = "none")+
    theme(
      legend.position = "right",
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_line(color = "#f0f0f0"),
      panel.grid.minor.x = element_blank(),
      axis.text.y = element_text(size = 11),
      axis.text.x = element_text(size = 10, vjust = 0,hjust = 0),
      axis.title.x = element_text(size = 12, margin = margin(t = 10))
    )
  ggsave(p.gene.accu, filename = "Figures/user_selected_GeneAccumultate-AllorTop10.pdf",width = 6.18,height = 5.25,units = "in",dpi = 600)
  write.table(gene.accu, "Figures_data/user_selected_GeneAccumultate-all.txt",quote = F,sep = "\t",row.names = F)
}

if(nrow(gene.accu) > 10 ){
  gene.accu.set <- gene.accu %>% slice(1:10)
  gene.accu.fun(gene.accu.set)
}else if (nrow(gene.accu) <= 10 & nrow(gene.accu) > 0 ){
  gene.accu.fun(gene.accu)
}


## ---- Disease Accumulation ----
disease.accu <- df %>% 
  select(variant_info, ClinVar_CLNDN) %>% 
  unique() %>% 
  filter(!is.na(ClinVar_CLNDN) & 
           !ClinVar_CLNDN %in% c("", "not provided", "not specified", "not_provided", "not_specified")) %>%
  separate_rows(ClinVar_CLNDN, sep = "(?<!_)&(?!_)") %>%
  group_by(ClinVar_CLNDN) %>% 
  summarise(count = n(), .groups = "drop") %>% 
  arrange(desc(count), ClinVar_CLNDN) %>%
  filter(!is.na(ClinVar_CLNDN) & 
           !ClinVar_CLNDN %in% c("", "not provided", "not specified", "not_provided", "not_specified")) %>%
  mutate(ClinVar_CLNDN = gsub("&_", " ", ClinVar_CLNDN),
         ClinVar_CLNDN = gsub("_"," ", ClinVar_CLNDN))

disease.accu.fun <- function(disease.accu.set){
  disease.accu.set$ClinVar_CLNDN <- factor(disease.accu.set$ClinVar_CLNDN, levels = rev(disease.accu.set$ClinVar_CLNDN))
  
  p.disease.accu <- ggplot(disease.accu.set, aes(x = count, y = ClinVar_CLNDN, color = count, size = count)) +
    geom_point()+
    #geom_bar(stat = "identity") +
    scale_color_viridis(option = "G", discrete = F,direction = 1, begin = 0.6, end = 0.35) +
    labs(x = "", y = "") +
    theme_bw() + guides(size = "none")+ 
    scale_size_continuous(range = c(3, 8))+
    theme(
      legend.position = "right",
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_line(color = "#f0f0f0"),
      panel.grid.minor.x = element_blank(),
      axis.text.y = element_text(size = 11),
      axis.text.x = element_text(size = 10, vjust = 0,hjust = 0),
      axis.title.x = element_text(size = 12, margin = margin(t = 10))
    )+
    scale_x_continuous(breaks = seq(min(disease.accu.set$count), max(disease.accu.set$count),2))
  
  width = max(disease.accu.set$count) - min(disease.accu.set$count)
  if(width > 10){width = width / 2}
  height = length(unique(disease.accu.set$ClinVar_CLNDN))
  if(height > 5) { height = height /1.5}
  
  ggsave(p.disease.accu, filename = "Figures/user_selected_DiseaseAccumultate-AllorTop10.pdf",width = width,height = height,units = "in",dpi = 600)
  write.table(disease.accu, "Figures_data/user_selected_DiseaseAccumultate-all.txt",quote = F,sep = "\t",row.names = F)
}

if(nrow(disease.accu) > 10 ){
  disease.accu.set <- disease.accu %>% slice(1:10)
  disease.accu.fun(disease.accu.set)
}else if (nrow(disease.accu) > 0 & nrow(disease.accu) <= 10){
  disease.accu.fun(disease.accu)
}


## ---- Mutation locations ----
mutation_df <- df %>% select(X.CHROM, POS, VARIANT_CLASS,ClinVar_CLNSIG,SYMBOL, REF, ALT) %>% 
  mutate(VARIANT_CLASS = case_when(VARIANT_CLASS == "insertion" ~ "INS", VARIANT_CLASS == "deletion" ~ "DEL", 
                                   VARIANT_CLASS == "SNV" ~ "SNV", T ~ "Other"),
         END = case_when(VARIANT_CLASS == "INS" ~ POS + nchar(REF) - nchar(ALT),
                         VARIANT_CLASS == "DEL" ~ POS + nchar(REF) - nchar(ALT), T ~ POS)) %>% 
  mutate(ClinVar_CLNSIG = sapply(strsplit(ClinVar_CLNSIG, split = "&"), `[`,1)) %>% unique() %>% 
  mutate(color = case_when(ClinVar_CLNSIG == "Pathogenic"~"red", 
                           ClinVar_CLNSIG == "Likely_pathogenic" ~ "orange",
                           ClinVar_CLNSIG == "Benign" ~ "darkgreen",
                           ClinVar_CLNSIG == "Likely_benign" ~ "darkgreen",
                           ClinVar_CLNSIG == "Uncertain_significant" ~ "darkgrey",
                           ClinVar_CLNSIG == "Conlicting_classifications_of_pathogenicity" ~ "darkgrey"))
mutation_colors <- c("SNV" = "#fce38a", "INS" = "#f38181", "DEL" = "#95e1d3", "Other" = "#9f9f9f")
clinical_colors <- c("Pathogenic" = "red", "Likely_pathogenic" = "orange")

if(nrow(mutation_df) > 0){
  png(filename = "Figures/ideogram.png",width = 18,height = 8,units = "in",res = 600)
  gr_mutations <- GRanges(
    seqnames = mutation_df$X.CHROM,
    ranges = IRanges(start = mutation_df$POS, 
                     end = mutation_df$POS),  # only SNV
    mutation_type = mutation_df$VARIANT_CLASS,
    clinical_significance = mutation_df$ClinVar_CLNSIG
  )
  
  kp <- plotKaryotype(genome = "hg38")
  snv <- mutation_df %>% filter(VARIANT_CLASS == "SNV")
  if(nrow(snv) > 0){
    kpPoints(kp, chr=snv$X.CHROM, x=snv$POS, y=rep(-0.2, nrow(snv)), col="#fce38a", cex=2)
    kpText(kp, chr=snv$X.CHROM, x=snv$POS, y=rep(0, nrow(snv)), labels=snv$SYMBOL, pos=3, cex=0.7, col = snv$color)
  }
  
  ins <- mutation_df %>% filter(VARIANT_CLASS == "INS")
  if(nrow(ins) > 0){
    kpSegments(kp, chr=ins$X.CHROM, x0=ins$POS, x1=ins$END, y0=rep(-0.2, nrow(ins)), y1=rep(-0.2, nrow(ins)), col="#f38181", lwd=14)
    kpText(kp, chr=ins$X.CHROM, x=ins$POS, y=rep(0, nrow(ins)), labels=ins$SYMBOL, pos=3, cex=0.7, col = ins$color)
  }
  
  del <- mutation_df %>% filter(VARIANT_CLASS == "DEL")
  if(nrow(del) > 0){
    kpSegments(kp, chr=del$X.CHROM, x0=del$POS, x1=del$END, y0=rep(-0.2, nrow(del)), y1=rep(-0.2, nrow(del)), col="#95e1d3", lwd=14)
    kpText(kp, chr=del$X.CHROM, x=del$POS, y=rep(0, nrow(del)), labels=del$SYMBOL, pos=3, cex=0.7, col = del$color)
  }
  
  # 统一左侧对齐位置
  x_text <- 0.85  # 你可以调节它的位置
  
  # 每个标志物的 Y 位置
  y_start <- 0.16
  y_gap <- 0.03
  
  # Pathogenic (红色字)
  text(x=x_text, y=y_start, labels="Pathogenic", col="red", cex=1, adj=c(0,-6))
  text(x=x_text, y=y_start - y_gap, labels="Likely pathogenic", col="orange", cex=1, adj=c(0,-6))
  text(x=x_text, y=y_start, labels="Benign", col="darkgreen", cex=1, adj=c(0,0))
  text(x=x_text, y=y_start - y_gap, labels="Likely benign", col="darkgreen", cex=1, adj=c(0,0))
  text(x=x_text, y=y_start, labels="Uncertain significance", col="darkgrey", cex=1, adj=c(0,-3))
  text(x=x_text, y=y_start - y_gap, labels="Conflicting classifications", col="darkgrey", cex=1, adj=c(0,-3))
  
  
  # SNV (黄色点)
  points(x=x_text - 0.02, y=y_start - 2*y_gap, col="#fce38a", pch=16, cex=1.5)
  text(x=x_text, y=y_start - 2*y_gap, labels="SNV", cex=1, adj=c(0,0))
  
  # INS (红色线)
  segments(x0=x_text - 0.025, y0=y_start - 3*y_gap, x1=x_text - 0.015, y1=y_start - 3*y_gap, col="#f38181", lwd=2)
  text(x=x_text, y=y_start - 3*y_gap, labels="INS", cex=1, adj=c(0,0))
  
  # DEL (绿色线)
  segments(x0=x_text - 0.025, y0=y_start - 4*y_gap, x1=x_text - 0.015, y1=y_start - 4*y_gap, col="#95e1d3", lwd=2)
  text(x=x_text, y=y_start - 4*y_gap, labels="DEL", cex=1, adj=c(0,0))
  
  dev.off()
}

