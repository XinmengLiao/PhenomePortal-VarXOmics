library(dplyr)
library(tidyverse)
library(data.table)
library(epigraphdb)
library(httr)
library(jsonlite)
library(maps)
library(g3viz)
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
input_file <- args[1]
output_dir = args[2]

# load files
result <- read.csv(input_file,header = T,sep = "\t") %>% unique()
colnames(result) <- gsub("am_","AlphaMissense_", colnames(result))

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
if(nrow(pqtl.df)> 0){
  write.table(pqtl.df,file = gsub(".txt",".snpmr.txt", input_file),quote = F,sep = "\t",row.names = F)
}


## ---- eqtl-mr (EpigraphDB) ----
rsid <- unlist(strsplit(as.character(result$Existing_variation), split = "&"))[1]
if(grepl("^rs", rsid)){
  query_snp_mr <- function(variant_id, qtl_type = "eQTL", pval_threshold = 0.00001) {
    
    # Base URL for variant-based queries
    base_url <- "https://api.epigraphdb.org/xqtl/single-snp-mr"
    
    # Parameters
    params <- list(
      variant = variant_id,
      qtl_type = qtl_type,
      pval_threshold = pval_threshold
    )
    
    # Make request
    response <- GET(base_url, query = params)
    
    # Check response
    if (status_code(response) == 200) {
      # Parse JSON
      result <- fromJSON(content(response, "text"))
      
      # Extract results if available
      if (result$metadata$empty_results) {
        message("No results found for variant: ", variant_id)
        return(NULL)
      } else {
        return(result$results)
      }
    } else {
      stop(paste("API request failed with status:", status_code(response)))
    }
  }
  results <- query_snp_mr(rsid, "eQTL", 1)
  if( !is.null(results)){
    snp.mr <- as.data.frame(results)
    write.table(snp.mr, gsub(".txt",".snpmr.txt", input_file),quote = F,sep = "\t",row.names = F)
  }
}

## ---- AF Map ----
world_map <- map_data("world")

country_to_region <- data.frame(
  region = c(
    # South Asia (SAS)
    "India", "Pakistan", "Bangladesh", "Sri Lanka", "Nepal", "Bhutan", "Afghanistan", "Maldives",
    
    # East Asia (EAS)
    "China", "South Korea", "North Korea", "Mongolia", "Taiwan",
    
    # Africa (AFR)
    "South Africa", "Nigeria", "Egypt", "Kenya", "Ghana", "Morocco", "Algeria", "Tunisia",
    "Libya", "Sudan", "Ethiopia", "Tanzania", "Uganda", "Cameroon", "Angola", "Mozambique",
    "Madagascar", "Burkina Faso", "Mali", "Niger", "Chad", "Somalia", "Zimbabwe", "Zambia",
    "Malawi", "Rwanda", "Burundi", "Togo", "Benin", "Guinea", "Sierra Leone", "Liberia",
    "Mauritania", "Senegal", "Gambia", "Guinea-Bissau", "Cape Verde", "Sao Tome and Principe",
    "Equatorial Guinea", "Gabon", "Republic of Congo", "Democratic Republic of the Congo",
    "Central African Republic", "Botswana", "Namibia", "Lesotho", "Swaziland", "Eritrea",
    "Djibouti", "Comoros", "Mauritius", "Seychelles", "Ivory Coast",
    
    # Americas (AMR)
    "USA", "Canada", "Mexico", "Brazil", "Argentina", "Chile", "Colombia", "Venezuela",
    "Peru", "Ecuador", "Bolivia", "Paraguay", "Uruguay", "Guyana", "Suriname",
    "French Guiana", "Guatemala", "Belize", "El Salvador", "Honduras", "Nicaragua",
    "Costa Rica", "Panama", "Cuba", "Jamaica", "Haiti", "Dominican Republic",
    "Trinidad and Tobago", "Barbados", "Grenada", "Saint Lucia", "Saint Vincent and the Grenadines",
    "Antigua and Barbuda", "Dominica", "Saint Kitts and Nevis", "Bahamas",
    
    # Northern/Western Europe (NFE)
    "UK", "France", "Germany", "Spain", "Italy", "Netherlands", "Belgium", "Portugal",
    "Switzerland", "Austria", "Denmark", "Sweden", "Norway", "Ireland", "Luxembourg",
    "Iceland", "Malta", "Cyprus", "Liechtenstein", "Monaco", "San Marino", "Vatican",
    "Andorra",
    
    # Middle East (MID)
    "Saudi Arabia", "Iran", "Iraq", "Turkey", "Israel", "Syria", "Lebanon", "Jordan",
    "Kuwait", "UAE", "Qatar", "Bahrain", "Oman", "Yemen", "Palestine",
    
    # Finland (FIN)
    "Finland",
    
    # Asia Japan (ASJ)
    "Japan"
  ),
  region_code = c(
    # SAS
    rep("SAS", 8),
    # EAS
    rep("EAS", 5),
    # AFR
    rep("AFR", 54),
    # AMR
    rep("AMR", 35),
    # NFE
    rep("NFE", 23),
    # MID
    rep("MID", 15),
    # FIN
    rep("FIN", 1),
    # ASJ
    rep("ASJ", 1)
  )
)


region_mapping <- data.frame(
  region_code = c("SAS", "EAS", "AFR", "AMR", "NFE", "MID", "FIN", "ASJ"),
  region_name = c("South Asia", "East Asia", "Africa", "Americas", 
                  "Northern/Western Europe", "Middle East", "Finland", "Asia Japan")
)

world_map_with_regions <- world_map %>%
  left_join(country_to_region, by = "region") %>%
  left_join(region_mapping, by = "region_code")

af.data <- result  %>% select(which(grepl("^gnom", names(.)))) %>% select(-gnomADe_AF, -gnomADg_AF) %>% mutate(across(everything(), ~replace_na(., 0)))
af.data <- af.data %>% t() %>% as.data.frame() %>% rownames_to_column(var = "pop") %>% 
  mutate(pop = gsub("gnomADe_|gnomADg_|_AF", "", pop))
colnames(af.data) <- c("pop","AF")
af.data <- af.data %>% group_by(pop) %>% summarise(max_af = max(AF)) %>% filter(pop !="REMAINING") %>% 
  mutate(pop = toupper(pop)) 

world <- world_map_with_regions %>% left_join(., af.data, by = c("region_code" = "pop")) %>% 
  mutate(AF = if_else(is.na(max_af), 0, max_af))

p.afmap <- ggplot(world, aes(x = long, y = lat, group = group)) +
  geom_polygon(aes(fill = AF), color = "black", size = 0.1) +
  scale_fill_gradient(low = "white", high = "red", 
                      name = "AF Value",
                      labels = scales::percent_format()) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5)
  ) +
  coord_fixed(1.3)+
  theme(legend.text = element_text(size = 5),legend.title = element_text(size = 5),legend.key.size =  unit(0.6, "cm"))

ggsave(p.afmap, filename = paste0(output_dir, "/AFMap.png"),dpi = 600,width = 4)


## ---- Lollipop plot ----
user.lol.data <- result %>% 
  select(SYMBOL, X.CHROM, POS,REF, ALT, Consequence,HGVSp) %>% unique() %>% 
  filter(HGVSp !="" & !is.na(HGVSp)) %>% 
  mutate(HGVSp = str_extract(HGVSp, "p\\..*")) %>% 
  rename(Protein_Change = HGVSp, Mutation_Class = Consequence) %>% 
  mutate(Protein_Change = URLdecode(Protein_Change)) %>%
  mutate(AA_Position = as.numeric(str_extract(Protein_Change, "\\d+"))) %>% 
  mutate(Mutation_Class = sapply(strsplit(split = "&", as.character(Mutation_Class)),`[`,1)) %>% unique()

plot.options <- g3Lollipop.theme(theme.name = "simple",title.text = user.lol.data$SYMBOL,
                                 y.axis.label = paste0(user.lol.data$SYMBOL, " mutations"),
                                 legend.title = "Consequence" )
chart <- g3Lollipop(
  mutation.dat = user.lol.data,gene.symbol = user.lol.data$SYMBOL,aa.pos.col = 'AA_Position',
  protein.change.col = "Protein_Change",
  gene.symbol.col = "SYMBOL",output.filename = "/Lollipop",
  plot.options =  plot.options
  
)
htmlwidgets::saveWidget(chart, paste0(output_dir, "/Lollipop.html"))






