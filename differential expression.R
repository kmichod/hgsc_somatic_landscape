#Prepare rna and epi/somatic files for differential expression analysis
#Author: Katheirne Lawson-Michod
#Date: October 18th, 2023

library(tidyverse)
library(DESeq2)

gene_list <- #genes identified by significantly mutated gene analysis or frequency based approaches in landscape paper 
  c("TP53", "USH2A", "LRP2", "MACF1", "RYR1", "WDFY4", "XIRP2", "APOB", "CSMD3", "DNAH7", "MGA", "OBSCN", "PIEZO2", "RYR3", "SYNE1", "UBR4", "BRCA1", "FAT3", "TACC2", "HMCN1", "DNAH3", "AHNAK", "HYDIN", "RYR2", "DST", "NF1", "BRCA2", "CDK12", "RB1", "KRAS", "GPR107", "GIMAP4", "ZNF681", "ANKRD30A")

directory <- "/Users/kayleighlawson-michod/Library/CloudStorage/OneDrive-UniversityofUtah/Doherty/AACES_Tissue_Grant/1.Writing/GitHub"
rna <- read.delim(file.path(directory, "data/rnaseq/salmon_raw_counts_for_way_pipeline_allCHR.tsv")) %>% mutate_at(vars(starts_with("Sample_")), round)
somatic_epi_merge <- readxl::read_xlsx(file.path(directory, "/data/all_data_rna.xlsx")) %>% mutate_at(vars(gene_list), as.factor) %>% filter(Study == "Black Schildkraut")

#remove "sample_" from column names in rna to map to suids
sample_columns <- grep("Sample_", colnames(rna))
new_colnames <- gsub("Sample_", "", colnames(rna)[sample_columns])
colnames(rna)[sample_columns] <- new_colnames
sampleidmapping <- read.csv(file.path(directory, "data/rnaseq/tissuegrant_epidata_07212023.csv")) %>% 
  mutate_at(vars("ID", "suid"), as.character) %>% filter(ID %in% new_colnames & race == "black")

for (i in 1:nrow(sampleidmapping)) {
  old_name <- as.character(sampleidmapping$ID[i])
  new_name <- as.character(sampleidmapping$suid[i])
  col_index <- which(names(rna) == old_name)
  if (length(col_index) > 0) {
    names(rna)[col_index] <- new_name
  }
}

#read in somatic/epi data and subset to cases included in rnaseq and tumor wes
somatic_epi_merge <- readxl::read_xlsx(file.path(directory, "/data/all_data_rna.xlsx")) %>% mutate_at(vars(gene_list), as.factor) %>% filter(Study == "Black Schildkraut")
suids <- intersect(names(rna), somatic_epi_merge$sample_id)
coldata <- somatic_epi_merge %>% filter(sample_id %in% suids) %>% mutate(TP53 = factor(TP53, levels = c(0, 1)))

#set up counts matrix for DESeq
rna <- rna %>% select(hgnc_symbol, suids)
rna_filter <- rna %>% filter(hgnc_symbol %in% gene_list) #filter to genes identified, to increase efficiency of DESeq, if you are interested in looking at other genes differntial expressed by somatic mutation remove filter
cts <- rna_filter[,-1]
rownames(cts) <- rna_filter[,1]

differential_expression <- function(genes){ #define function that applies DESeqDataSetFromMatrix() from DESeq2 for all genes identified 
  result_list <- list() #define empty list to add results to before binding 
  for(gene in genes){ 
    design_formula <- formula(paste("~", gene))
    dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = design_formula)
    dds <- DESeq(dds)
    res <- results(dds)
    res_df <- as.data.frame(res)
    res_df <- tibble::rownames_to_column(res_df, "Gene")
    res_df <- res_df %>% filter(Gene == gene)
    result_list[[gene]] <- res_df }
  merged_df <- bind_rows(result_list, .id = "Gene")
  return(merged_df)
}

differential_expression_results <- differential_expression(gene_list)
