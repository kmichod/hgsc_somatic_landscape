#Differential expression for genes identified in the somatic landscape paper
#Author: Katherine Lawson-Michod
#Date: October 18th, 2023

#BiocManager::install("DESeq2")
library(DESeq2)
library(tidyverse)

rna <- read.delim(file.path(directory, "/data/salmon_raw_counts_for_way_pipeline_allCHR_rename.tsv")) #read in salmon raw counts, the file sample names have been renamed to match suids in the somatic/epi data, and the counts have been rounded
coldata <- readxl::read_xlsx(file.path(directory, "/data/all_data_rna_filtered.xlsx")) #read in somatic/epi data, will use for defining coniditon, has or does not have somatic mutation in gene of interest

rna_filter <- rna %>% filter(hgnc_symbol %in% gene_list) #filter to genes identified, to increase efficiency of DESeq, if you are interested in looking at other genes differntial expressed by somatic mutation remove filter
cts <- rna_filter[,-1]
rownames(cts) <- rna_filter[,1]

gene_list <- #genes identified by significantly mutated gene analysis or frequency based approaches in landscape paper 
  c("USH2A", "LRP2", "MACF1", "RYR1", "WDFY4", "XIRP2", "APOB", "CSMD3", "DNAH7", "MGA", "OBSCN", "PIEZO2", "RYR3", "SYNE1", "UBR4", "BRCA1", "FAT3", "TACC2", "HMCN1", "DNAH3", "AHNAK", "HYDIN", "RYR2", "DST", "NF1", "BRCA2", "CDK12", "RB1", "KRAS", "GPR107", "GIMAP4", "ZNF681", "ANKRD30A")

differential_expression <- function(genes){ #define function that applies DESeqDataSetFromMatrix() from DESeq2 for all genes identified 
  result_list <- list() #define empty list to add results to before binding 
  for(gene in genes){ 
    design_formula <- formula(paste("~", gene))
    dds <- DESeqDataSetFromMatrix(countData = count_mat, colData = coldata, design = design_formula)
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
