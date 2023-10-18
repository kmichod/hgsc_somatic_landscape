#Prepare rna and epi/somatic files for differential expression analysis
#Author: Katheirne Lawson-Michod
#Date: October 18th, 2023

library(tidyverse)
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

somatic_epi_merge <- readxl::read_xlsx(file.path(directory, "/data/all_data_rna.xlsx")) %>% mutate_at(vars(gene_list), as.factor) %>% filter(Study == "Black Schildkraut")
suids <- intersect(names(rna), somatic_epi_merge$sample_id)
rna <- rna %>% select(hgnc_symbol, suids)
rna$duplicate <- duplicated(rna$hgnc_symbol)
rna <- rna %>% filter(duplicate != TRUE) %>% select(-c(duplicate))

somatic_epi_merge <- somatic_epi_merge %>% filter(sample_id %in% suids)
somatic_epi_merge$TP53 <- factor(somatic_epi_merge$TP53, levels = c(0, 1))
write.delim(rna, file.path(directory, "/data/salmon_raw_counts_for_way_pipeline_allCHR_rename.tsv"))
writexl::write_xlsx(somatic_epi_merge, file.path(directory, "/data/all_data_rna_filtered.xlsx"))

