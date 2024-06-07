#Process cosmic refitting results
#Author: Katherine Lawson-Michod

library(readxl)
library(dplyr)
library(tidyr)

directory <- "~/hgsc_somatic_landscape_paper/" #set working directory
source(file.path(directory, "code/2.utils.R"))

# COSMIC results for analytic samples ran using SigAssignment following matrix generation with MatrixGenerator, following FFPE adjustment
tissue_black_corrected_sbs_96 <- read.delim(file.path(directory, "programs/cosmic/SigProfilerAssignment/SigProfilerAssignment/schildkraut_black_ffpecorrected_n191/Assignment_Solution/Activities/Assignment_Solution_Activities.txt"))
tissue_black_corrected_sbs_96 <- tissue_black_corrected_sbs_96 %>% filter(Samples %in% tissue_black_suids_final)
tissue_black_corrected_sbs_96$Samples <- as.factor(tissue_black_corrected_sbs_96$Samples)

# Read in TCGA cosmic signatures for both the reprocessed data using SigAssignment and the reuslts of the Nature paper
tcga_nature_sbs_96 <- read.csv(file.path(directory, "data/tcga_cosmic/TCGA_WES_Ovary_sigProfiler_SBS_signatures_reported_from_paper.csv"))
tcga_nature_sbs_96 <- tcga_nature_sbs_96 %>% select(-c("Cancer.Types"))
tcga_nature_sbs_96 <- tcga_nature_sbs_96 %>% mutate(., Sample.Names = gsub('.{16}$', '', Sample.Names))
tcga_white_nature_sbs_96 <- tcga_nature_sbs_96 %>% filter(., Sample.Names %in% tcga_white_cases_under_79_ids)
tcga_white_nature_sbs_96 <- tcga_white_nature_sbs_96 %>% rename(Samples = Sample.Names) %>% select(-c(Accuracy))
length(unique(tcga_white_nature_sbs_96$Samples)) #176 cases

# Define levels for the different signatures and recode signature with eitiology information 
eitiologySBS96.levels <- list(
  Age = c('SBS1'),
  MMRD = c('SBS6', 'SBS14', 'SBS15', 'SBS20', 'SBS21', 'SBS26', 'SBS44'),
  POLE = c('SBS10a', 'SBS10b', 'SBS10c', 'SBS10d', 'SBS28'),
  HRD = c('SBS3'), 
  BER = c('SBS30', 'SBS36'), 
  Treatment_Chemotherapy = c("SBS11", "SBS25", "SBS31", "SBS35", "SBS86", "SBS87", "SBS90"),
  Treatment_Immunosuppressants = c("SBS32"),
  APOBEC = c("SBS2", "SBS13"),
  Tobacco = c("SBS4", "SBS29", "SBS92"),
  UV = c("SBS7a", "SBS7b", "SBS7c", "SBS7d", "SBS38"),
  AA = c("22"),
  Colibactin = c("88"),
  Artifact = c("SBS27", "SBS43", "SBS45", "SBS46", "SBS47", "SBS48", "SBS49", "SBS50", "SBS51", "SBS52", "SBS53", "SBS54", "SBS55", "SBS56", "SBS57", "SBS58", "SBS59", "SBS60", "SBS95"),
  Lymphoid = c("SBS9", "SBS84", "SBS85")
)

recode_eitiology_sbs <- function(value) {
  case_when(
    value %in% eitiologySBS96.levels$Age ~ "Age",
    value %in% eitiologySBS96.levels$MMRD ~ "MMRD",
    value %in% eitiologySBS96.levels$POLE ~ "POLE",
    value %in% eitiologySBS96.levels$HRD ~ "HRD",
    value %in% eitiologySBS96.levels$BER ~ "BER",
    value %in% eitiologySBS96.levels$Treatment_Chemotherapy ~ "Chemotherapy",
    value %in% eitiologySBS96.levels$Treatment_Immunosuppressants ~ "Immunosuppressants",
    value %in% eitiologySBS96.levels$APOBEC ~ "APOBEC",
    value %in% eitiologySBS96.levels$Tobacco ~ "Tobacco",
    value %in% eitiologySBS96.levels$UV ~ "UV",
    value %in% eitiologySBS96.levels$AA ~ "AA",
    value %in% eitiologySBS96.levels$Colibactin ~ "Colibactin",
    value %in% eitiologySBS96.levels$Artifact ~ "Artifact",
    value %in% eitiologySBS96.levels$Lymphoid ~ "Lymphoid",
    TRUE ~ NA_character_  # Return NA for unmatched values
  )
}

# Write function that takes counts of mutations attributed to each signature and calculate the proportion of the total tumor mutation burden attributed to each signature and annotate with known eitiology information   
sbs96_longer = function(cosmic_counts_df){
  cols = names(cosmic_counts_df)[-1]
  df_long <- 
    cosmic_counts_df %>% 
    pivot_longer(cols, names_to = "signature", values_to = "activity") %>%
    mutate(eitiology = signature) %>%
    mutate(eitiology = recode_eitiology_sbs(eitiology)) %>%
    #mutate(eitiology = coalesce(eitiology, signature))
  return(df_long)
}

sbs96_counts_to_tmb = function(df_long){
  study <-  gsub("_long","", deparse(substitute(df_long)))
  df_eitiology <-
    df_long %>% 
    group_by(eitiology) %>% 
    summarize(mutation_count = sum(activity)) %>% 
    mutate(total_tmb = sum(mutation_count)) %>%
    mutate(., perc = (mutation_count/total_tmb)*100)
  df_eitiology <- df_eitiology %>%
    rename(., !!paste0("mutation_count_", study) := mutation_count, !!paste0("total_tmb_", study) := total_tmb, !!paste0("perc_", study) := perc)
  return(df_eitiology)
}

sbs96_counts_to_tmb_signature = function(df_long){
  study <-  gsub("_long","", deparse(substitute(df_long)))
  df_signature <-
    df_long %>% 
    group_by(signature) %>% 
    summarize(mutation_count = sum(activity)) %>% 
    mutate(total_tmb = sum(mutation_count)) %>%
    mutate(., perc = (mutation_count/total_tmb)*100)
  df_signature <- df_signature %>%
    rename(., !!paste0("mutation_count_", study) := mutation_count, !!paste0("total_tmb_", study) := total_tmb, !!paste0("perc_", study) := perc)
  return(df_signature)
}

# Apply function  
sbs_tissue_black_long <- sbs96_longer(tissue_black_corrected_sbs_96)
sbs_tissue_black_long_exclude_mmrd <- sbs_tissue_black_long %>% filter(!signature %in% c('SBS6', 'SBS14', 'SBS15', 'SBS20', 'SBS21', 'SBS26', 'SBS44'))
sbs_tissue_black <- sbs96_counts_to_tmb(sbs_tissue_black_long_exclude_mmrd)
sbs_tissue_black_signature <- sbs96_counts_to_tmb_signature(sbs_tissue_black_long_exclude_mmrd)

sbs_tcga_white_nature_long <- sbs96_longer(tcga_white_nature_sbs_96)
sbs_tcga_white_nature_long_exclude_mmrd <- sbs_tcga_white_nature_long %>% filter(!signature %in% c('SBS6', 'SBS14', 'SBS15', 'SBS20', 'SBS21', 'SBS26', 'SBS44'))
sbs_tcga_white_nature <- sbs96_counts_to_tmb(sbs_tcga_white_nature_long_exclude_mmrd)
sbs_tcga_white_nature_signature <- sbs96_counts_to_tmb_signature(sbs_tcga_white_nature_long_exclude_mmrd)

chemo <- sbs_tissue_black_long %>% filter(eitiology == "Chemotherapy" & activity > 0)
chemo_neo <- chemo %>% filter(Samples %in% suids_neo_yes)
chemo_no_neo <- chemo %>% filter(Samples %in% suids_neo_no)
chemo_no_neo <- chemo %>% filter(Samples %in% suids_neo_missing)

chemo_sum <- chemo %>% group_by(Samples, eitiology) %>% summarise(activity = sum(activity))
hist(chemo_no_neo$activity)
hist(chemo_no_neo$activity, 
     main = "Chemotherapy Signature Distribution in Schildkraut-B Missing", 
     xlab = "Number of Somatic Mutations Per Sample", 
     ylab = "Count of Samples")

hist(chemo_no_neo$activity)
hist(chemo_neo$activity)

chemo_no_neo <- chemo_no_neo %>%
  group_by(signature) %>% summarise(median(activity), 
                                            min(activity), 
                                            max(activity))
chemo_neo <- chemo_neo %>%
  group_by(signature) %>% summarise(median(activity), 
                                    min(activity), 
                                    max(activity))
writexl::write_xlsx(chemo_no_neo, file.path(directory, "/neoadj/chemo_no_neo.xlsx"))
writexl::write_xlsx(chemo_neo, file.path(directory, "/neoadj/chemo_neo.xlsx"))
writexl::write_xlsx(sbs_tissue_black_long, file.path(directory, "/results/cosmic/sbs_tissue_black_long.xlsx"))
writexl::write_xlsx(sbs_tissue_black, file.path(directory, "/results/cosmic/sbs_tissue_black.xlsx"))
writexl::write_xlsx(sbs_tcga_white_nature_long, file.path(directory, "/results/cosmic/sbs_tcga_white_nature_long.xlsx"))
writexl::write_xlsx(sbs_tcga_white_nature, file.path(directory, "/results/cosmic/sbs_tcga_white_nature.xlsx"))

# Define binary variable for whether or not someone has signature based on 10 mutations attributed to a given signature threshold
tissue_cosmic <- 
  sbs_tissue_black_long %>% 
  filter(activity > 9) %>% 
  distinct(Samples, eitiology) %>%
  group_by(Samples, eitiology) %>% 
  tally() %>%
  pivot_wider(names_from = eitiology, values_from = n) %>%
  mutate_all(~ifelse(is.na(.), 0, .)) %>%
  mutate(Study = "Schildkraut-B") %>%
  rename("sample_id" = Samples) %>%
  mutate_at(vars("sample_id", "Study"), as.character)
tissue_cosmic_signature <- 
  sbs_tissue_black_long %>% 
  filter(activity > 9) %>% 
  distinct(Samples, signature) %>%
  group_by(Samples, signature) %>% 
  tally() %>%
  pivot_wider(names_from = signature, values_from = n) %>%
  mutate_all(~ifelse(is.na(.), 0, .)) %>%
  mutate(Study = "Schildkraut-B") %>%
  rename("sample_id" = Samples) %>%
  mutate_at(vars("sample_id", "Study"), as.character)
tissue_cosmic <- merge(tissue_cosmic, tissue_cosmic_signature, by = c("sample_id", "Study"))
tissue_add <- setdiff(unique(sbs_tissue_black_long$Samples), tissue_cosmic$sample_id)
tissue_cosmic_names <- names(tissue_cosmic)
col_names <- setdiff(tissue_cosmic_names, c("sample_id", "Study"))
length(col_names)
zeros_list <- rep(0, 63)

# Create a new dataframe with the desired rows
tissue_add_df <- data.frame(
  sample_id = tissue_add,
  Study = "Schildkraut-B",
  setNames(data.frame(matrix(zeros_list, ncol = length(col_names))), col_names)
)

colnames(tissue_add_df) <- colnames(tissue_cosmic)
tissue_cosmic_final <- rbind(tissue_cosmic, tissue_add_df)

tcga_cosmic <- 
  sbs_tcga_white_nature_long %>% 
  filter(activity > 9) %>% 
  distinct(Samples, eitiology) %>%
  group_by(Samples, eitiology) %>% 
  tally() %>%
  pivot_wider(names_from = eitiology, values_from = n) %>%
  mutate_all(~ifelse(is.na(.), 0, .)) %>%
  mutate(Study = "TCGA-W") %>%
  rename("sample_id" = Samples) %>%
  mutate_at(vars("sample_id", "Study"), as.character)
tcga_cosmic_signature <- 
  sbs_tcga_white_nature_long %>% 
  filter(activity > 9) %>% 
  distinct(Samples, signature) %>%
  group_by(Samples, signature) %>% 
  tally() %>%
  pivot_wider(names_from = signature, values_from = n) %>%
  mutate_all(~ifelse(is.na(.), 0, .)) %>%
  mutate(Study = "TCGA-W") %>%
  rename("sample_id" = Samples) %>%
  mutate_at(vars("sample_id", "Study"), as.character)
tcga_cosmic <- merge(tcga_cosmic, tcga_cosmic_signature, by = c("sample_id", "Study"))
tcga_add <- setdiff(unique(sbs_tcga_white_nature_long$Samples), tcga_cosmic$sample_id)
tcga_cosmic_names <- names(tcga_cosmic)
col_names <- setdiff(tcga_cosmic_names, c("sample_id", "Study"))
length(col_names)
zeros_list <- rep(0, 38)

# Create a new dataframe with the desired rows
tcga_add_df <- data.frame(
  sample_id = tcga_add,
  Study = "TCGA-W",
  setNames(data.frame(matrix(zeros_list, ncol = length(col_names))), col_names)
)
colnames(tcga_add_df) <- colnames(tcga_cosmic)
tcga_cosmic_final <- rbind(tcga_cosmic, tcga_add_df)

# Create a new dataframe with the desired rows
all_columns <- unique(c(colnames(tissue_cosmic_final), colnames(tcga_cosmic_final)))


# Add missing columns and fill with 0 in tissue_cosmic_final
missing_columns_tissue <- setdiff(all_columns, colnames(tissue_cosmic_final))
tissue_cosmic_final[, missing_columns_tissue] <- 0

# Add missing columns and fill with 0 in tcga_cosmic_final
missing_columns_tcga <- setdiff(all_columns, colnames(tcga_cosmic_final))
tcga_cosmic_final[, missing_columns_tcga] <- 0

cosmic <- rbind(tissue_cosmic_final, tcga_cosmic_final)
writexl::write_xlsx(cosmic, file.path(directory, "/results/cosmic/cosmic.xlsx"))

