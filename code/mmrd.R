#Look at MMRD signature 
#Author: Katherine Lawson-Michod
#Date: October 20, 2023 

directory <- "/Users/kayleighlawson-michod/Library/CloudStorage/OneDrive-UniversityofUtah/Doherty/AACES_Tissue_Grant/1.Writing/GitHub"
source(file.path(directory, "code/2.utils.R")) 

library(gtsummary)
library(gt)

sbs_tissue_black_long <- readxl::read_xlsx(file.path(directory, "/results/cosmic/sbs_tissue_black_long.xlsx"))
sbs_tissue_black <- readxl::read_xlsx(file.path(directory, "/results/cosmic/sbs_tissue_black.xlsx"))
sbs_tcga_white_nature_long <- readxl::read_xlsx(file.path(directory, "/results/cosmic/sbs_tcga_white_nature_long.xlsx"))
sbs_tcga_white_nature <- readxl::read_xlsx(file.path(directory, "/results/cosmic/sbs_tcga_white_nature.xlsx"))

# Look at MMRD signature 
sbs_tcga_white_nature_long %>% filter(eitiology == "MMRD" & activity > 0)
sbs_tissue_black_long %>% filter(eitiology == "MMRD" & activity > 0 & !signature %in% signatures_missing_tcga) %>% group_by(signature) %>% summarize(activity = sum(activity)) %>% arrange(desc(activity))
sbs_tissue_black_long %>% filter(eitiology == "MMRD" & activity > 0) %>% pull(Samples) %>% unique()

MMRD_SUIDS <- sbs_tissue_black_long %>% filter(eitiology == "MMRD" & activity > 9) %>% pull(Samples)
MMRD_neg_SUIDs <- setdiff(unique(sbs_tissue_black_long$Samples), MMRD_SUIDS)

hgsc_38_high_confidence <- hgsc_38_high_confidence %>% #Use the hgsc_38_high_confidence because it does not have intronic variants filtered
  mutate(MMRD_Sig = 
           ifelse(SUID %in% MMRD_SUIDS, "Positive",
                           ifelse(SUID %in% MMRD_neg_SUIDs, "Negtive", NA)))

MMRD_reuslts <- hgsc_38_high_confidence %>% 
  filter(Hugo_Symbol %in% c("MLH1", "MSH2", "MSH6", "PMS2")) %>% 
  group_by(MMRD_Sig) %>% 
  distinct(SUID, Hugo_Symbol, Variant_Classification) %>%
  count(Hugo_Symbol, Variant_Classification) 

write.xlsx(MMRD_reuslts, file.path(directory,"results/supplemental_materials/mmrd_comparison/somatic_mutations.xlsx"))
