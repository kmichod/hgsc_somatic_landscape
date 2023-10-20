#HGSC somatic landscape paper supplemental table 1 
#Author: Katherine Lawson-Michod
#Date: July 6th, 2023
#Last updated: October 3rd, 2023

directory <- "/Users/kayleighlawson-michod/Library/CloudStorage/OneDrive-UniversityofUtah/Doherty/AACES_Tissue_Grant/1.Writing/GitHub"
source(file.path(directory, "code/2.utils.R"))
##########################################################VCF Count - Supplemental Table 1A#########################################################

vc_count <- function(variant_df, study_name){
  variant_count_filtered_df <- 
    variant_df %>%
    group_by(Variant_Classification) %>% 
    tally()
  
  total_mutations <- sum(variant_count_filtered_df$n)
  
  variant_count_filtered_df <- variant_count_filtered_df %>%
    mutate(perc = round((n/total_mutations)*100, 2))
  return(variant_count_filtered_df)
}

tissue_vc_unfiltered <- vc_count(tissue_black_snp_filtered_somatic_mutations) 
tissue_vc_unfiltered <- tissue_vc_unfiltered %>% rename(n_tissue = n, perc_tissue = perc)
tcga_vc_unfiltered <- vc_count(tcga_white_somatic_mutations_filtered) 
tcga_vc_unfiltered <- tcga_vc_unfiltered %>% rename(n_tcga= n, perc_tcga = perc)
vc_unfiltered_df <- merge(tissue_vc_unfiltered, tcga_vc_unfiltered, all.x = T, all.y = T)
tissue_vc <- vc_count(tissue_black_somatic_mutations) 
tissue_vc <- tissue_vc %>% rename(n_tissue = n, perc_tissue = perc)
tcga_vc <- vc_count(tcga_white_somatic_mutations)
tcga_vc <- tcga_vc %>% rename(n_tcga= n, perc_tcga = perc)
vc_filtered_df <- merge(tissue_vc, tcga_vc, all.x = T, all.y = T)

writexl::write_xlsx(vc_unfiltered_df, file.path(directory, "/results/S1A.xlsx"))
writexl::write_xlsx(vc_filtered_df, file.path(directory, "/results/S2A.xlsx"))
