#Merge Epi, Somatic, RNASeq data
#Author: Katherine Lawson-Michod

directory <- "~/hgsc_somatic_landscape_paper/" #set working directory
source(file.path(directory, "code/2.utils.R")) 


###Somatic
#Define binary variable for whether or no someone has a mutation in a gene
distinct_mut_tissue <- tissue_black_somatic_mutations %>% 
  distinct(SUID, Hugo_Symbol) %>% #Only interested in any or none, need to remove if you are interested in multiple mutations in the same gene in the same person
  group_by(SUID, Hugo_Symbol) %>% 
  tally() %>%
  pivot_wider(names_from = Hugo_Symbol, values_from = n) %>%
  mutate_all(~ifelse(is.na(.), 0, .)) %>%
  mutate(Study = "Schildkraut-B") %>%
  rename("sample_id" = SUID) %>%
  mutate_at(vars("sample_id", "Study"), as.character)

distinct_mut_tcga <- tcga_white_somatic_mutations %>% 
  distinct(case_submitter_id, Hugo_Symbol) %>% 
  group_by(case_submitter_id, Hugo_Symbol) %>% 
  tally() %>%
  pivot_wider(names_from = Hugo_Symbol, values_from = n) %>%
  mutate_all(~ifelse(is.na(.), 0, .)) %>%
  mutate(Study = "TCGA-W") %>%
  rename("sample_id" = case_submitter_id) %>%
  mutate_at(vars("sample_id", "Study"), as.character)

somatic <- rbind(distinct_mut_tissue, distinct_mut_tcga) #462

###COSMIC 
cosmic <- readxl::read_xlsx(file.path(directory, "/results/cosmic/cosmic.xlsx")) %>%
  mutate_at(vars("sample_id", "Study"), as.character) #366

###EPI
## Schildkraut-B  
# Recode clinical variables for table 1
tissue_epi <- tissue_black_epi_final %>% 
  filter(suid %in% tissue_black_suids_final) %>%
  select(., c(suid, dobyear, diagyear, vitalstatus, cancersite, 
              histology, behavior, stage, refage, hispanic, birthplace, 
              education, married, race, timelastfu, neoadj_treat)) %>%
  mutate_at(., c('vitalstatus', 'cancersite', 'histology', 'behavior', 'stage', 'hispanic', 'birthplace', 'education', 'married', 'neoadj_treat'), as.factor) %>%
  filter(stage != 1) %>% # remove stage 1 because excluded from tcga-w
  mutate(Study = "Schildkraut-B", 
         vitalstatus=recode(vitalstatus, `1`="Alive",`2`="Deceased"),
         cancersite=recode(cancersite, `1`="Ovarian",`2`="Tubal", `3`="Peritoneal",`4`="Ovarian or tubal"), 
         histology=recode(histology, `1`="Serous",`2`="Endometrioid"),
         behavior=recode(behavior, `2`="Invasive"),
         hispanic=recode(hispanic, `1`="Hispanic",`2`="Not Hispanic", `99`="Missing"), 
         birthplace=recode(birthplace, `1`="Born in the US",`2`="Born outside of the US", `99`="Missing"), 
         education=recode(education, `1`="High school graduate/GED or less",`2`="Some college", `3`="College graduate", `4`="Graduate/professional school", `99`="Missing"),
         married=recode(married, `1`="Single/never married",`2`="Married/living as married", `3`="Divorced/separated", `4`="Widowed", `99`="Missing"),
         Race=recode(race, `2`="Black"),
         "Neoadjuvant treatment" = recode(neoadj_treat, `1` = "Yes", `2` = "No"),
         stage_b = ifelse(stage %in% c(1,2), "FIGO II", ifelse(stage %in% c(3), "FIGO III and IV", ifelse(stage == 99, "Missing", NA))),
         dobyear_10yr = ifelse(dobyear < 1940, "< 1940",
                               ifelse(dobyear >= 1940 & dobyear < 1950, "1940-1949",
                                      ifelse(dobyear >= 1950 & dobyear < 1960, "1950-1959",
                                             ifelse(dobyear >= 1960, "≥ 1960", NA)))),
         diagyear_10yr = ifelse(diagyear < 2000, "< 2000", 
                                ifelse(diagyear > 1999 & diagyear < 2005, "2000-2004",
                                       ifelse(diagyear > 2004 & diagyear < 2010, "2005-2009",
                                              ifelse(diagyear > 2009, "2010-2015", NA)))), 
         refage_5_year = ifelse(refage > 34 & refage < 40, '35-39', 
                                ifelse(refage > 39 & refage < 45, '40-44', 
                                       ifelse(refage > 44 & refage < 50, '45-49', 
                                              ifelse(refage > 49 & refage < 55, '50-54', 
                                                     ifelse(refage > 54 & refage < 60, '55-59', 
                                                            ifelse(refage > 59 & refage < 65, '60-64', 
                                                                   ifelse(refage > 64 & refage < 70, '65-69',
                                                                          ifelse(refage > 69, '≥ 70', NA))))))))) %>%
  mutate_at(., c('dobyear_10yr', 'diagyear_10yr', 'stage_b'), as.factor) %>%
  rename("Vital status" = vitalstatus, 
         "Cancer site" = cancersite, 
         "Histology" = histology, 
         "Behavior" = behavior, 
         "Age at diagnosis" = refage, 
         "Ethnicity" = hispanic, 
         "Birth place" = birthplace, 
         "Education" = education, 
         "Married" = married, 
         "Year of birth" = dobyear_10yr,
         "Year of diagnosis" = diagyear_10yr, 
         "Stage" = stage_b,
         sample_id = suid) %>% 
  select(., -c(stage, Married, Education, Behavior, `Birth place`, race, "Cancer site", "Histology", neoadj_treat)) #remove extraneous variables prior to tcga merge

#TCGA-W
tcga_epi <- tcga_white_cases_under_79 %>% 
  select(., c(days_to_death, case_submitter_id, figo_stage, age_at_diagnosis, year_of_diagnosis, ethnicity, race, vital_status, year_of_birth, days_to_last_follow_up, figo_stage, site_of_resection_or_biopsy, tissue_or_organ_of_origin, primary_diagnosis)) %>%
  mutate('Age at diagnosis' = round(age_at_diagnosis/365, 0), 
         'Vital status'=recode(vital_status, `Alive`="Alive",`Dead`="Deceased"),
         Ethnicity=recode(ethnicity, `hispanic or latino`="Hispanic",`not hispanic or latino`="Not Hispanic", `not reported`= "Missing"), 
         # 'Cancer site'=recode(site_of_resection_or_biopsy, `Specified parts of peritoneum`="Peritoneal",`Ovary`="Ovarian"), 
         'Race'=recode(race, `white`="White"), 
         # 'Histology'=recode(primary_diagnosis, `Serous cystadenocarcinoma, NOS`="Serous"),
         'Stage'=ifelse(figo_stage %in% c("Stage IIIA", "Stage IIIB", "Stage IIIC", "Stage IV"), "FIGO III and IV",
                        ifelse(figo_stage %in% c("Stage IIA", "Stage IIB", "Stage IIC"), "FIGO II", "Missing")),
         'Study' = "TCGA-W",
         'Neoadjuvant treatment' = as.factor("No"),
         days_to_last_follow_up = recode(days_to_last_follow_up, "'--" = NA_character_)) %>%
  mutate_at(vars(vital_status, ethnicity, race, figo_stage, site_of_resection_or_biopsy, tissue_or_organ_of_origin), as.factor) %>%
  mutate_at(vars(days_to_last_follow_up, year_of_birth, year_of_diagnosis, age_at_diagnosis), as.numeric) %>%
  mutate(diagyear_10yr = ifelse(year_of_diagnosis < 2000, "< 2000",
                                ifelse(year_of_diagnosis > 1999 & year_of_diagnosis < 2005, "2000-2004",
                                       ifelse(year_of_diagnosis > 2004 & year_of_diagnosis < 2010, "2005-2009",
                                              ifelse(year_of_diagnosis > 2009, "2010-2015", NA)))), 
         dobyear_10yr = ifelse(year_of_birth < 1940, "< 1940",
                               ifelse(year_of_birth >= 1940 & year_of_birth < 1950, "1940-1949",
                                      ifelse(year_of_birth >= 1950 & year_of_birth < 1960, "1950-1959",
                                             ifelse(year_of_birth >= 1960, "≥ 1960", NA)))), 
         refage_5_year = ifelse(`Age at diagnosis` > 34 & `Age at diagnosis` < 40, '35-39', 
                                ifelse(`Age at diagnosis` > 39 & `Age at diagnosis` < 45, '40-44', 
                                       ifelse(`Age at diagnosis` > 44 & `Age at diagnosis` < 50, '45-49', 
                                              ifelse(`Age at diagnosis` > 49 & `Age at diagnosis` < 55, '50-54', 
                                                     ifelse(`Age at diagnosis` > 54 & `Age at diagnosis` < 60, '55-59', 
                                                            ifelse(`Age at diagnosis` > 59 & `Age at diagnosis` < 65, '60-64', 
                                                                   ifelse(`Age at diagnosis` > 64 & `Age at diagnosis` < 70, '65-69',
                                                                          ifelse(`Age at diagnosis` > 69, '≥ 70', NA))))))))) %>%
  mutate_at(vars(diagyear_10yr, dobyear_10yr, refage_5_year), as.factor) %>%
  rename('Year of diagnosis' = diagyear_10yr, 'Year of birth' = dobyear_10yr, diagyear = year_of_diagnosis, dobyear = year_of_birth) %>%
  select(-c(age_at_diagnosis, ethnicity, figo_stage, site_of_resection_or_biopsy, 
            race, tissue_or_organ_of_origin, primary_diagnosis)) %>% 
  filter(`Age at diagnosis` < 79 | is.na(`Age at diagnosis`)) %>% 
  rename(sample_id = case_submitter_id) %>% 
  mutate(timelastfu = ifelse(vital_status == "Alive", days_to_last_follow_up,
                             ifelse(vital_status == "Dead", days_to_death, NA))) %>% 
  select(-c("days_to_last_follow_up", "days_to_death", "vital_status"))
epi <- rbind(tissue_epi, tcga_epi) #462

###RNA
#Read in Waypipeline cluster assignments
rnaseq_raw <- read_csv(file.path(directory, "/data/rnaseq/FullClusterMembership.csv")) %>% 
  filter(Dataset %in% c("aaces.rnaseq.eset", "TCGA")) %>% 
  select(c("Dataset", "...1", "ClusterK4_kmeans"))

#Read in dataframe for mapping SUIds to IDs
sampleidmapping <- read.csv(file.path(directory, "data/rnaseq/tissuegrant_epidata_07212023.csv")) %>% select(c("ID", "suid"))

#restructure sample ID for TCGA to match somatic and epi files and map schildkraut-B to suids
rnaseq_rename <- rnaseq_raw %>% 
  mutate(sample_id = ifelse(Dataset == "TCGA", gsub("\\.", "-", `...1`), 
                           ifelse(Dataset == "aaces.rnaseq.eset", gsub("^Sample_", "", `...1`), `...1`))) %>%
  left_join(sampleidmapping, by = c("sample_id" = "ID")) %>%
  mutate(sample_id = ifelse(
    Dataset == "aaces.rnaseq.eset", suid, sample_id), 
    ClusterK4_kmeans = recode(ClusterK4_kmeans, "1" = "C1.MES", "2" = "C5.PRO", "3" = "C2.IMM", "4" = "C4.DIF"),
    Dataset = recode(Dataset, "aaces.rnaseq.eset" = "Schildkraut-B", "TCGA" = "TCGA-W")) %>%
  rename(Study = Dataset) %>%
  select(c(Study, ClusterK4_kmeans, sample_id)) %>%
  mutate_at(vars("Study", "sample_id"), as.character)


somatic_select <- somatic %>% select(Study, sample_id, c("TP53", "USH2A", "LRP2", "MACF1", "RYR1", "WDFY4", "XIRP2", "APOB", "CSMD3", "DNAH7", "MGA", "OBSCN", "PIEZO2", "RYR3", "SYNE1", "UBR4", "BRCA1", "FAT3", "TACC2", "HMCN1", "DNAH3", "AHNAK", "HYDIN", "RYR2", "DST", "NF1", "BRCA2", "CDK12", "RB1", "KRAS", "GPR107", "GIMAP4", "ZNF681", "ANKRD30A"))
merge1 <- merge(epi, somatic_select)
merge2 <- merge(merge1, cosmic, all.x = TRUE)
merge_3 <- merge(merge2, rnaseq_rename, by = c("Study", "sample_id"), all.x = TRUE, all.y = FALSE) 
merge_3$duplicate <- duplicated(merge_3$sample_id)
all_data <- merge_3 %>% filter(duplicate == "FALSE")

writexl::write_xlsx(all_data, file.path(directory, "results/all_results.xlsx"))