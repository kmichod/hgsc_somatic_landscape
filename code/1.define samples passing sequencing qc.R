#Define samples passing sequencing qc
#Author: Katherine Lawson-Michod

library(tidyverse)
library(dplyr)

directory <- "~/hgsc_somatic_landscape_paper/" #set working directory

#Read in qc stats for tumor WES
qc_stats_black_cases <- readxl::read_xlsx(file.path(directory, "reference_data/qc_stats/qcStats17Feb2022_new_old_master.xlsx"))
qc_stats_black_cases <- qc_stats_black_cases %>% mutate_at(vars(`Coverage at 0.95 of Target BPs`), as.numeric) %>% mutate_at(vars(`Concordance Check`), as.factor)

qc_stats_black_cases <- qc_stats_black_cases %>% mutate(qc_coverage_pass = factor(ifelse(`Coverage at 0.95 of Target BPs` > 11, "PASS", "FAIL"))) #Define variable for whether coverage at 95% was equal to or greater than 12X
qc_stats_black_cases <- qc_stats_black_cases %>% mutate(., `Concordance Check`=recode(`Concordance Check`, 'FAIL, OK'='PASS', 'OK'='PASS')) #Recode concordance check okay = pass
qc_stats_black_cases <- qc_stats_black_cases %>% mutate(., `qc_pass_overall`= as.factor(ifelse(qc_coverage_pass == "PASS" & `Concordance Check` == "PASS", "PASS", "FAIL"))) #Define overall pass variable which is pass if both qc coverage and concordance check are pass, otherwise samples fail

#Define indicator variable for duplicate cases (cases that were re-sequenced and have two rounds of qc)
qc_stats_black_cases$duplicate <- as.character(duplicated(qc_stats_black_cases$SUID)) 
qc_stats_black_cases_duplicate <- qc_stats_black_cases %>% filter(duplicate==TRUE)
duplicate_suids <- c(as.character(unique(qc_stats_black_cases_duplicate$SUID))) #these are the suids for individuals who were re-sequenced n=214 following the full set

#Define sample IDs for pass and fail as lists of factors
duplicates_df <- as.data.frame(qc_stats_black_cases %>% group_by(SUID, qc_pass_overall) %>% summarise(count=n()))
suids_duplicate_failed <- duplicates_df %>% filter(qc_pass_overall=='FAIL' & count == 2) %>% pull(SUID) %>% unique() %>% as.character() #define suid list for individuals who were resequenced and failed both rounds
suids_duplicate_pass <- qc_stats_black_cases_duplicate %>% filter(!SUID %in% suids_duplicate_failed) %>% pull(SUID) %>% unique() %>% as.character() #suids for individuals who were resequenced and past at least one round
suids_not_duplicate_failed <- qc_stats_black_cases %>% filter(!SUID %in% duplicate_suids & qc_pass_overall=='FAIL') %>% pull(SUID) %>% unique() %>% as.character()
suids_not_duplicate_pass <- qc_stats_black_cases %>% filter(!SUID %in% duplicate_suids & qc_pass_overall=='PASS') %>% pull(SUID) %>% unique() %>% as.character()
pass_sample_qc_tumor <- append(suids_duplicate_pass, suids_not_duplicate_pass) #suids for people passing tumor WES qc
fail_sample_qc_tumor <- append(suids_duplicate_failed, suids_not_duplicate_failed) #suids for people failing tumor WES qcs

fail_20 <- qc_stats_black_cases %>% filter(SUID %in% fail_sample_qc_tumor)
duplicates <- qc_stats_black_cases %>% filter(SUID %in% duplicate_suids)
original <- qc_stats_black_cases %>% filter(`Ori or New` == "Ori")
new <- qc_stats_black_cases %>% filter(`Ori or New` == "New")
hist(original$`Coverage at 0.95 of Target BPs`)
hist(new$`Coverage at 0.95 of Target BPs`)

#read in pilot data to look at overlap with full sequencing and resequencing
pilot <- readxl::read_excel(file.path(directory, "reference_data/qc_stats/Suids for tissue WES 12112019.xlsx"), sheet = "Pilot Cases - WES complete")
pilot_qc <- readxl::read_xlsx(file.path(directory, "reference_data/qc_stats/readme_qcReport_Stats_withComments.xlsx"), sheet = "Stats V2 design - input")
pilot_suids <- as.character(unique(pilot$suid))
full_seq <- append(pass_sample_qc_tumor, fail_sample_qc_tumor)
length(unique(full_seq))
all_seq <- unique(append(as.character(full_seq), as.character(pilot_suids)))
setdiff(pilot_suids, full_seq)

#Read in qc stats for normal WES
normal <- readxl::read_excel(file.path(directory, "reference_data/qc_stats/qcReport_Stats_TN_Ovarian_original.xlsx"), sheet = "Merged With Concordance")
normal <- normal %>% filter(Source == "Normal") %>% 
  mutate_at(vars('# Fastq Reads', '# Unfiltered Alignments', 'Fraction Duplicate', 'Mean Insert Size', 'Fraction Overlapping BPs', '# Unique BPs', 
                 'Fraction Q20 BPs', 'Fraction Q30 BPs', 'Mean on Target Coverage', 'Coverage at 0.9 of Target BPs', 'Coverage at 0.95 of Target BPs'), as.numeric) %>%
  mutate_at(vars('Bam Concordance Check', 'SUID for TN Pairs', 'Source'), as.factor)
normal <- normal %>% filter(!is.na(`SUID for TN Pairs`))
length(normal$`SUID for TN Pairs`)

normal <- normal %>% mutate(qc_coverage_pass = as.factor(ifelse(`Coverage at 0.9 of Target BPs` > 11, "PASS", "FAIL"))) #Define variable for whether coverage at 90% was greater than or equal to 12X
normal <- normal %>% mutate(., `Concordance Check`=recode(`Bam Concordance Check`, 'FAIL, OK'='PASS', 'OK'='PASS')) #Recode concordance check so that the ok = pass - confirm with david 
normal <- normal %>% mutate(., `qc_pass_overall`=as.factor(ifelse(qc_coverage_pass == "PASS" & `Concordance Check` == "PASS", "PASS", "FAIL"))) #Define overall pass variable which is pass if both qc coverage and concordance check are pass, otherwise they fail

#Define reasons for exclusion - matches supplemental figure 1
cord_fail <- qc_stats_black_cases %>% filter(., SUID %in% fail_sample_qc_tumor & `Concordance Check` == "FAIL") %>% pull (SUID) %>% unique() %>% as.character()
cov_fail_tumor <- qc_stats_black_cases %>% filter(., SUID %in% fail_sample_qc_tumor & qc_coverage_pass == "FAIL") %>% pull (SUID) %>% unique() %>% as.character()
cov_fail_normal <- normal %>% filter(., qc_coverage_pass == "FAIL") %>% pull ("SUID for TN Pairs") %>% unique() %>% as.character()
temp <- append(as.character(cord_fail), as.character(cov_fail_tumor))
temp_2 <- append(as.character(temp), as.character(cov_fail_normal))
hypermutated_suid <- c(150015, 110209) #sample 150015 and 110209 showed an extreme hyper mutated phenotype remove from analytic files
tissue_fail_suids <- unique(append(temp_2, as.character(hypermutated_suid)))
tissue_fail_suids_df <- as.data.frame(tissue_fail_suids)
tissue_pass_suids <- unique(setdiff(all_seq, tissue_fail_suids))
tissue_pass_suids_df <- as.data.frame(tissue_pass_suids)
length(tissue_pass_suids) + length(tissue_fail_suids)
   
#Counts for each qc exclusion reason in supplemental figure 1
setdiff(cord_fail, append(cov_fail_tumor, cov_fail_normal))
intersect(cord_fail, cov_fail_tumor)
intersect(cord_fail, cov_fail_normal)
length(setdiff(cov_fail_tumor, append(cord_fail, cov_fail_normal)))
intersect(cov_fail_tumor, cov_fail_normal)
length(setdiff(cov_fail_normal, append(cord_fail, cov_fail_tumor)))

#Save suids pass/fail files to read into utilities file
write_delim(tissue_fail_suids_df, file.path(directory, "reference_data/fail_qc_suids_somatic_landscape.txt"))
write_delim(tissue_pass_suids_df, file.path(directory, "reference_data/pass_qc_suids_somatic_landscape.txt"))
