#Cohort characteristics
#Author: Katherine Lawson-Michod
#Date: July 5th, 2023
#Last Updated: Updated October 20th, 2023

directory <- "/Users/kayleighlawson-michod/Library/CloudStorage/OneDrive-UniversityofUtah/Doherty/AACES_Tissue_Grant/1.Writing/GitHub"
source(file.path(directory, "code/2.utils.R")) 

library(gtsummary)
library(gt)

all_results <- readxl::read_xlsx(file.path(directory, "results/all_results.xlsx"))

#Subset to TCGA cases that were included in COSMIC analyses 
tcga_sbs_96 <- read.delim("/Users/kayleighlawson-michod/Library/CloudStorage/OneDrive-UniversityofUtah/Doherty/AACES_Tissue_Grant/1.Analyses/Programs/SigProfilerAssignment/SigProfilerAssignment/tcga_cosmic3.3/Assignment_Solution/Activities/Assignment_Solution_Activities.txt")
tcga_white_sbs_96 <- tcga_sbs_96 %>% filter(., Samples %in% tcga_white_cases_under_79_ids)
tcga_hgsc_white_wes_cosmic_ids <- c(unique(as.character(tcga_white_sbs_96$Samples)))
tcga_final_cosmic <- all_results %>% filter(Study == "TCGA-W" & sample_id %in% tcga_hgsc_white_wes_cosmic_ids)
tcga_final_cosmic <- tcga_final_cosmic %>% mutate(Study = recode(Study, 'TCGA-W' = 'TCGA-W Cosmic')) 
table_1_data <- rbind(all_results, tcga_final_cosmic)

table_1_data_final <- table_1_data %>% 
  mutate_at(vars(dobyear, diagyear), as.numeric, 
            vars(refage_5_year, Race, Study), as.factor) %>% 
  select(-c(dobyear, diagyear, 'Age at diagnosis')) %>% 
  rename('Age at diagnosis'= refage_5_year,
         'Gene Expression Subtype' = ClusterK4_kmeans)

table_1_data_final <- table_1_data_final %>%
  mutate_all(as.character) %>%
  mutate_all(~ ifelse(. == "Missing", NA, .)) %>%
  mutate_all(as.factor)

table_1_data_final <- 
  table_1_data_final %>%
  mutate(`Age at diagnosis` = factor(`Age at diagnosis`, levels = c("35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "≥ 70")), 
         `Year of diagnosis` = factor(`Year of diagnosis`, levels = c("< 2000", "2000-2004", "2005-2009", "2010-2015")), 
         `Year of birth` = factor(`Year of birth`, levels = c("< 1940", "1940-1949", "1950-1959", "≥ 1960")), 
         Study = factor(Study, levels = c("Schildkraut-B", "TCGA-W", "TCGA-W Cosmic")), 
         `Gene Expression Subtype` = factor(`Gene Expression Subtype`, levels = c("C1.MES", "C5.PRO", "C2.IMM", "C4.DIF")))

table_1 <- table_1_data_final %>%
  tbl_summary(by = "Study",
             missing = "no", 
             statistic = list(all_categorical() ~ "{n} ({p}%)"),
             digits = list(all_categorical() ~ c(0, 0)),
             type = all_dichotomous() ~ "categorical",
             include=c("Age at diagnosis", "Year of diagnosis", "Stage", "Ethnicity", "Year of birth", "Neoadjuvant treatment", 'Gene Expression Subtype')) %>%
  modify_header(all_stat_cols() ~ "**{level}** , N = ({n})") %>%
#  modify_spanning_header(all_stat_cols() ~ "**Study Population**") %>%
  modify_table_styling(all_stat_cols(), footnote = "N (%)
                       missing values are excluded") %>%
  bold_labels() %>%
  modify_header(label = "") 

table_1_gt <- table_1 %>%
  as_gt() %>%
  # gt::cols_align(align = "left") %>%
  # gt::cols_width(c(2, 3) ~ px(20)) %>%
  gt::tab_options(table.font.names = "Arial", table.font.size = 12, table.font.color = "black", 
                  table_body.hlines.style = "none",
                  data_row.padding = px(1.5)) %>%
  gt::tab_footnote(
    footnote = "FIGO: International Federation of Gynecology and Obstetrics")

gt::gtsave(table_1_gt, file.path(directory, "results/tables/table_1.png"))
gtsummary::as_hux_xlsx(table_1, file.path(directory, "results/tables/table_1.xlsx"))
