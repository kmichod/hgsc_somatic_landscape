#HGSC somatic landscape paper figure 1
#Author: Katherine Lawson-Michod

library(maftools)
library(grid)
library(tidyverse)
library(glue)
library(ggpubr)
library(ggtext)

#change directory to location of somatic_landscape folder
directory <- "~/hgsc_somatic_landscape_paper" #set working directory
source(file.path(directory, "code/2.utils.R"))

#Load all data
somatic_epi_merge <- readxl::read_xlsx(file.path(directory, "/results/all_results.xlsx"))
# Rename to merge as clinical data with maf tools
tissue_black_somatic_mutations <- tissue_black_somatic_mutations %>% 
  select(-c(Tumor_Sample_Barcode)) %>% rename("Tumor_Sample_Barcode" = "SUID")
tcga_white_somatic_mutations <- tcga_white_somatic_mutations %>% 
  select(-c(Tumor_Sample_Barcode)) %>% rename("Tumor_Sample_Barcode" = "case_submitter_id")

#Define lists of SUIDs for each gene expression subtype
unique(somatic_epi_merge$ClusterK4_kmeans)
tissue_mes_ids <- somatic_epi_merge %>% filter(Study == "Schildkraut-B", ClusterK4_kmeans == "C1.MES") %>% pull(sample_id)
mes_tissue_black_somatic_mutations <- tissue_black_somatic_mutations %>% filter(Tumor_Sample_Barcode %in% tissue_mes_ids)
tissue_pro_ids <- somatic_epi_merge %>% filter(Study == "Schildkraut-B", ClusterK4_kmeans == "C5.PRO") %>% pull(sample_id)
pro_tissue_black_somatic_mutations <- tissue_black_somatic_mutations %>% filter(Tumor_Sample_Barcode %in% tissue_pro_ids)
tissue_immuno_ids <- somatic_epi_merge %>% filter(Study == "Schildkraut-B", ClusterK4_kmeans == "C2.IMM") %>% pull(sample_id)
immuno_tissue_black_somatic_mutations <- tissue_black_somatic_mutations %>% filter(Tumor_Sample_Barcode %in% tissue_immuno_ids)
tissue_diff_ids <- somatic_epi_merge %>% filter(Study == "Schildkraut-B", ClusterK4_kmeans == "C4.DIF") %>% pull(sample_id)
diff_tissue_black_somatic_mutations <- tissue_black_somatic_mutations %>% filter(Tumor_Sample_Barcode %in% tissue_diff_ids)

tcga_mes_ids <- somatic_epi_merge %>% filter(Study == "TCGA-W", ClusterK4_kmeans == "C1.MES") %>% pull(sample_id)
mes_tcga_white_somatic_mutations <- tcga_white_somatic_mutations %>% filter(Tumor_Sample_Barcode %in% tcga_mes_ids)
tcga_pro_ids <- somatic_epi_merge %>% filter(Study == "TCGA-W", ClusterK4_kmeans == "C5.PRO") %>% pull(sample_id)
pro_tcga_white_somatic_mutations <- tcga_white_somatic_mutations %>% filter(Tumor_Sample_Barcode %in% tcga_pro_ids)
tcga_immuno_ids <- somatic_epi_merge %>% filter(Study == "TCGA-W", ClusterK4_kmeans == "C2.IMM") %>% pull(sample_id)
immuno_tcga_white_somatic_mutations <- tcga_white_somatic_mutations %>% filter(Tumor_Sample_Barcode %in% tcga_immuno_ids)
tcga_diff_ids <- somatic_epi_merge %>% filter(Study == "TCGA-W", ClusterK4_kmeans == "C4.DIF") %>% pull(sample_id)
diff_tcga_white_somatic_mutations <- tcga_white_somatic_mutations %>% filter(Tumor_Sample_Barcode %in% tcga_diff_ids)

#Calculate gene frequency and define lists of the highest frequency genes in tcga and shildkraut
gene_frequency <- function(variants_df, id_column, gene_column)
{
  study <-  gsub("_maf","", deparse(substitute(variants_df)))
  variants_df <- variants_df %>% filter()
  count_df <- variants_df %>% 
    distinct({{id_column}}, {{gene_column}}) %>% 
    group_by({{gene_column}}) %>% tally()
  if (study == 'tcga_white') {count_df <- count_df %>% mutate(sample_size = 272)}
  if (study == 'tissue_black') {count_df <- count_df %>% mutate(sample_size = 191)}
  count_df <- count_df %>%
    mutate(perc = (n/sample_size)*100, 0) %>%
    arrange(desc(perc)) %>% 
    rename(., gene={{gene_column}}, !!paste0("n_", study) := n, !!paste0("sample_size_", study) := sample_size, !!paste0("perc_", study) := perc) %>%
    mutate(gene = as.factor(gene))
  return(count_df)
}

tissue_black_maf <- tissue_black_somatic_tcga_mutations_maf@data
tcga_white_maf <- tcga_white_somatic_mutations_maf@data
tissue_black_maf_freq <- gene_frequency(tissue_black_maf, Tumor_Sample_Barcode, Hugo_Symbol)
tissue_black_maf_freq <- tissue_black_maf_freq %>% filter(., !gene %in% c("CNTNAP5", "HUWE1", "NEB", "PAK3", "SPDYE5", "ZNF729", "ALG1", "LILRA5", "LILRA1")) #filter genes not covered by TCGA 
tcga_white_maf_freq <- gene_frequency(tcga_white_maf, Tumor_Sample_Barcode, Hugo_Symbol)

tissue_black_maf_freq_top <- tissue_black_maf_freq %>% mutate(perc_tissue_black = round(perc_tissue_black, 0)) %>% filter(perc_tissue_black > 3)
tissue_black_maf_freq_top_list <- unique(tissue_black_maf_freq_top$gene)
tcga_white_maf_freq_top <- tcga_white_maf_freq %>% mutate(perc_tcga_white = round(perc_tcga_white, 0)) %>% filter(perc_tcga_white > 3)
tcga_white_maf_freq_top_list <- unique(tcga_white_maf_freq_top$gene)

#Set colors for plotting
vc_cols = c("#A6CEE3", "#1F78B4", "#000000", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#5D3FD3", "#808080", "#B2DF8A", "#FB9A99", "#CAB2D6", "#6A3D9A", "#E6AB02", "#A6761D", "#666666", "#1B9E77")
names(vc_cols) = c("Silent", "Frame shift del", "Multi hit", "Missense", "In frame ins", "Nonsense", "In frame del", "Splice site", "Frame shift ins", "Splice region", "Nonstop", "RNA", "5'UTR", "5'Flank", "3'Flank", "3'UTR", "Intron", "Translation start site")

#Define function to prepare data frame for plotting 
to_counts_for_merge <- function(variants_df, id_column, full_df_input)
{
  study <-  gsub("_somatic_mutations","", deparse(substitute(variants_df)))
  one_mut_per_person <- variants_df %>% 
    distinct({{id_column}}, Hugo_Symbol, Variant_Classification)
  count_df <- one_mut_per_person %>% 
    filter(!Variant_Classification %in% c("Intron")) %>%
    group_by(Hugo_Symbol, Variant_Classification) %>% tally()
  if (study == 'tcga_white') {count_df <- count_df %>% mutate(sample_size = 272)}
  if (study == 'mes_tcga_white') {count_df <- count_df %>% mutate(sample_size = 71)}
  if (study == 'pro_tcga_white') {count_df <- count_df %>% mutate(sample_size = 42)}
  if (study == 'immuno_tcga_white') {count_df <- count_df %>% mutate(sample_size = 48)}
  if (study == 'diff_tcga_white') {count_df <- count_df %>% mutate(sample_size = 91)}
  if (study == 'tissue_black' | study == 'tissue_black_tcga_filtered') {count_df <- count_df %>% mutate(sample_size = 191)}
  if (study == 'mes_tissue_black') {count_df <- count_df %>% mutate(sample_size = 45)}
  if (study == 'pro_tissue_black') {count_df <- count_df %>% mutate(sample_size = 40)}
  if (study == 'immuno_tissue_black') {count_df <- count_df %>% mutate(sample_size = 56)}
  if (study == 'diff_tissue_black') {count_df <- count_df %>% mutate(sample_size = 2)}
  count_df <- count_df %>%
    mutate(perc = (n/sample_size)*100) %>%
    mutate(Variant_Classification=recode(Variant_Classification, 
                                         `Missense_Mutation`="Missense",
                                         `Nonsense_Mutation`="Nonsense",
                                         `Frame_Shift_Del`="Frame shift del", 
                                         `Splice_Site`="Splice site", 
                                         `In_Frame_Del`="In frame del", 
                                         `Frame_Shift_Ins`="Frame shift ins", 
                                         `Splice_Region`="Splice region", 
                                         `In_Frame_Ins`="In frame ins", 
                                         `Nonstop_Mutation`="Nonstop",
                                         `Translation_Start_Site`="Translation start site"))
  if (study == 'tissue_black' | study == 'tcga_white') {count_df <- count_df %>% 
    rename(., `Variant Classification` = `Variant_Classification`, !!paste0("n_", study) := n, !!paste0("sample_size_", study) := sample_size, !!paste0("perc_", study) := perc)}
  if (study == 'mes_tissue_black' | study == 'pro_tissue_black' | study == 'immuno_tissue_black' | study == 'diff_tissue_black') {count_df <- count_df %>% 
    rename(., `Variant Classification` = `Variant_Classification`, n_tissue_black := n, sample_size_tissue_black := sample_size, perc_tissue_black := perc)}
  if (study == 'mes_tcga_white' | study == 'pro_tcga_white' | study == 'immuno_tcga_white' | study == 'diff_tcga_white') {count_df <- count_df %>% 
    rename(., `Variant Classification` = `Variant_Classification`, n_tcga_white := n, sample_size_tcga_white := sample_size, perc_tcga_white := perc)}
  
  one_mut_per_person_full <- full_df_input %>% 
    distinct({{id_column}}, Hugo_Symbol, Variant_Classification)
  
  if (study == 'tissue_black' | study == 'mes_tissue_black' | study == 'pro_tissue_black' | study == 'immuno_tissue_black' | study == 'diff_tissue_black') {
    full_df <- one_mut_per_person_full %>% 
      filter(!Variant_Classification %in% c("Intron")) %>%
      mutate(Variant_Classification=recode(Variant_Classification, 
                                           `Missense_Mutation`="Missense",
                                           `Nonsense_Mutation`="Nonsense",
                                           `Frame_Shift_Del`="Frame shift del", 
                                           `Splice_Site`="Splice site", 
                                           `In_Frame_Del`="In frame del", 
                                           `Frame_Shift_Ins`="Frame shift ins", 
                                           `Splice_Region`="Splice region", 
                                           `In_Frame_Ins`="In frame ins", 
                                           `Nonstop_Mutation`="Nonstop",
                                           `Translation_Start_Site`="Translation start site")) %>%
      rename(`Variant Classification` = Variant_Classification) %>% 
      group_by(Hugo_Symbol, `Variant Classification`) %>% 
      summarise(n_tissue_black = 0, sample_size_tissue_black = 191, perc_tissue_black = 0.000000001) %>% 
      mutate_at(vars(n_tissue_black), as.integer)
    count_df_replace <-  anti_join(full_df, count_df, by = c("Hugo_Symbol", "Variant Classification"))
    count_df_final <- bind_rows(count_df, count_df_replace)} #Add back in signatures that have a frequency of 0 in TCGA for plotting
  if (study == 'tcga_white' | study == 'mes_tcga_white' | study == 'pro_tcga_white' | study == 'immuno_tcga_white' | study == 'diff_tcga_white') {
    full_df <- one_mut_per_person_full %>% 
      filter(!Variant_Classification %in% c("Intron")) %>%
      mutate(Variant_Classification=recode(Variant_Classification, 
                                           `Missense_Mutation`="Missense",
                                           `Nonsense_Mutation`="Nonsense",
                                           `Frame_Shift_Del`="Frame shift del", 
                                           `Splice_Site`="Splice site", 
                                           `In_Frame_Del`="In frame del", 
                                           `Frame_Shift_Ins`="Frame shift ins", 
                                           `Splice_Region`="Splice region", 
                                           `In_Frame_Ins`="In frame ins", 
                                           `Nonstop_Mutation`="Nonstop",
                                           `Translation_Start_Site`="Translation start site")) %>%
      rename(`Variant Classification` = Variant_Classification) %>% 
      group_by(Hugo_Symbol, `Variant Classification`) %>% 
      summarise(n_tcga_white = 0, sample_size_tcga_white = 191, perc_tcga_white = 0.000000001) %>% 
      mutate_at(vars(n_tcga_white), as.integer)
    count_df_replace <-  anti_join(full_df, count_df, by = c("Hugo_Symbol", "Variant Classification"))
    count_df_final <- bind_rows(count_df, count_df_replace)} #Add back in signatures that have a frequency of 0 in TCGA for plotting
  return(count_df_final)
}


#Apply function 
TCGA_merge <- to_counts_for_merge(tcga_white_somatic_mutations, Tumor_Sample_Barcode, tcga_white_somatic_mutations)
Schildkraut_merge <- to_counts_for_merge(tissue_black_somatic_mutations, Tumor_Sample_Barcode, tissue_black_somatic_mutations)
tissue_black_somatic_mutations_mes_merge <- to_counts_for_merge(mes_tissue_black_somatic_mutations, Tumor_Sample_Barcode, tissue_black_somatic_mutations)
tissue_black_somatic_mutations_pro_merge <- to_counts_for_merge(pro_tissue_black_somatic_mutations, Tumor_Sample_Barcode, tissue_black_somatic_mutations)
tissue_black_somatic_mutations_immuno_merge <- to_counts_for_merge(immuno_tissue_black_somatic_mutations, Tumor_Sample_Barcode, tissue_black_somatic_mutations)
tissue_black_somatic_mutations_diff <- to_counts_for_merge(diff_tissue_black_somatic_mutations, Tumor_Sample_Barcode, tissue_black_somatic_mutations)
tcga_white_somatic_mutations_mes_merge <- to_counts_for_merge(mes_tcga_white_somatic_mutations, Tumor_Sample_Barcode, tcga_white_somatic_mutations)
tcga_white_somatic_mutations_pro_merge <- to_counts_for_merge(pro_tcga_white_somatic_mutations, Tumor_Sample_Barcode, tcga_white_somatic_mutations)
tcga_white_somatic_mutations_immuno_merge <- to_counts_for_merge(immuno_tcga_white_somatic_mutations, Tumor_Sample_Barcode, tcga_white_somatic_mutations)
tcga_white_somatic_mutations_diff <- to_counts_for_merge(diff_tcga_white_somatic_mutations, Tumor_Sample_Barcode, tcga_white_somatic_mutations)

#Define function for plotting 
plot_merge_2 <- function(df_1, df_2, gene_list, y_min, y_max)
{
  df_1_filtered <- df_1 %>% filter(., Hugo_Symbol %in% {{gene_list}})
  df_2_filtered <- df_2 %>% filter(., Hugo_Symbol %in% {{gene_list}})
  df_merge <- merge(df_2_filtered, df_1_filtered, all.x=TRUE, all.y=TRUE)
  ggplot(df_merge, aes(x=x) ) +
    geom_col( aes(x = Hugo_Symbol, y = -perc_tissue_black, fill = `Variant Classification`)) +
    geom_col( aes(x = Hugo_Symbol, y = perc_tcga_white, fill = `Variant Classification`)) +
    scale_x_discrete(limits = rev({{gene_list}})) + 
    scale_fill_manual(values = vc_cols) +
    theme_bw() +
    ylab("") +
    xlab("") + 
    geom_hline(yintercept = 0, color = "black") + 
    theme(axis.text.x = element_text(size = 10, color = "black", hjust = 1), 
       #   axis.text.y = element_text(size = 10, face = "italic", color = "black", hjust = 1),
          axis.title.x = element_text(size = 10, color = "black"),
          plot.title = element_text(size = 10, color = "black", face = "bold"),
       #  legend.position = c(.25, .3), 
          legend.text=element_text(size=10), 
          legend.title = element_text(size=10, face = "bold"), 
          legend.key.size = unit(.25, "cm"),
          legend.background=element_rect(fill = alpha("white", 0)), ) +
    scale_y_continuous(limits = c({{y_min}}, {{y_max}}), labels = c(paste({{y_min}}*-1), paste({{y_min}}/-2), "0", paste({{y_max}}/2), paste({{y_max}})), expand = c(0,0)) +
    coord_flip()}

###Figure 2
tissue_black_maf_freq_top_list <- as.character(tissue_black_maf_freq_top_list)
tcga_white_maf_freq_top_list <- as.character(tcga_white_maf_freq_top_list)
freq_list <- append(tissue_black_maf_freq_top_list, tcga_white_maf_freq_top_list)
freq_list_1 <- append(freq_list, tcga_genes)
freq_list_2 <- append(freq_list_1, tissue_genes)
freq_list_3_tp <- unique(freq_list_2)
freq_list_3 <- freq_list_3_tp[-1]
tissue_genes_exclude_tp53 <- c("BRCA1", "KRAS", "GPR107", "GIMAP4", "ZNF681", "RB1", "ANKRD30A")
gene_lists <- list(tcga_genes, tissue_genes, tissue_black_maf_freq_top_list, tcga_white_maf_freq_top_list)
freq_list <- unique(freq_list)
freq_list_exclude_tp53 <- freq_list[-1]

# Function to find which lists each gene is in
find_gene_lists <- function(gene_lists) {
  all_genes <- unique(unlist(gene_lists))
  gene_in_lists <- lapply(all_genes, function(gene) {
    lists_containing_gene <- sapply(gene_lists, function(lst) gene %in% lst)
    which(lists_containing_gene)
  })
  names(gene_in_lists) <- all_genes
  return(gene_in_lists)
}

# Call the function
gene_lists_result <- find_gene_lists(gene_lists)

# Print the result
print(gene_lists_result)

all <- c("TP53", "BRCA1")  #000000
tcga_both_tissue_freq <- c("CSMD3", "FAT3") #555555
tcga_and_tissue_smg_only <- c("RB1") #808080
tissue_and_tcga_freq_only <- c("USH2A", "LRP2", "APOB") #aaaaaa

tcga_both_only <- c("NF1") #b31529
tcga_smg_only <- c("BRCA2", "CDK12", "GABRA6") #d75f4c
tcga_freq_only <- c("HMCN1", "DNAH3", "AHNAK", "HYDIN", "RYR2", "DST") #f6a482

tissue_smg_only <- c("KRAS", "GPR107", "GIMAP4", "ZNF681", "ANKRD30A") #1065ab
tissue_freq_only <- c("MACF1", "RYR1", "WDFY4", "XIRP2", "DNAH7", "MGA", "OBSCN", "PIEZO2", "RYR3", "SYNE1", "UBR4", "TACC2") #3a93c3

#########################################
df_1_filtered <- TCGA_merge %>% filter(., Hugo_Symbol %in% freq_list_3 & `Variant Classification` %in% names(vc_cols)) 
df_2_filtered <- Schildkraut_merge %>% filter(., Hugo_Symbol %in% freq_list_3 & `Variant Classification` %in% names(vc_cols))
df_merge <- left_join(df_2_filtered, df_1_filtered)
df_merge <- df_merge %>%
  mutate(
    Hugo_Symbol_Color = ifelse(Hugo_Symbol %in% all, 
                               glue::glue("<strong><i style='color:#000000;font-size:10px'>{Hugo_Symbol}</i></strong>"),
                               ifelse(Hugo_Symbol %in% tcga_both_tissue_freq, glue::glue("<strong><i style='color:#555555;font-size:10px'>{Hugo_Symbol}</i></strong>"),
                                      ifelse(Hugo_Symbol %in% tcga_and_tissue_smg_only, glue::glue("<strong><i style='color:#808080;font-size:10px'>{Hugo_Symbol}</i></strong>"),
                                             ifelse(Hugo_Symbol %in% tissue_and_tcga_freq_only, glue::glue("<strong><i style='color:#aaaaaa;font-size:10px'>{Hugo_Symbol}</i></strong>"),
                                                    ifelse(Hugo_Symbol %in% tcga_both_only, glue::glue("<strong><i style='color:#990013;font-size:10px'>{Hugo_Symbol}</i></strong>"),
                                                           ifelse(Hugo_Symbol %in% tcga_smg_only, glue::glue("<strong><i style='color:#d75f4c;font-size:10px'>{Hugo_Symbol}</i></strong>"),
                                                                  ifelse(Hugo_Symbol %in% tcga_freq_only, glue::glue("<strong><i style='color:#f6a482;font-size:10px'>{Hugo_Symbol}</i></strong>"),
                                                                         ifelse(Hugo_Symbol %in% tissue_smg_only, glue::glue("<strong><i style='color:#1065ab;font-size:10px'>{Hugo_Symbol}</i></strong>"),
                                                                                glue::glue("<strong><i style='color:#3a93c3;font-size:10px'>{Hugo_Symbol}</i></strong>"))))))))))
df_perc_total <- df_merge %>% 
  group_by(Hugo_Symbol_Color, Hugo_Symbol) %>% 
  summarise(perc_tissue_black_total = sum(perc_tissue_black), 
            perc_tcga_white_total = sum(perc_tcga_white), 
            mean_sum = sum(perc_tissue_black_total, perc_tcga_white_total)) %>%
  arrange(desc(perc_tissue_black_total), desc(mean_sum)) %>%
  mutate(Hugo_Symbol_Color = factor(Hugo_Symbol_Color, levels = unique(Hugo_Symbol_Color)), 
         Hugo_Symbol = factor(Hugo_Symbol_Color, levels = unique(Hugo_Symbol_Color)))

Hugo_Symbol_Levels <- levels(df_perc_total$Hugo_Symbol_Color)
Hugo_Symbol_ONLY_Levels <- levels(df_perc_total$Hugo_Symbol)

#################################################

plot_merge_3 <- function(df_1, df_2, gene_list, y_min, y_max, input_tissue, input_tcga)
{
  df_1_filtered <- df_1 %>% filter(., Hugo_Symbol %in% {{gene_list}} & `Variant Classification` %in% names(vc_cols)) 
  df_2_filtered <- df_2 %>% filter(., Hugo_Symbol %in% {{gene_list}} & `Variant Classification` %in% names(vc_cols))
  df_merge <- left_join(df_2_filtered, df_1_filtered)
  df_merge <- df_merge %>%
    mutate(Hugo_Symbol_Color = glue::glue("<strong><i style='color:#000000;font-size:10px'>{Hugo_Symbol}</i></strong>"))
  df_perc_total <- df_merge %>% 
    group_by(Hugo_Symbol_Color) %>% 
    summarise(perc_tissue_black_total = sum(perc_tissue_black), 
              perc_tcga_white_total = sum(perc_tcga_white), 
              mean_sum = sum(perc_tissue_black_total, perc_tcga_white_total)) %>%
    arrange(desc(perc_tissue_black_total), desc(mean_sum)) %>%
    mutate(Hugo_Symbol_Color = factor(Hugo_Symbol_Color, levels = unique(Hugo_Symbol_Color)))
    
  Hugo_Symbol_Levels <- levels(df_perc_total$Hugo_Symbol_Color)
  
  ggplot(df_merge, aes(x=x) ) +
    geom_col( aes(x = Hugo_Symbol_Color, y = -{{input_tissue}}, fill = `Variant Classification`)) +
    geom_col( aes(x = Hugo_Symbol_Color, y = {{input_tcga}}, fill = `Variant Classification`)) +
    scale_x_discrete(limits = rev(Hugo_Symbol_Levels)) +
    scale_fill_manual(values = vc_cols) +
    theme_bw() +
    ylab("") +
    xlab("") + 
    geom_hline(yintercept = 0, color = "black") + 
    theme(axis.text.x = element_text(size = 10, color = "black", hjust = 1),
          axis.text.y = element_markdown(hjust = 1),
          axis.title.x = element_text(size = 10, color = "black", face = "bold"),
          plot.title = element_text(size = 10, color = "black", face = "bold"),
          legend.position = "none", 
          legend.text=element_text(size=10), 
          legend.title = element_text(size=10, face = "bold"), 
          legend.key.size = unit(.25, "cm"),
          legend.background=element_rect(fill = alpha("white", 0))) +
    scale_y_continuous(limits = c({{y_min}}, {{y_max}}), labels = c(paste({{y_min}}*-1), paste({{y_min}}/-2), "0", paste({{y_max}}/2), paste({{y_max}})), expand = c(0,0)) +
    coord_flip()}


freq_genes_plot <- plot_merge_3(TCGA_merge, Schildkraut_merge, freq_list_3, -20, 20, perc_tissue_black, perc_tcga_white) + ylab("Percent With Mutation") 
freq_genes_plot_mes <- plot_merge_3(tcga_white_somatic_mutations_mes_merge, tissue_black_somatic_mutations_mes_merge, freq_list_3, -20, 20, perc_tissue_black, perc_tcga_white) + theme(axis.text.y = element_blank()) + ylab("Percent With Mutation")
freq_genes_plot_pro <- plot_merge_3(tcga_white_somatic_mutations_pro_merge, tissue_black_somatic_mutations_pro_merge, freq_list_3, -20, 20, perc_tissue_black, perc_tcga_white) + theme(axis.text.y = element_blank()) +ylab("Percent With Mutation")
freq_genes_plot_immuno <- plot_merge_3(tcga_white_somatic_mutations_immuno_merge, tissue_black_somatic_mutations_immuno_merge, freq_list_3, -20, 20, perc_tissue_black, perc_tcga_white) + theme(axis.text.y = element_blank()) +ylab("Percent With Mutation")
freq_genes_plot_diff <- plot_merge_3(tcga_white_somatic_mutations_diff, tissue_black_somatic_mutations_diff, freq_list_3, -20, 20, perc_tissue_black, perc_tcga_white) + theme(axis.text.y = element_blank()) +ylab("Percent With Mutation")

tcga_white_somatic_mutations_mes_merge <- tcga_white_somatic_mutations_mes_merge %>% mutate(gene_expression = "MES")
tcga_white_somatic_mutations_pro_merge <- tcga_white_somatic_mutations_pro_merge %>% mutate(gene_expression = "PRO")
tcga_white_somatic_mutations_immuno_merge <- tcga_white_somatic_mutations_immuno_merge %>% mutate(gene_expression = "IMM")
tcga <- rbind(tcga_white_somatic_mutations_mes_merge, tcga_white_somatic_mutations_pro_merge, tcga_white_somatic_mutations_immuno_merge)
tcga <- tcga %>% filter(Hugo_Symbol %in% freq_list_3) %>% group_by(gene_expression, Hugo_Symbol) %>% summarise(sum(perc_tcga_white)) 
tcga <- tcga %>% arrange((Hugo_Symbol))
writexl::write_xlsx(tcga, file.path(directory, "results/Figure3_tcga.xlsx"))

tissue_black_somatic_mutations_mes_merge <- tissue_black_somatic_mutations_mes_merge %>% mutate(gene_expression = "MES")
tissue_black_somatic_mutations_pro_merge <- tissue_black_somatic_mutations_pro_merge %>% mutate(gene_expression = "PRO")
tissue_black_somatic_mutations_immuno_merge <- tissue_black_somatic_mutations_immuno_merge %>% mutate(gene_expression = "IMM")
tissue <- rbind(tissue_black_somatic_mutations_mes_merge, tissue_black_somatic_mutations_pro_merge, tissue_black_somatic_mutations_immuno_merge)
tissue <- tissue %>% filter(Hugo_Symbol %in% freq_list_3) %>% group_by(gene_expression, Hugo_Symbol) %>% summarise(sum(perc_tissue_black)) 
tissue <- tissue %>% arrange((Hugo_Symbol))
writexl::write_xlsx(tissue, file.path(directory, "results/Figure3_tissue.xlsx"))


tcga_white_somatic_mutations_mes_merge %>% filter(Hugo_Symbol %in% freq_list_3)
freq_genes_plot <- freq_genes_plot +
  ggtitle("         Schildkraut-B               TCGA-W") + 
  theme(plot.margin = margin(0.5, 0, 0, 0, "cm"))

freq_genes_plot_mes <- freq_genes_plot_mes + 
  ggtitle("C1.MES") + 
  theme(plot.margin = margin(0.5, 0, 0, 0, "cm")) + 
  theme(axis.title.x = element_blank())
freq_genes_plot_pro <- freq_genes_plot_pro +  
  ggtitle("C5.PRO") + 
  theme(plot.margin = margin(0.5, 0, 0, 0, "cm")) + 
  theme(axis.text.y = element_blank(), axis.title.x = element_blank())
freq_genes_plot_immuno <- freq_genes_plot_immuno +
  ggtitle("C2.IMM") + 
  theme(plot.margin = margin(0.5, 0, 0, 0, "cm")) + 
  theme(axis.text.y = element_blank(), axis.title.x = element_blank())
freq_genes_plot_diff <- freq_genes_plot_diff +
  ggtitle("C4.DIF") + 
  theme(plot.margin = margin(0.5, 0, 0, 0, "cm")) + 
  theme(axis.text.y = element_blank(), axis.title.x = element_blank())
Figure3_A <- ggarrange(freq_genes_plot) 
Figure3_A <- annotate_figure(Figure3_A, top = textGrob("", gp = gpar(cex = 1.3, fontsize = 8)))
Figure3_B <- ggarrange(freq_genes_plot_mes, freq_genes_plot_pro, freq_genes_plot_immuno, freq_genes_plot_diff, ncol = 4, widths = c(1.5, 1.5, 1.5, 1.5))
Figure3_B <- annotate_figure(Figure3_B, bottom = textGrob(expression(bold("          Percent With Mutation")), gp = gpar(cex = 1.3, fontsize = 8)), top = textGrob("Schildkraut-B is shown on left and TCGA-W is shown on the right for all graphs", gp = gpar(cex = 1.3, fontsize = 8)))
Figure3 <- ggarrange(Figure3_A, Figure3_B, ncol = 2, labels = c("A", "B"), widths = c(1, 2))

tiff(file.path(directory, "/results/figures/figure3.tiff"), units="in", width=11, height=6, res=300)
Figure3
dev.off()

