#HGSC somatic landscape paper figure 1 and 2
#Author: Katherine Lawson-Michod

library(maftools)
library(ggpubr)
library("ggVennDiagram")

directory <- "~/hgsc_somatic_landscape_paper/" #set working directory
source(file.path(directory, "/code/2.utils.R"))

# Figure 1
# Define variables and values in Figure 1A Table
tissue_black_tcga_filtered_somatic_mutations <- tissue_black_somatic_mutations
genes_black_n <- length(unique(tissue_black_tcga_filtered_somatic_mutations$Hugo_Symbol))
coding_black_n <- nrow(tissue_black_tcga_filtered_somatic_mutations)
count_per_gene_black <- tissue_black_tcga_filtered_somatic_mutations %>% group_by(SUID) %>% tally()
min_black <- min(count_per_gene_black$n)
max_black <- max(count_per_gene_black$n)
median_black <- median(count_per_gene_black$n)
coding_per_case_black <- paste0(median_black, " (", min_black, "-", max_black, ")")
tmb_black <- maftools::tmb(tissue_black_somatic_tcga_mutations_maf)
median_tmb_black <- median(tmb_black$total_perMB)

genes_white_n <- length(unique(tcga_white_somatic_mutations$Hugo_Symbol))
coding_white_n <- nrow(tcga_white_somatic_mutations) #16102
count_per_gene_tcga <- tcga_white_somatic_mutations %>% group_by(Tumor_Sample_Barcode) %>% tally()
min_white <- min(count_per_gene_tcga$n)
max_white <- max(count_per_gene_tcga$n)
median_white <- median(count_per_gene_tcga$n)
coding_per_case_white<- paste0(median_white, " (", min_white, "-", max_white, ")")
tmb_white <- maftools::tmb(tcga_white_somatic_mutations_maf)
median_tmb_white <- 0.80

#Make Figure 1A Table
Figure1_TableA <- data.frame(Study = c("Schildkraut", "TCGA"),
                             N = c("191", "272"), 
                             Genes=c(genes_black_n, genes_white_n),
                             "Coding Variants"=c(coding_black_n, coding_white_n),
                             "Coding Variants Per Case (Median, Range)"= c(coding_per_case_black, coding_per_case_white), 
                             "Tumor Mutational Burden (TMB)"= c(median_tmb_black, median_tmb_white))


Figure1_TableA <- Figure1_TableA %>% dplyr::rename(
  "Samples (N)" = "N", 
  "Genes Mutated (N)" = "Genes", 
  "Coding Variants (N)" = Coding.Variants, 
                                            "Median Coding Variants (Range)" = Coding.Variants.Per.Case..Median..Range., 
                                            "TMB" = Tumor.Mutational.Burden..TMB.)
rownames(Figure1_TableA)
Figure1_TableA <- Figure1_TableA %>% select(-c("Study"))
Figure1_TableA <- as_tibble(t(Figure1_TableA))
Figure1_TableA <- Figure1_TableA %>% 
  rename('Schildkraut-B' = V1, 'TCGA-W' = V2)

Figure1A <- ggtexttable(Figure1_TableA, 
                        rows = c("Samples", "Mutated Genes", "Total Coding Variants", "Median Coding Variants Per Sample (Range)", "Tumor Mutation Burden"), 
                        cols = colnames(Figure1_TableA),
                        theme = ttheme("light", base_size = 10, rownames.style = rownames_style(size = 10, face = "bold")))


#Make histograms for panel B and C
gene_frequency <- function(variants_df, id_column, gene_column) #define function for calculating gene frequency
{
  study <-  gsub("_maf","", deparse(substitute(variants_df)))
  sample_size_temp <- variants_df %>% pull(Tumor_Sample_Barcode) %>% unique()
  count_df <- variants_df %>% 
    distinct({{id_column}}, {{gene_column}}) %>% 
    group_by({{gene_column}}) %>% tally()
  if (study == 'tcga_white') {count_df <- count_df %>% mutate(sample_size = 272)}
  if (study == 'tissue_black') {count_df <- count_df %>% mutate(sample_size = 191)}
  count_df <- count_df %>%
    mutate(perc = (n/sample_size)*100) %>%
    arrange(desc(perc)) %>% 
    rename(., gene={{gene_column}}, !!paste0("n_", study) := n, !!paste0("sample_size_", study) := sample_size, !!paste0("perc_", study) := perc) %>%
    mutate(gene = as.factor(gene))
  return(count_df)
}

tissue_black_maf <- tissue_black_somatic_tcga_mutations_maf@data #pull somatic data from maf file
tissue_black_maf <- tissue_black_maf %>% filter(., !Hugo_Symbol %in% c("CNTNAP5", "HUWE1", "NEB", "PAK3", "SPDYE5", "ZNF729", "ALG1", "LILRA5", "LILRA1")) #filter genes not covered by TCGA 
tcga_white_maf <- tcga_white_somatic_mutations_maf@data #pull somatic data from maf file
tissue_black_maf_freq <- gene_frequency(tissue_black_maf, Tumor_Sample_Barcode, Hugo_Symbol)
tcga_white_maf_freq <- gene_frequency(tcga_white_maf, Tumor_Sample_Barcode, Hugo_Symbol)

tissue_black_maf_freq_top <- tissue_black_maf_freq %>% mutate(perc_tissue_black = round(perc_tissue_black, 0)) %>% filter(perc_tissue_black > 3)
tissue_black_maf_freq_top_list <- unique(tissue_black_maf_freq_top$gene)
tcga_white_maf_freq_top <- tcga_white_maf_freq %>% mutate(perc_tcga_white = round(perc_tcga_white, 0)) %>% filter(perc_tcga_white > 3)
tcga_white_maf_freq_top_list <- unique(tcga_white_maf_freq_top$gene)

tissue_unique <- setdiff(tissue_black_maf_freq_top_list, tcga_white_maf_freq_top_list) 
tcga_unique <- setdiff(tcga_white_maf_freq_top_list, tissue_black_maf_freq_top_list)
tcga_white_maf_freq_tissue_genes <- tcga_white_maf_freq %>% filter(gene %in% tissue_unique)
tissue_black_maf_freq_tcga_genes <- tissue_black_maf_freq %>% filter(gene %in% tcga_unique)
tissue_top_comparison <- merge(tissue_black_maf_freq_top, tcga_white_maf_freq_tissue_genes)
tissue_top_comparison <- tissue_top_comparison %>% mutate(diff = perc_tissue_black-perc_tcga_white)
tcga_top_comparison <- merge(tcga_white_maf_freq_top, tissue_black_maf_freq_tcga_genes)
tcga_top_comparison <- tcga_top_comparison %>% mutate(diff = perc_tissue_black-perc_tcga_white)

all_genes <- merge(tissue_black_maf_freq, tcga_white_maf_freq)
all_genes <- all_genes %>% mutate(diff = perc_tissue_black - perc_tcga_white)
tissue_black_maf_freq_stack <- tissue_black_maf_freq %>% rename(perc = perc_tissue_black) %>% mutate(study = "Schildkraut-B") %>% select(gene, perc, study)
tcga_white_maf_freq_stack <- tcga_white_maf_freq %>% rename(perc = perc_tcga_white) %>% mutate(study = "TCGA-W") %>% select(gene, perc, study)
all_genes_stack <- rbind(tissue_black_maf_freq_stack, tcga_white_maf_freq_stack)
all_genes_stack$study <- factor(all_genes_stack$study, levels = c("Schildkraut-B", "TCGA-W"))
all_genes_stack %>% filter(gene != "TP53") %>% pull(perc) %>% max()

density_plot <- 
  all_genes_stack %>% 
  rename(Study = study) %>%
  ggplot(., aes(perc, fill=Study)) + 
  geom_density(position = "dodge", adjust = 3.2, alpha = 0.5) + 
  xlim(0, 7) +
  ylim(0, 2) +
  theme_bw() +
  ylab("Probability Density") +
  scale_fill_manual(values=c("#00BFC4", "#F8766D")) + 
  xlab("Percent With Mutation") +
  theme(axis.text.x = element_text(size = 8, color = "black", hjust = 1), 
#        axis.text.y = element_blank(),
#        axis.ticks.y = element_blank(), 
        axis.title.x = element_text(size = 8, color = "black", face = "bold"), 
        axis.title.y = element_text(size = 8, color = "black", face = "bold"),
        legend.key.size = unit(.25, "cm"),
        legend.background=element_rect(fill = alpha("white", 0)), 
        legend.position = c(.7, .89), 
        legend.text=element_text(size=8, face = "bold"), 
        legend.title = element_blank())

diff_plot <- 
  all_genes %>% 
  mutate(diff_abs = abs(diff)) %>%
  ggplot(aes(diff_abs)) + 
  geom_density(adjust = 2) + 
  xlim(0, 7) + 
  ylim(0, 2) +
  theme_bw() +
  ylab(" ") +
  xlab("Abs Diff in % With Mut") + 
  theme(axis.text.x = element_text(size = 8, color = "black", hjust = 1), 
#        axis.text.y = element_blank(),
#        axis.ticks.y = element_blank(), 
        axis.title.x = element_text(size = 8, color = "black", face = "bold"), 
        axis.title.y = element_text(size = 8, color = "black", face = "bold"),
        legend.text=element_text(size=8), 
        legend.title = element_text(size=8, face = "bold"))

#MutSig plots for panel B
tcga_exclude_tp53 <- tcga_genes[-1]
tissue_exclude_tp53 <- tissue_genes[-1]

#Set colors for plotting
vc_cols = c("#A6CEE3", "#1F78B4", "#000000", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#5D3FD3", "#808080", "#B2DF8A", "#FB9A99", "#CAB2D6", "#6A3D9A", "#E6AB02", "#A6761D", "#666666", "#1B9E77")
names(vc_cols) = c("Silent", "Frame shift del", "Multi hit", "Missense", "In frame ins", "Nonsense", "In frame del", "Splice site", "Frame shift ins", "Splice region", "Nonstop", "RNA", "5'UTR", "5'Flank", "3'Flank", "3'UTR", "Intron", "Translation start site")

to_counts_for_merge <- function(variants_df, id_column)
{
  study <-  gsub("_somatic_mutations","", deparse(substitute(variants_df)))
  one_mut_per_person <- variants_df %>% 
    distinct({{id_column}}, Hugo_Symbol, Variant_Classification)
  count_df <- one_mut_per_person %>% 
    filter(!Variant_Classification %in% c("Intron")) %>%
    group_by(Hugo_Symbol, Variant_Classification) %>% tally()
  if (study == 'tcga_white') {count_df <- count_df %>% mutate(sample_size = 272)}
  if (study == 'tissue_black' | study == 'tissue_black_tcga_filtered') {count_df <- count_df %>% mutate(sample_size = 191)}
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
                                         `Translation_Start_Site`="Translation start site")) %>%
    rename(., `Variant Classification` = `Variant_Classification`, !!paste0("n_", study) := n, !!paste0("sample_size_", study) := sample_size, !!paste0("perc_", study) := perc)
  return(count_df)
}

TCGA_merge <- to_counts_for_merge(tcga_white_somatic_mutations, Tumor_Sample_Barcode)
Schildkraut_merge <- to_counts_for_merge(tissue_black_somatic_mutations, Tumor_Sample_Barcode)

plot_merge_2 <- function(df_1, df_2, gene_list, y_min, y_max)
{
  df_1_filtered <- df_1 %>% filter(., Hugo_Symbol %in% {{gene_list}})
  df_2_filtered <- df_2 %>% filter(., Hugo_Symbol %in% {{gene_list}})
  df_merge <- merge(df_2_filtered, df_1_filtered, all.x = TRUE, all.y = TRUE)
  ggplot(df_merge, aes(x=x) ) +
    geom_col( aes(x = Hugo_Symbol, y = -perc_tissue_black, fill = `Variant Classification`)) +
    geom_col( aes(x = Hugo_Symbol, y = perc_tcga_white, fill = `Variant Classification`)) +
    scale_x_discrete(limits = rev({{gene_list}})) + 
    scale_fill_manual(values = vc_cols) +
    theme_bw() +
    ylab("Percent With Mutation") +
    xlab("") + 
    geom_hline(yintercept = 0, color = "black") + 
    theme(axis.text.x = element_text(size = 10, color = "black", hjust = 1), 
          axis.text.y = element_text(size = 10, face = "bold.italic", color = "black", hjust = 1),
          axis.title.x = element_text(size = 10, color = "black", face = "bold"),
          plot.title = element_text(size = 10, color = "black", face = "bold"),
       #  legend.position = c(.25, .3), 
          legend.text=element_text(size=10), 
          legend.title = element_text(size=10, face = "bold"), 
          legend.key.size = unit(.25, "cm"),
          legend.background=element_rect(fill = alpha("white", 0)), ) +
 #   labs(fill="Variant\nClass") +
    scale_y_continuous(limits = c({{y_min}}, {{y_max}}), labels = c(paste({{y_min}}*-1), paste({{y_min}}/-2), "0", paste({{y_max}}/2), paste({{y_max}})), expand = c(0,0)) +
    coord_flip()}

tcga_genes_plot <- plot_merge_2(TCGA_merge, Schildkraut_merge, tcga_genes, -100, 100)
tcga_genes_plot_excludeTP53 <- plot_merge_2(TCGA_merge, Schildkraut_merge, tcga_exclude_tp53, -10, 10) + theme(legend.position = "none")

tcga_genes_plot <- tcga_genes_plot + theme(legend.position = "bottom", legend.justification = "right") + ggtitle("            Schildkraut-B                  TCGA-W") + guides(fill = guide_legend(title.position = "top"))
tcga_genes_plot_excludeTP53 <- tcga_genes_plot_excludeTP53 + ggtitle("                    Schildkraut-B                                         TCGA-W")

density_plot <- density_plot + theme(plot.margin = margin(0.75, 0, 0, 0, "cm"))
diff_plot <- diff_plot + theme(plot.margin = margin(0.75, 0, 0, 0, "cm"))
tcga_genes_plot <- tcga_genes_plot + theme(plot.margin = margin(0.5, 0, 0, 0, "cm"))
tcga_genes_plot_excludeTP53 <- tcga_genes_plot_excludeTP53 + theme(plot.margin = margin(0.5, 0, 0, 0, "cm"))

Figure1_A_plots <- ggarrange(Figure1A, density_plot, diff_plot, ncol=3, labels = c("A", "B", "C"), widths = c(2.5, 1, 1))
Figure1_B <- ggarrange(tcga_genes_plot, tcga_genes_plot_excludeTP53, ncol=2, labels =c("D", ""), widths = c(1, 1.40))
Figure1 <- ggarrange(Figure1_A_plots, Figure1_B, nrow = 2, heights = c(1, 1.5))

tiff(file.path(directory, "results/figures/figure_1.tiff"), units="in", width=10, height=5, res=300)
Figure1 
dev.off() 


######################################################################FIGURE 2######################################################################

tissue_genes_plot <- plot_merge_2(TCGA_merge, Schildkraut_merge, tissue_genes_exclude_tp53, -5, 5)
tissue_genes_plot <- tissue_genes_plot + theme(legend.position = "bottom", legend.justification = "right") + ggtitle("      Schildkraut-B                      TCGA-W") + guides(fill = guide_legend(title.position = "top"))

tissue_black_maf_freq_top_list_exclude_tp53 <- tissue_black_maf_freq_top_list[-1]
tcga_white_maf_freq_top_list_exclude_tp53 <- tcga_white_maf_freq_top_list[-1]

#define levels order for the frequency only based analyses - schildkrautB
schildkrautB_freq_only_data <- Schildkraut_merge %>% 
  filter(., Hugo_Symbol %in% tissue_black_maf_freq_top_list_exclude_tp53) %>%
  mutate(Hugo_Symbol_Color = glue::glue("<strong><i style='color:#000000;font-size:10px'>{Hugo_Symbol}</i></strong>"))
schildkrautB_freq_only_data_total <- schildkrautB_freq_only_data %>% 
  group_by(Hugo_Symbol_Color) %>% 
  summarise(perc_tissue_black_total = sum(perc_tissue_black)) %>%
  arrange(desc(perc_tissue_black_total)) %>%
  mutate(Hugo_Symbol_Color = factor(Hugo_Symbol_Color, levels = unique(Hugo_Symbol_Color)))
Hugo_Symbol_Levels_SchildkrautB_freq_only <- levels(schildkrautB_freq_only_data_total$Hugo_Symbol_Color)

schildkrautB_freq_only_data$Hugo_Symbol <- 
  factor(schildkrautB_freq_only_data$Hugo_Symbol, levels=Hugo_Symbol_ONLY_Levels)

schildkrautB_freq_only_plot <- schildkrautB_freq_only_data %>%
  ggplot() +
  geom_col( aes(x = Hugo_Symbol_Color, y = perc_tissue_black, fill = `Variant Classification`)) +
  scale_x_discrete(limits = rev(Hugo_Symbol_Levels_SchildkrautB_freq_only)) +
  scale_fill_manual(values = vc_cols) +
  theme_bw() +
  ylab("Percent With Mutation") +
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
  scale_y_continuous(limits = c(0, 10)) +
  coord_flip()

#define levels order for the frequency only based analyses - tcgaW
TCGA_merge_freq_only_data <- TCGA_merge %>% 
  filter(., Hugo_Symbol %in% tcga_white_maf_freq_top_list_exclude_tp53) %>%
  mutate(Hugo_Symbol_Color = glue::glue("<strong><i style='color:#000000;font-size:10px'>{Hugo_Symbol}</i></strong>"))
TCGA_merge_freq_only_data_total <- TCGA_merge_freq_only_data %>% 
  group_by(Hugo_Symbol_Color) %>% 
  summarise(perc_tcga_white_total = sum(perc_tcga_white)) %>%
  arrange(desc(perc_tcga_white_total)) %>%
  mutate(Hugo_Symbol_Color = factor(Hugo_Symbol_Color, levels = unique(Hugo_Symbol_Color)))
Hugo_Symbol_Levels_tcgaW_freq_only <- levels(TCGA_merge_freq_only_data_total$Hugo_Symbol_Color)

tcgaW_freq_only_plot <- TCGA_merge_freq_only_data %>%
  ggplot() +
  geom_col( aes(x = Hugo_Symbol_Color, y = perc_tcga_white, fill = `Variant Classification`)) +
  scale_x_discrete(limits = rev(Hugo_Symbol_Levels_tcgaW_freq_only)) +
  scale_fill_manual(values = vc_cols) +
  theme_bw() +
  ylab("Percent With Mutation") +
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
  scale_y_continuous(limits = c(0, 10)) + 
  coord_flip()


tissue_genes_plot_margins <- tissue_genes_plot + theme(plot.margin = margin(0.75, 0, 0, 0, "cm"))
schildkrautB_freq_only_plot_margins <- schildkrautB_freq_only_plot + theme(plot.margin = margin(0.75, 0, 0, 0, "cm")) + ggtitle("Schildkraut-B")
tcgaW_freq_only_plot_margins <- tcgaW_freq_only_plot + theme(plot.margin = margin(0.75, 0, 0, 0, "cm")) + ggtitle("TCGA-W")

tiff(file.path(directory, "results/figures/figure_2.tiff"), units="in", width=12, height=6, res=300)
ggarrange(tissue_genes_plot_margins, schildkrautB_freq_only_plot_margins, tcgaW_freq_only_plot_margins, ncol=3, labels = c("A", "B", "C"), widths = c(1, 1, 1))
dev.off()
