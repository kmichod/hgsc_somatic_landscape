#HGSC somatic landscape paper figure 1
#Author: Katheirne Lawson-Michod
#Date: July 6th, 2023
#Last updated: October 3rd, 2023

library(maftools)

directory <- "/Users/kayleighlawson-michod/Library/CloudStorage/OneDrive-UniversityofUtah/Doherty/AACES_Tissue_Grant/1.Writing/GitHub"
source(file.path(directory, "code/figure_2.R"))
tissue_clin <- somatic_epi_merge %>% 
  filter(Study == "Black Schildkraut") %>% 
  rename("Tumor_Sample_Barcode" = "sample_id", 
         "Gene Expression Subtype" = "ClusterK4_kmeans")
  
tcga_clin <- somatic_epi_merge %>% 
  filter(Study == "White TCGA") %>% 
  rename("Tumor_Sample_Barcode" = "sample_id",
         "Gene Expression Subtype" = "ClusterK4_kmeans")

tissue_black_somatic_tcga_mutations_maf <- read.maf("/Users/kayleighlawson-michod/Library/CloudStorage/OneDrive-UniversityofUtah/Doherty/AACES_Tissue_Grant/1.Analyses/Data/Somatic/MAF/hg38/high_confidence/tissue_black_tcga_filtered_somatic_mutations.maf", clinicalData = tissue_clin)
tcga_white_somatic_mutations_maf <- read.maf("/Users/kayleighlawson-michod/Library/CloudStorage/OneDrive-UniversityofUtah/Doherty/AACES_Tissue_Grant/1.Analyses/Data/TCGA/somatic/tcga_white_somatic_mutations_filtered.maf", clinicalData = tcga_clin)


## Supplemental Figure 4

tiff(file.path(directory, "/results/figures/s4_tissue_balck_tiv.tiff"), units="in", width=3, height=2, res=300)
maftools::titv(tissue_black_somatic_tcga_mutations_maf)
dev.off()

tiff(file.path(directory, "/results/figures/s4_tcga_white_tiv.tiff"), units="in", width=8, height=7, res=300)
maftools::titv(tcga_white_somatic_mutations_maf)
dev.off()

## Supplemental Figure 5

tissue_exclude_tp53 <- tissue_genes[-1]
tissue_genes_plot_excludeTP53 <- plot_merge_2(TCGA_merge, Schildkraut_merge, tissue_exclude_tp53, -5, 5)

tissue_genes_plot_excludeTP53 <- tissue_genes_plot_excludeTP53 + 
  theme(legend.position = "bottom", legend.justification = "left") + 
  ggtitle("                       Schildkraut-B                                         TCGA-W    ") + guides(fill = guide_legend(title.position = "top")) + ylab("Population (%) With Mutation") 

tiff(file.path(directory, "/results/figures/s5_smg_tissue.tiff"), units="in", width=6, height=10, res=300)
ggarrange(tissue_genes_plot_excludeTP53, tissue_genes_plot_excludeTP53, nrow = 2, labels = c("A", "B"))
dev.off()

## Supplemental Figure 6

#Define colors for the gene expression subtypes in the oncoplots 
subtype_cols = c("#89CFF0", "#EE4B2B", "#00FF00", "#FF00FF", "#36454F")
names(subtype_cols) = c("Mesenchymal", "Proliferative", "Immunoreactive", "Differentiated", "Not included in RNASeq")
subtype_cols = list('Gene_Expression_Subtype' = subtype_cols)


tiff(file.path(directory, "/results/figures/s6_tissue_landscape_overview_191.tiff"), units="in", width=10, height=8, res=300)
maftools::oncoplot(
  tissue_black_somatic_tcga_mutations_maf, 
  genes = freq_list_3_tp, 
  draw_titv = FALSE, 
  gene_mar = 7, 
 # colors = vc_cols, 
  titleText = "Schildkraut-B", 
  clinicalFeatures = 'Gene_Expression_Subtype',
  annotationColor = subtype_cols,
  sortByAnnotation = TRUE, 
  annotationOrder = c("Mesenchymal", "Proliferative", "Immunoreactive", "Differentiated", "NA"))
dev.off()

tiff(file.path(directory, "/results/figures/s6_tcga_landscape_overview.tiff"), units="in", width=10, height=8, res=300)
oncoplot_tcga <- maftools::oncoplot(
  tcga_white_somatic_mutations_maf, 
  draw_titv = FALSE, 
  genes = freq_list_3_tp, 
  gene_mar = 7, 
  #colors = vc_cols, 
  titleText = "TCGA-W", 
  clinicalFeatures = 'Gene_Expression_Subtype',
  annotationColor = subtype_cols,
  sortByAnnotation = TRUE, annotationOrder = c("Mesenchymal", "Proliferative", "Immunoreactive", "Differentiated", "NA"))
dev.off()

tiff("/Users/kayleighlawson-michod/Library/CloudStorage/OneDrive-UniversityofUtah/Doherty/AACES_Tissue_Grant/1.Writing/somatic_landscape/figures/figure_2.tiff", units="in", width=5, height=5, res=300)
tissue_genes_plot_excludeTP53
dev.off()


