# Functional Analysis for Natalia Neonatal vs Adult vs Stimulated
# Author: Keshava Prasad Gubbi

# Libraries
suppressPackageStartupMessages(library("clusterProfiler"))
suppressPackageStartupMessages(library("AnnotationDbi"))
suppressPackageStartupMessages(library("readxl"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("org.Mm.eg.db"))
suppressPackageStartupMessages(library("ggrepel"))
suppressPackageStartupMessages(library("corrplot"))
suppressPackageStartupMessages(library("GO.db"))
suppressPackageStartupMessages(library("GOstats"))
suppressPackageStartupMessages(library("enrichplot"))
suppressPackageStartupMessages(library("glue"))
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("rrvgo"))
suppressPackageStartupMessages(library("qpcR") )
here::i_am("natalia_functional_analysis.R")
here()

# Input Data
DE_genes_Neo_R848vsNeo_PBS_qDC2 <- as.data.frame(read_excel("DE genes_Neo R848 vs. Neo PBS_Neo vs. Adult.xlsx", sheet = "Neo R848 vs. Neo PBS_qDC2"))
DE_genes_AdultvsNeo_qDC2 <- as.data.frame(read_excel("DE genes_Neo R848 vs. Neo PBS_Neo vs. Adult.xlsx", sheet = "Adult vs. Neo_qDC2"))

# Function to save generic plots
saveplot <- function(plot, name) {
  # Function to save the plots
  ggsave(filename = file.path(Comparison_path, glue("/{name}_{Comparison}.jpeg")),
         plot = plot, dpi = 300, width = 20, height = 20, units = "in")
}

# RRVGO
rrvgo_function <- function(GOTERM_object){
  simMatrix <- calculateSimMatrix(GOTERM_object$ID,
                                  orgdb = "org.Mm.eg.db",
                                  ont = "BP",
                                  method = "Rel")
  scores <- setNames(-log10(GOTERM_object$qvalue),
                     GOTERM_object$ID)
  reducedTerms <- reduceSimMatrix(simMatrix = simMatrix,
                                  scores,
                                  threshold = 0.7,
                                  orgdb = "org.Mm.eg.db")
  object_list <- list(simMatrix, reducedTerms)
  return(object_list)

}

# Actual Table at play
significant_UPReg_Table <- DE_genes_Neo_R848vsNeo_PBS_qDC2
significant_DOWNReg_Table <- DE_genes_Neo_R848vsNeo_PBS_qDC2
#################Neo_R848vsNeo_PBS_qDC2 ####################################3

# Define the Comparison and also create the folder for saving all plots and results to be saved as per the comparison
Comparison <- "Neo_R848vsNeo_PBS_qDC2"
# Comparison <- "AdultvsNeo_qDC2"
Comparison_path <- file.path(here(), glue("{Comparison}"))
if (!dir.exists(here(Comparison_path))) {dir.create(here(Comparison_path))}

# Impose the respective restrictions on FoldChange and Padj. # Create Up and Down Regulated Tables first, that are also statistically significant.
# # All values of Avg_logFC above 0 is being considered Up Regulated.
significant_UPReg_Table <- filter(significant_UPReg_Table,
                                  (significant_UPReg_Table$avg_logFC > 0) & (significant_UPReg_Table$p_val_adj < 0.05) )

significant_DOWNReg_Table <- filter(significant_DOWNReg_Table,
                                    (significant_DOWNReg_Table$avg_logFC < 0) & (significant_DOWNReg_Table$p_val_adj < 0.05) )

### GO over-representation analysis for UP Regulated Genes
GO_UPRegResults <- enrichGO(gene = significant_UPReg_Table$gene,
                            OrgDb = "org.Mm.eg.db",
                            keyType = "SYMBOL",
                            ont = "BP",
                            pAdjustMethod = "BH",
                            readable = TRUE
)
head(GO_UPRegResults)
write.csv(as.data.frame(GO_UPRegResults),
          file.path(Comparison_path ,glue("GO_UpReg_results_{Comparison}.csv")))
# Plots
# Up Regulated plots
GO_UPReg_Barplot <- plot(barplot(GO_UPRegResults,
                                 showCategory = 25,
                                 font.size = 35,
                                 title = glue(" Up Regulated GO Terms_{Comparison}"),
                                 label_format = 50))
saveplot(GO_UPReg_Barplot, glue("GO_UPReg_Barplot_{Comparison}"))

GO_UPReg_Dotplot <- plot(dotplot(GO_UPRegResults,
                                 showCategory = 25,
                                 font.size = 35,
                                 title = glue(" Up Regulated GO Terms_{Comparison}"),
                                 label_format = 50))
saveplot(GO_UPReg_Dotplot, glue("GO_UPReg_Dotplot_{Comparison}"))


GO_UPReg_Cnetplot <- plot(cnetplot(GO_UPRegResults,
                                   showCategory = 25,
                                   font.size = 25,
                                   label_format = 55,
                                   colorEdge = FALSE
))
saveplot(GO_UPReg_Cnetplot, glue("GO_UPReg_Cnetplot_{Comparison}"))

GO_UPReg_Heatplot <- plot(heatplot(GO_UPRegResults, label_format = 30))
saveplot(GO_UPReg_Heatplot, glue("GO_UPReg_Heatplot_{Comparison}"))

edox2 <- pairwise_termsim(GO_UPRegResults)
GO_UPReg_enrichtreeplot <- plot(treeplot(edox2, color = "p.adjust", fontsize = 8))
saveplot(GO_UPReg_enrichtreeplot, glue("GO_UPReg_enrichtreeplot_{Comparison}"))

termSim_values <- tibble(edox2@termsim)
write.csv(termSim_values, file.path(Comparison_path , glue("termSimValues_UPReg_{Comparison}.csv")))

GO_UPReg_emapplot <- plot(emapplot(edox2, repel = TRUE))
saveplot(GO_UPReg_emapplot, glue("GO_UPReg_emapplot_{Comparison}"))

# GSEA Analysis for Up Regulated List

geneList_UP <- significant_UPReg_Table$avg_logFC
names(geneList_UP) <- as.character(significant_UPReg_Table$gene)
geneList_UP <- sort(geneList_UP, decreasing = TRUE)

gseGO_UP <- gseGO(geneList     = geneList_UP,
              OrgDb = "org.Mm.eg.db",
              keyType = "SYMBOL",
              ont = "BP",
              pAdjustMethod = "BH",
              #readable = TRUE,
              pvalueCutoff = 0.05,
              verbose      = TRUE
              )
gseGO_UP_df <- as.data.frame(gseGO_UP)
write.csv(as.data.frame(gseGO_UP_df),
          file.path(Comparison_path , glue("gseGO_UP_df_{Comparison}.csv")))

if (nrow(gseGO_UP_df > 0)){
  gseGO_UP_gseaplot <- gseaplot2(gseGO_UP, geneSetID = 1:10, base_size = 30,
                                 pvalue_table = TRUE, ES_geom = "line")
  saveplot(gseGO_UP_gseaplot, glue("UP_gseaplot_{Comparison}"))
  break
}

### GO over-representation analysis for DOWN Regulated Genes
GO_DOWNReg_Results <- enrichGO(gene = significant_DOWNReg_Table$gene,
                            OrgDb = "org.Mm.eg.db",
                            keyType = "SYMBOL",
                            ont = "BP",
                            pAdjustMethod = "BH",
                            readable = TRUE
)
write.csv(as.data.frame(GO_DOWNReg_Results),
          file.path(Comparison_path , glue("GO_DOWNReg_Results_{Comparison}.csv")))

# Plots
# Down Regulated plots
GO_DOWNReg_Barplot <- plot(barplot(GO_DOWNReg_Results,
                                 showCategory = 25,
                                 font.size = 35,
                                 title = glue(" Down Regulated GO Terms_{Comparison}"),
                                 label_format = 50))
saveplot(GO_DOWNReg_Barplot, glue("GO_DOWNReg_Barplot_{Comparison}"))

GO_DOWNReg_Dotplot <- plot(dotplot(GO_DOWNReg_Results,
                                 showCategory = 25,
                                 font.size = 35,
                                 title = glue(" Down Regulated GO Terms_{Comparison}"),
                                 label_format = 50))
saveplot(GO_DOWNReg_Dotplot, glue("GO_DOWNReg_Dotplot_{Comparison}"))

GO_DOWNReg_Cnetplot <- plot(cnetplot(GO_DOWNReg_Results,
                                   showCategory = 25,
                                   font.size = 25,
                                   label_format = 55,
                                   colorEdge = FALSE
))
saveplot(GO_DOWNReg_Cnetplot, glue("GO_DOWNReg_Cnetplot_{Comparison}"))

GO_DOWNReg_Heatplot <- plot(heatplot(GO_DOWNReg_Results, label_format = 30))
saveplot(GO_DOWNReg_Heatplot, glue("GO_DOWNReg_Heatplot_{Comparison}"))

edox2_down <- pairwise_termsim(GO_DOWNReg_Results)
GO_DOWNReg_enrichtreeplot <- plot(treeplot(edox2_down, color = "p.adjust", fontsize = 8))
saveplot(GO_DOWNReg_enrichtreeplot, glue("GO_DOWNReg_enrichtreeplot_{Comparison}"))

termSim_values_DOWNReg <- tibble(edox2_down@termsim)
write.csv(termSim_values_DOWNReg, file.path(Comparison_path , glue("termSimValues_DOWNReg_{Comparison}.csv")))

GO_DOWNReg_emapplot <- plot(emapplot(edox2_down, repel = TRUE))
saveplot(GO_DOWNReg_emapplot, glue("GO_DOWNReg_emapplot_{Comparison}"))

# GSEA Analysis for Down Regulated List

geneList_DOWN <- significant_DOWNReg_Table$avg_logFC
names(geneList_DOWN) <- as.character(significant_DOWNReg_Table$gene)
geneList_DOWN <- sort(geneList_DOWN, decreasing = TRUE)

gseGO_DOWN <- gseGO(geneList     = geneList_DOWN,
                  OrgDb = "org.Mm.eg.db",
                  keyType = "SYMBOL",
                  ont = "BP",
                  pAdjustMethod = "BH",
                  #readable = TRUE,
                  pvalueCutoff = 0.05,
                  verbose      = TRUE
)

gseGO_DOWN_df <- as.data.frame(gseGO_DOWN)
gseGO_DOWN_gseaplot <- gseaplot2(gseGO_DOWN, geneSetID = 1:10, base_size = 30,
                                 pvalue_table = TRUE, ES_geom = "line")
saveplot(gseGO_DOWN_gseaplot, glue("DOWN_gseaplot_{Comparison}"))

#################AdultvsNeo_qDC2 ####################################3

# Define the Comparison and also create the folder for saving all plots and results to be saved as per the comparison
Comparison <- "AdultvsNeo_qDC2"
Comparison_path <- file.path(here(), glue("{Comparison}"))
if (!dir.exists(here(Comparison_path))) {dir.create(here(Comparison_path))}

# Change the tables of UpReg and DownReg
significant_UPReg_Table <- DE_genes_AdultvsNeo_qDC2
significant_DOWNReg_Table <- DE_genes_AdultvsNeo_qDC2

# Impose the respective restrictions on FoldChange and Padj. # Create Up and Down Regulated Tables first, that are also statistically significant.
# # All values of Avg_logFC above 0 is being considered Up Regulated.
significant_UPReg_Table <- filter(significant_UPReg_Table, (significant_UPReg_Table$avg_logFC > 0) & (significant_UPReg_Table$p_val_adj < 0.05) )

significant_DOWNReg_Table <- filter(significant_DOWNReg_Table,
                                    (significant_DOWNReg_Table$avg_logFC < 0) & (significant_DOWNReg_Table$p_val_adj < 0.05) )

# Function to save generic plots
saveplot <- function(plot, name) {
  # Function to save the plots
  ggsave(filename = file.path(Comparison_path, glue("/{name}_{Comparison}.png")),
         plot = plot, dpi = 300, width = 25, height = 25, units = "in")
}

### GO over-representation analysis for UP Regulated Genes
GO_UPRegResults <- enrichGO(gene = significant_UPReg_Table$gene,
                            OrgDb = "org.Mm.eg.db",
                            keyType = "SYMBOL",
                            ont = "BP",
                            pAdjustMethod = "BH",
                            readable = TRUE
)
# head(GO_UPRegResults)
write.csv(as.data.frame(GO_UPRegResults),
          file.path(Comparison_path ,glue("GO_UpReg_results_{Comparison}.csv")))

# Plots
# Up Regulated plots
GO_UPReg_Barplot <- plot(barplot(GO_UPRegResults,
                                 showCategory = 25,
                                 font.size = 35,
                                 title = glue(" Up Regulated GO Terms_{Comparison}"),
                                 label_format = 50))
saveplot(GO_UPReg_Barplot, glue("GO_UPReg_Barplot_{Comparison}"))

GO_UPReg_Dotplot <- plot(dotplot(GO_UPRegResults,
                                 showCategory = 25,
                                 font.size = 35,
                                 title = glue(" Up Regulated GO Terms_{Comparison}"),
                                 label_format = 50))
saveplot(GO_UPReg_Dotplot, glue("GO_UPReg_Dotplot_{Comparison}"))

GO_UPReg_Cnetplot <- plot(cnetplot(GO_UPRegResults,
                                   showCategory = 25,
                                   font.size = 25,
                                   label_format = 55,
                                   colorEdge = FALSE
))
saveplot(GO_UPReg_Cnetplot, glue("GO_UPReg_Cnetplot_{Comparison}"))

GO_UPReg_Heatplot <- plot(heatplot(GO_UPRegResults, label_format = 30))
saveplot(GO_UPReg_Heatplot, glue("GO_UPReg_Heatplot_{Comparison}"))

edox2 <- pairwise_termsim(GO_UPRegResults)
GO_UPReg_enrichtreeplot <- plot(treeplot(edox2, color = "p.adjust", fontsize = 8))
saveplot(GO_UPReg_enrichtreeplot, glue("GO_UPReg_enrichtreeplot_{Comparison}"))

termSim_values <- tibble(edox2@termsim)
write.csv(termSim_values, file.path(Comparison_path , glue("termSimValues_UPReg_{Comparison}.csv")))

GO_UPReg_emapplot <- plot(emapplot(edox2, repel = TRUE))
saveplot(GO_UPReg_emapplot, glue("GO_UPReg_emapplot_{Comparison}"))

# GSEA Analysis for Up Regulated List

geneList_UP <- significant_UPReg_Table$avg_logFC
names(geneList_UP) <- as.character(significant_UPReg_Table$gene)
geneList_UP <- sort(geneList_UP, decreasing = TRUE)

gseGO_UP <- gseGO(geneList     = geneList_UP,
                  OrgDb = "org.Mm.eg.db",
                  keyType = "SYMBOL",
                  ont = "BP",
                  pAdjustMethod = "BH",
                  #readable = TRUE,
                  pvalueCutoff = 0.05,
                  verbose      = TRUE
)
gseGO_UP_df <- as.data.frame(gseGO_UP)
if (nrow(gseGO_UP_df > 0)){
  gseGO_UP_gseaplot <- gseaplot2(gseGO_UP, geneSetID = 1:10, base_size = 30,
                                 pvalue_table = TRUE, ES_geom = "line")
  saveplot(gseGO_UP_gseaplot, glue("UP_gseaplot_{Comparison}"))
  break
}

### GO over-representation analysis for DOWN Regulated Genes
GO_DOWNReg_Results <- enrichGO(gene = significant_DOWNReg_Table$gene,
                               OrgDb = "org.Mm.eg.db",
                               keyType = "SYMBOL",
                               ont = "BP",
                               pAdjustMethod = "BH",
                               readable = TRUE
)
write.csv(as.data.frame(GO_DOWNReg_Results),
          file.path(Comparison_path , glue("GO_DOWNReg_Results_{Comparison}.csv")))

# Plots
# Down Regulated plots
GO_DOWNReg_Barplot <- plot(barplot(GO_DOWNReg_Results,
                                   showCategory = 25,
                                   font.size = 35,
                                   title = glue(" Down Regulated GO Terms_{Comparison}"),
                                   label_format = 50))
saveplot(GO_DOWNReg_Barplot, glue("GO_DOWNReg_Barplot_{Comparison}"))

GO_DOWNReg_Dotplot <- plot(dotplot(GO_DOWNReg_Results,
                                   showCategory = 25,
                                   font.size = 35,
                                   title = glue(" Down Regulated GO Terms_{Comparison}"),
                                   label_format = 50))
saveplot(GO_DOWNReg_Dotplot, glue("GO_DOWNReg_Dotplot_{Comparison}"))

GO_DOWNReg_Cnetplot <- plot(cnetplot(GO_DOWNReg_Results,
                                     showCategory = 25,
                                     font.size = 25,
                                     label_format = 55,
                                     colorEdge = FALSE
))
saveplot(GO_DOWNReg_Cnetplot, glue("GO_DOWNReg_Cnetplot_{Comparison}"))

GO_DOWNReg_Heatplot <- plot(heatplot(GO_DOWNReg_Results, label_format = 30))
saveplot(GO_DOWNReg_Heatplot, glue("GO_DOWNReg_Heatplot_{Comparison}"))

edox2_down <- pairwise_termsim(GO_DOWNReg_Results)
GO_DOWNReg_enrichtreeplot <- plot(treeplot(edox2_down, color = "p.adjust", fontsize = 8))
saveplot(GO_DOWNReg_enrichtreeplot, glue("GO_DOWNReg_enrichtreeplot_{Comparison}"))

termSim_values_DOWNReg <- tibble(edox2_down@termsim)
write.csv(termSim_values_DOWNReg, file.path(Comparison_path , glue("termSimValues_DOWNReg_{Comparison}.csv")))

GO_DOWNReg_emapplot <- plot(emapplot(edox2_down, repel = TRUE))
saveplot(GO_DOWNReg_emapplot, glue("GO_DOWNReg_emapplot_{Comparison}"))

# GSEA Analysis for Up Regulated List

geneList_DOWN <- significant_DOWNReg_Table$avg_logFC
names(geneList_DOWN) <- as.character(significant_DOWNReg_Table$gene)
geneList_DOWN <- sort(geneList_DOWN, decreasing = TRUE)

gseGO_DOWN <- gseGO(geneList     = geneList_DOWN,
                    OrgDb = "org.Mm.eg.db",
                    keyType = "SYMBOL",
                    ont = "BP",
                    pAdjustMethod = "BH",
                    #readable = TRUE,
                    pvalueCutoff = 0.05,
                    verbose      = TRUE
)

gseGO_DOWN_df <- as.data.frame(gseGO_DOWN)
gseGO_DOWN_gseaplot <- gseaplot2(gseGO_DOWN, geneSetID = 1:10, base_size = 30,
                                 pvalue_table = TRUE, ES_geom = "line")
saveplot(gseGO_DOWN_gseaplot, glue("DOWN_gseaplot_{Comparison}"))


############################################################
############################################################
# 1. Finding the common Terms between two GO lists
# Input Data
#
GO_DOWNReg_Results_AdultvsNeo_qDC2 <- as.data.frame(read.csv("~/R/Natalia_FunctionalAnalysis/AdultvsNeo_qDC2/GO_DOWNReg_Results_AdultvsNeo_qDC2.csv", stringsAsFactors = FALSE))

GO_DOWNReg_Results_Neo_R848vsNeo_PBS_qDC2 <- as.data.frame(read.csv("~/R/Natalia_FunctionalAnalysis/Neo_R848vsNeo_PBS_qDC2/GO_DOWNReg_Results_Neo_R848vsNeo_PBS_qDC2.csv", stringsAsFactors = FALSE))

#common_terms <- inner_join(df1, df2, by = "ID")
common_terms_adult <- merge(GO_DOWNReg_Results_AdultvsNeo_qDC2, GO_DOWNReg_Results_Neo_R848vsNeo_PBS_qDC2, by.x = "ID", by.y = "ID", stringsAsFactors = FALSE)
common_terms_stimulated <- merge(GO_DOWNReg_Results_Neo_R848vsNeo_PBS_qDC2, GO_DOWNReg_Results_AdultvsNeo_qDC2, by.x = "ID", by.y = "ID", stringsAsFactors = FALSE)

write.csv(as.data.frame(common_terms_adult),
          file.path(Comparison_path ,glue("common_terms_adult.csv")))
write.csv(as.data.frame(common_terms_stimulated),
          file.path(Comparison_path ,glue("common_terms_stimulated.csv")))

#2 . compare clusters via cluster profiler
# Lets first take teh two genelists and put them together into single dataframe . This is a dataframe of unequal length so must be created in an unual way.

genelist_adult <- (DE_genes_AdultvsNeo_qDC2$gene)
genelist_stimulated <- (DE_genes_Neo_R848vsNeo_PBS_qDC2$gene)
max_length <- max(c(length(genelist_adult), length(genelist_stimulated)))    # Find out maximum length
genelist_df <- data.frame(col1 = c(genelist_adult),
                         col2 = c(genelist_stimulated,
                                  rep(NA, max_length - length(genelist_stimulated))))
colnames(genelist_df)   <- c("genelist_adult", "genelist_stimulated")

genelist_compare_adultVsStimulated <- compareCluster(geneClusters = genelist_df,
                          fun = "enrichGO",
                          keyType = "SYMBOL",
                          OrgDb = org.Mm.eg.db)

genelist_compare_adultVsStimulated_Dotplot <- plot(dotplot(genelist_compare_adultVsStimulated,
                                                   showCategory = 25,
                                                   font.size = 25,
                                                   title = glue(" genelist_compare_adultVsStimulated"),
                                                   label_format = 50))
saveplot(genelist_compare_adultVsStimulated_Dotplot, "genelist_compare_adultVsStimulated_DotPlot")

genelist_compare_adultVsStimulated_Cnetplot <- plot(cnetplot(genelist_compare_adultVsStimulated,
                                                             showCategory = 25,
                                                             font.size = 25,
                                                             label_format = 55,
                                                             colorEdge = FALSE
))
saveplot(genelist_compare_adultVsStimulated_Cnetplot, "genelist_compare_adultVsStimulated_Cnetplot")

# 3. Simplify GO enrichment results to Apply to multiple lists of GO IDs
#### https://www.bioconductor.org/packages/devel/bioc/vignettes/simplifyEnrichment/inst/doc/simplifyEnrichment.html

# lt <- as.list(GO_DOWNReg_Results_AdultvsNeo_qDC2, GO_DOWNReg_Results_Neo_R848vsNeo_PBS_qDC2)
#
# simplifyGOFromMultipleLists(GO_DOWNReg_Results_AdultvsNeo_qDC2$ID,
#                             padj_cutoff = 0.05,
#                             ont = "BP")
# 4. rrvgo : Calculate similarity matrix and also perform visualisation


rrvgo_GO_stimulated <- rrvgo_function(GO_DOWNReg_Results_Neo_R848vsNeo_PBS_qDC2)

heatmap_plot_stimulated <- heatmapPlot(simMatrix = rrvgo_GO_stimulated[[1]],
                                        reducedTerms = rrvgo_GO_stimulated[[2]],
                                        annotateParent = TRUE,
                                        annotationLabel = "parentTerm",
                                        fontsize = 6)
saveplot(heatmap_plot_stimulated, "heatmap_GO_DOWNReg_Results_Neo_R848vsNeo_PBS_qDC2")
scatterplot_stimulated <- scatterPlot(rrvgo_GO_stimulated[[1]], rrvgo_GO_stimulated[[2]])
saveplot(scatterplot_stimulated, "scatterplot_GO_DOWNReg_Results_Neo_R848vsNeo_PBS_qDC2")


rrvgo_GO_adult <- rrvgo_function(GO_DOWNReg_Results_Neo_R848vsNeo_PBS_qDC2)
heatmap_plot_adult <- heatmapPlot(simMatrix = rrvgo_GO_adult[[1]],
                                  reducedTerms = rrvgo_GO_adult[[2]],
                                  annotateParent = TRUE,
                                  annotationLabel = "parentTerm",
                                  fontsize = 6)
scatterplot_adult <- scatterPlot(rrvgo_GO_stimulated[[1]], rrvgo_GO_stimulated[[2]])
saveplot(heatmap_plot_adult, glue("heatmap_rrvgo_GO_DOWNReg_Results_AdultvsNeo_qDC2"))
saveplot(scatterplot_adult, glue("scatterplot_rrvgo_GO_DOWNReg_Results_AdultvsNeo_qDC2"))

# Creating Compare clusters again for natalia since she asked it to be done separately for Up and Down Regulated
# Need to make two separate df -Up and Down Reg
#
#
significant_UPReg_Stimulated <- filter(DE_genes_Neo_R848vsNeo_PBS_qDC2,
                                       (DE_genes_Neo_R848vsNeo_PBS_qDC2$avg_logFC > 0) &
                                       (DE_genes_Neo_R848vsNeo_PBS_qDC2$p_val_adj < 0.05) )

significant_DOWNReg_Stimulated <- filter(DE_genes_Neo_R848vsNeo_PBS_qDC2,
                                         (DE_genes_Neo_R848vsNeo_PBS_qDC2$avg_logFC < 0) &
                                           (DE_genes_Neo_R848vsNeo_PBS_qDC2$p_val_adj < 0.05) )

significant_UPReg_adult <- filter(DE_genes_AdultvsNeo_qDC2,
                                  (DE_genes_AdultvsNeo_qDC2$avg_logFC > 0) &
                                    (DE_genes_AdultvsNeo_qDC2$p_val_adj < 0.05) )

significant_DOWNReg_adult <- filter(DE_genes_AdultvsNeo_qDC2,
                                         (DE_genes_AdultvsNeo_qDC2$avg_logFC < 0) &
                                           (DE_genes_AdultvsNeo_qDC2$p_val_adj < 0.05) )

# Compare Cluster comparing Up Reg
genelist_UpReg_adult <- (significant_UPReg_adult$gene)
length(genelist_UpReg_adult)
genelist_UpReg_stimulated <- (significant_UPReg_Stimulated$gene)
length(genelist_UpReg_stimulated)
max_length <- max(c(length(genelist_UpReg_adult), length(genelist_UpReg_stimulated)))    # Find out maximum length
print(max_length)
length(genelist_UpReg_adult) <- max_length
length(genelist_UpReg_stimulated) <- max_length
genelist_UpReg_df <- as.data.frame(cbind(genelist_UpReg_adult, genelist_UpReg_stimulated))
colnames(genelist_UpReg_df)   <- c("genelist_adult_UpRegulated", "genelist_stimulated_UpRegulated")

genelist_compare_UpReg_adultVsStimulated <- compareCluster(geneClusters = genelist_UpReg_df,
                                                     fun = "enrichGO",
                                                     keyType = "SYMBOL",
                                                     OrgDb = org.Mm.eg.db)

genelist_compare_UpReg_adultVsStimulated_Dotplot <- plot(dotplot(genelist_compare_UpReg_adultVsStimulated,
                                                           showCategory = 25,
                                                           font.size = 30,
                                                           title = glue(" genelist_compare_adultVsStimulated"),
                                                           label_format = 60))
saveplot(genelist_compare_UpReg_adultVsStimulated_Dotplot, "genelist_compare_UpReg_adultVsStimulated_Dotplot")

genelist_compare_UpReg_adultVsStimulated_Cnetplot <- plot(cnetplot(genelist_compare_UpReg_adultVsStimulated,
                                                                     showCategory = 25,
                                                                     font.size = 25,
                                                                     label_format = 55,
                                                                     colorEdge = FALSE
))
saveplot(genelist_compare_UpReg_adultVsStimulated_Cnetplot, "genelist_compare_UpReg_adultVsStimulated_Cnetplot")

# Compare Cluster comparing Down Reg
genelist_DOWNReg_adult <- (significant_DOWNReg_adult$gene)
length(genelist_DOWNReg_adult)
genelist_DOWNReg_stimulated <- (significant_DOWNReg_Stimulated$gene)
length(genelist_DOWNReg_stimulated)
max_length <- max(c(length(genelist_DOWNReg_adult), length(genelist_DOWNReg_stimulated)))    # Find out maximum length
print(max_length)
length(genelist_DOWNReg_adult) <- max_length
length(genelist_DOWNReg_stimulated) <- max_length
genelist_DownReg_df <- as.data.frame(cbind(genelist_DOWNReg_adult, genelist_DOWNReg_stimulated))
colnames(genelist_DownReg_df)   <- c("genelist_adult_DownRegulated", "genelist_stimulated_DownRegulated")

genelist_compare_DownReg_adultVsStimulated <- compareCluster(geneClusters = genelist_DownReg_df,
                                                           fun = "enrichGO",
                                                           keyType = "SYMBOL",
                                                           OrgDb = org.Mm.eg.db)

genelist_compare_DownReg_adultVsStimulated_Dotplot <- plot(dotplot(genelist_compare_DownReg_adultVsStimulated,
                                                                 showCategory = 25,
                                                                 font.size = 25,
                                                                 title = glue(" genelist_compare_adultVsStimulated"),
                                                                 label_format = 60))
saveplot(genelist_compare_DownReg_adultVsStimulated_Dotplot, "genelist_compare_DownReg_adultVsStimulated_Dotplot")

genelist_compare_DownReg_adultVsStimulated_Cnetplot <- plot(cnetplot(genelist_compare_DownReg_adultVsStimulated,
                                                             showCategory = 25,
                                                             font.size = 25,
                                                             label_format = 55,
                                                             colorEdge = FALSE
))
saveplot(genelist_compare_DownReg_adultVsStimulated_Cnetplot, "genelist_compare_DownReg_adultVsStimulated_Cnetplot")
