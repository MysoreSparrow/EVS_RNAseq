# # Author : Keshav Prasad Gubbi
# # Title: "Eva_RNAseq Pipeline"
#
# # Libraries
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("pheatmap"))
suppressPackageStartupMessages(library("PoiClaClu"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("PoiClaClu"))
suppressPackageStartupMessages(library("vsn"))
suppressPackageStartupMessages(library("EnhancedVolcano"))
suppressPackageStartupMessages(library("gplots"))
suppressPackageStartupMessages(library("org.Mm.eg.db"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("genefilter"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("org.Mm.eg.db"))
suppressPackageStartupMessages(library("AnnotationDbi"))
suppressPackageStartupMessages(library("apeglm"))
suppressPackageStartupMessages(library("ComplexHeatmap"))
suppressPackageStartupMessages(library("clusterProfiler"))
suppressPackageStartupMessages(library("ggrepel"))
suppressPackageStartupMessages(library("corrplot"))
suppressPackageStartupMessages(library("GO.db"))
suppressPackageStartupMessages(library("edgeR"))
suppressPackageStartupMessages(library("GOstats"))
suppressPackageStartupMessages(library("pathview"))
suppressPackageStartupMessages(library("gage"))
suppressPackageStartupMessages(library("gageData"))
suppressPackageStartupMessages(library("GOSemSim"))
suppressPackageStartupMessages(library("DOSE"))
suppressPackageStartupMessages(library("enrichplot"))
suppressPackageStartupMessages(library("ggnewscale"))
suppressPackageStartupMessages(library("glue"))
suppressPackageStartupMessages(library("ggupset"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("stats"))
suppressPackageStartupMessages(library("FactoMineR"))
suppressPackageStartupMessages(library("factoextra"))
suppressPackageStartupMessages(library("pcaExplorer"))
suppressPackageStartupMessages(library("here"))

# File Path Declarations
here::i_am("Eva_TotalRNAseq.R")
here()

# Also Create a comparison Variable: That Could be used later for all other comparison titles using a
# glue Variable. Define the Comparison and also create the folder for saving all plots and results to be
# saved as per the comparison

Comparison <- "adult_spfVsd7_spf"
# Comparison <- "d7_GFVsd7_spF"
# Comparison <- "adult_GFVsd7_GF"
# Comparison <- "d7_WTVsd7_spF"
# Comparison <- "d7_WTVsadult_WT"
# Comparison <- "adult_GFVsd7_GF"

# Folder Paths for Different Comparisons
Comparison_path <- file.path(here(), glue("{Comparison}"))
if (!dir.exists(here(Comparison_path))) { dir.create(here(Comparison_path)) } else { print("Dir already exists!") }
paste0(Comparison_path)

# for UninfectedVSInfected is the volcano plot notation: for Deseq2 results,
# the infected is Numerator.and uninfected is Denominator.

# counts Matrix processing
## Input the Count Matrix
countsmatrix <- read.csv(file.path(here(), "featurecounts_Eva_totalRNAseq.csv"))
### Clean up the Count Matrix
rownames(countsmatrix) <- countsmatrix[, 2] # converting first column of gene ID into rownames,
# to be used for sanity check later
nrow(countsmatrix)
# GenderGenes Filtering}
gendergenes_biomart <- data.frame(read.csv("~/R/Eva_TotalRNASeq/20221017_y-chromosomal genes_biomart.txt",
                                           stringsAsFactors = FALSE))
## Remove rows with only empty cells
gendergenes_biomart <- gendergenes_biomart[!apply(gendergenes_biomart == " ", 1, all),]
nrow(gendergenes_biomart)
colnames(gendergenes_biomart)[1] <- "EnsemblID" # Just trying to get the colnames of ensembl column to be same.

# We can use the anti_join() function to return all rows in the first data frame that do not have a matching row
# in the second data frame
countsmatrix <- anti_join(countsmatrix, gendergenes_biomart, by = 'EnsemblID')
## Removal of Gender Genes from ENSEMBL ID column itself
# Filter out the other 5 known gender genes as well from other RNAseq Projects.
countsmatrix <- countsmatrix %>% filter(countsmatrix$EnsemblID != "ENSMUSG00000086503",
                                        countsmatrix$EnsemblID != "ENSMUSG00000097571",
                                        countsmatrix$EnsemblID != "ENSMUSG00000086370",
                                        countsmatrix$EnsemblID != "ENSMUSG00000031329",
                                        countsmatrix$EnsemblID != "ENSMUSG00000030057")
nrow(countsmatrix)

# It is IMPORTANT to keep the names of the genes in the rownames
countsmatrix <- subset(countsmatrix, select = -c(X, EnsemblID)) # dropping the X column
## Display the column names
colnames(countsmatrix)

### Annotating and Exporting ENSEMBL ID into Gene Symbols

# Adding genes annotated from ENSEMBL ID to Gene symbols and ENTREZ Id to countsmatrix table.
# Will be keeping the symbols and entrez columns to be added later into results table as it is for later use
cm_row <- rownames(countsmatrix)
# Mapping the ENSEMBL ID to Symbol and ENTREZ ID
symbols <- mapIds(org.Mm.eg.db, keys = cm_row, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

# Create symbols column
symbols <- symbols[!is.na(symbols)]
symbols <- symbols[match(rownames(countsmatrix), names(symbols))]

# Creating a new column called genename and putting in the symbols and entrez columns into count matrix
countsmatrix$genename <- symbols

# Removing all rows with NA values for genenames, so that those rows are filtered out.
countsmatrix <- unique(countsmatrix[rowSums(is.na(countsmatrix)) == 0, ]) # Apply rowSums & is.na
# Moving the ENSEMBL ID from rownames into separate column for itself.
countsmatrix <- tibble::rownames_to_column(countsmatrix, "E_ID")
# Removing the duplicated genes so that then these genes can be made into rownames for countsmatrix
countsmatrix <- distinct(countsmatrix[!duplicated(countsmatrix$genename), ])
# Now make the ganename column into rownames of count matrix
rownames(countsmatrix) <- countsmatrix[, "genename"]

# Dropping the column E_ID, genenames so that only numeric values are present in it as an input of DESEq Object.
countsmatrix <- subset(countsmatrix, select = -c(genename, E_ID))

# Changing countsmatrix into Matrix of numeric values so that only numeric values are present in it as an input of DESEq Object.
countsmatrix <- as.matrix(countsmatrix)
class(countsmatrix) <- "numeric"

## Develop ColData by importing metadata
# Read the csv file and change the column name. the samples.csv is a list of sample names, ie, the names of bam files.
coldata <- as.data.frame(read.csv(file.path(here(), "Samples_eva.csv")))
colnames(coldata) <- c("Sample_Name", "MouseType", "condition") # change name of one of the columns
rownames(coldata) <- coldata[, 1] # move the sample names into row names
# The metadata can be found in a df called coldata!
## the elements from Sample_Name from coldata must the the colnames of countsmatrix
colnames(countsmatrix) <- coldata$Sample_Name

# # *****************Now to the comparisons*************
### Reduce larger Matrix to smaller one - based on comparison
paste0(Comparison )
switch(Comparison,
       "adult_spfVsd7_spf" = {(coldata <- coldata[c(9, 10, 11, 18, 19, 20), ]) },
       "d7_WTVsd7_spF" = {(coldata <- coldata[c(1, 2, 4, 9, 10, 11), ]) }
)
switch(Comparison,
       "adult_spfVsd7_spf" = {(countsmatrix <- countsmatrix[, c(9, 10, 11, 18, 19, 20)]) },
       "d7_WTVsd7_spF" = {(countsmatrix <- countsmatrix[, c(1, 2, 4, 9, 10, 11)]) }
)
# **********************FUNCTIONS**************************************************************************************
# Function to save generic plots
saveplot <- function(plot, plotname) {
  # Function to save the plots
  extension <- ".jpeg"
  ggsave(filename = file.path(Comparison_path, paste(plotname, glue("_{Comparison}_"), extension),sep=""),
         plot = plot, dpi = 300, width = 10, height = 10, units = "in")
  dev.off()
}
# **********************DESeq Analysis********************************
# Sanity Check for DDS
all(rownames(coldata) %in% colnames(countsmatrix))
ncol(countsmatrix) == nrow(coldata)
dim(countsmatrix)

## Creating the DESeq Data set Object
dds <- DESeqDataSetFromMatrix(
  countData = countsmatrix,
  colData = coldata,
  #design = ~MouseType
  design = ~condition
)
# Further filtering of low count genes
keep <- rowSums(counts(dds)) > 10
dds <- dds[keep, ]
nrow(dds)
## Applying VST transformation
vsd <- vst(dds, blind = FALSE)
# head(assay(vsd), 3)
colData(vsd)
vsd_coldata <- colData(vsd)
dds <- estimateSizeFactors(dds)
### Euclidean Distance between samples
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$Sample_Name
colnames(sampleDistMatrix) <- vsd$Sample_Name
colors <- colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(255)
(EuclideanDistanceHeatmap <- pheatmap(sampleDistMatrix,
                                      clustering_distance_rows = sampleDists,
                                      clustering_distance_cols = sampleDists,
                                      main = glue("Euclidean Distance of Samples: {Comparison}"),
                                      col = colors
))
### Poisson Distance between Samples
poisd <- PoissonDistance(t(counts(dds))) # raw counts or unnormalised data
samplePoisDistMatrix <- as.matrix(poisd$dd)
rownames(samplePoisDistMatrix) <- dds$Sample_Name
colnames(samplePoisDistMatrix) <- dds$Sample_Name
colors <- colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(255)
(poisson_dist_plot <- pheatmap(samplePoisDistMatrix,
                               clustering_distance_rows = poisd$dd,
                               clustering_distance_cols = poisd$dd,
                               main = glue("Poisson Distance of Samples: {Comparison}"),
                               col = colors
))
#  **************************PCA Plot**********************************
# Functions for Plot aethetics and saving PCA Plots
color_values <- c("black", "black", "red", "red", "red", "red", "red", "red", "red", "red", "black", "black")
# The basic set of common aesthetic settings for PCA plots,
theme.my.own <- list(theme_bw(),
                      geom_point(size = 3),
                      coord_fixed(),
                      scale_y_continuous(breaks = seq(-100, 100, 10),
                                         sec.axis = sec_axis(~ . * 1, labels = NULL, breaks = NULL)),
                      scale_x_continuous(breaks = seq(-50, 50, 10),
                                         sec.axis = sec_axis(~ . * 1, labels = NULL, breaks = NULL)),
                      theme_classic(),
                      geom_hline(yintercept = 0, color = "gray", linewidth = 1),
                      geom_vline(xintercept = 0, color = "gray", linewidth = 1),
                      theme(text = element_text(size = 15),
                            axis.text = element_text(size = 15),
                            legend.position = "bottom",
                            aspect.ratio = 1),
                      #geom_text(size = 4, hjust = 0, vjust = 0),
                      geom_text_repel(size = 4, min.segment.length = 0.1)
)
# PCA Plot Calculation
# Calculating all PCA Values
plotPCA_local <- function(object, intgroup = "condition", ntop = 500, returnData = TRUE, nPC = 4) {
  # calculate the variance for each gene
  rv <- rowVars(assay(object))
  ntop <- 500
  # select the ntop genes by variance
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select, ]))
  # summary(pca)
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = ":"))
  } else {
    colData(object)[[intgroup]]
  }
  # assembly the data for the plot
  d <- cbind(
    pca$x[, seq_len(min(nPC, ncol(pca$x))), drop = FALSE],
    data.frame(group = group, intgroup.df, name = colnames(object))
  )
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:nPC]
    # l <- list(pca,d)
    # return(l)
    return(d)
  }
}
## PCA Plot with VST Data
### Function for calculating percentvar for different variables
percentvar_calculation <- function(pcaData_variable) {
  percentvar_variable <- round(100 * attr(pcaData_variable, "percentVar"), digits = 3)
  return(percentvar_variable)
}

pcaData <- plotPCA_local(vsd, intgroup = c("condition", "Sample_Name"), returnData = T)
#pcaData <- plotPCA_local(vsd, intgroup = c("MouseType", "Sample_Name"), returnData = T)
pcaData
percentVar <- percentvar_calculation(pcaData)
percentVar

# PC Plot: PC1 vs PC2
(PCAplot_vst <- ggplot(pcaData, aes(x = PC1, y = PC2, color = Sample_Name, label = Sample_Name)) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle(glue("PCA: {Comparison}")) +
  #scale_colour_manual(values = color_values) +
  theme.my.own)
saveplot(PCAplot_vst, plotname = "PCA_PC1vsPC2")

# PCA Plot : PC3 vs PC4
(PCAplot_pc34 <- ggplot(
  pcaData,
  aes(x = PC3,y = PC4, color = Sample_Name, label = Sample_Name)) +
  xlab(paste0("PC3: ", percentVar[3], "% variance")) +
  ylab(paste0("PC4: ", percentVar[4], "% variance")) +
  ggtitle(glue("PCA: {Comparison}")) +
  scale_colour_manual(values = color_values) +
  theme.my.own)
saveplot(PCAplot_pc34, plotname = "PCA_PC3vsPC4")

# ************************FactoExtra************************
# calculate the variance for top 500 gene

rv <- rowVars(assay(vsd))
ntop <- 500
# select the ntop genes by variance
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
res.pca <- PCA(t(assay(vsd)[select, ]), graph = FALSE, scale.unit = FALSE)
summary.PCA(res.pca)

# Visualize eigenvalues/variances
fviz_screeplot(res.pca, addlabels = TRUE)
eig.val <- get_eigenvalue(res.pca)
var <- get_pca_var(res.pca)

## Genes + PCA Biplots
heat.colors <- brewer.pal(6, "RdYlBu")
## Genes + PCA Biplots
(Genes_Biplot <- fviz_pca_biplot(res.pca, repel = TRUE))
saveplot(Genes_Biplot, "Genes_Biplot")

(Genes_contributions_Biplot <- fviz_pca_var(res.pca, col.var = "contrib", repel = TRUE,
                                            gradient.cols = c("Gray", "blue", "pink", "yellow",
                                                              "orange", "green", "red", "black")))
saveplot(Genes_contributions_Biplot, "Genes_contributions_Biplot")
# Contributions of variables to PC2
(top25_genes_dim2 <- fviz_contrib(res.pca, choice = "var", axes = 2, top = 25))
saveplot(top25_genes_dim2, "top25_genes_dim2")
# # Contributions of variables to PC1
(top25_genes_dim1 <- fviz_contrib(res.pca, choice = "var", axes = 1, top = 25))
saveplot(top25_genes_dim1, "top25_genes_dim1")