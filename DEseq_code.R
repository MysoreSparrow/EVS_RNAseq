# Author : Keshav Prasad Gubbi

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
suppressPackageStartupMessages(library("glmpca"))
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

## Input the Count Matrix

countsmatrix <- read.csv("~/Documents/Eva_RNAseq/featurecounts_Eva_totalRNAseq.csv",
  stringsAsFactors = TRUE
)

### Clean up the Count Matrix
rownames(countsmatrix) <- countsmatrix[, 2] # converting first column of gene ID into rownames, to be used for sanity check later

# It is IMPORTANT to keep the names of the genes in the rownames
countsmatrix <- subset(countsmatrix, select = -c(X, EnsemblID)) # dropping the X column

## Display the column names
colnames(countsmatrix)

### Annotating and Exporting ENSEMBL ID into Gene Symbols

# Adding genes annotated from ENSEMBL ID to Gene symbols and ENTREZ Id to countsmatrix table.
# Will be keeping the symbols and entrez columsn to be added later into results table as it is for later use

cm_row <- rownames(countsmatrix)
head(cm_row)

# Mapping the ENSEMBL ID to Symbol and ENTREZ ID
symbols <- mapIds(
  org.Mm.eg.db,
  keys = cm_row,
  column = c("SYMBOL"),
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Create symbols column
symbols <- symbols[!is.na(symbols)]
symbols <- symbols[match(rownames(countsmatrix), names(symbols))]

# Creating a new column called genename and putting in the symbols and entrez columns into count matrix
countsmatrix$genename <- symbols

# Removing all rows with NA values for genenames, so that those rows are filtered out.
countsmatrix <- unique(countsmatrix[rowSums(is.na(countsmatrix)) == 0, ]) # Apply rowSums & is.na

nrow(countsmatrix)

# Moving the ENSEMBL ID from rownames into separate column for itself.
countsmatrix <- tibble::rownames_to_column(countsmatrix, "E_ID")
# Removing the duplicated genes so that then these genes can be made into rownames for countsmatrix
countsmatrix <- distinct(countsmatrix[!duplicated(countsmatrix$genename), ])

# Now make the ganename column into rownames of count matrix
rownames(countsmatrix) <- countsmatrix[, "genename"]

# Dropping the column E_ID, genenames so that only numeric values are present in it as an input of DESEq Object.
countsmatrix <- subset(countsmatrix, select = -c(genename, E_ID)) #

# Changing countsmatrix into Matrix of numeric values so that only numeric values are present in it as an input of DESEq Object.
countsmatrix <- as.matrix(countsmatrix)
class(countsmatrix) <- "numeric"


### The Count Matrix is:
head(countsmatrix, 10)


## Develop ColData by importing metadata

# Read the csv file and change the column name. the samples.csv is a list of sample names, ie, the names of bam files.
sample_ID <- read.csv("~/Documents/Eva_RNAseq/Samples.csv")

coldata <- as.data.frame(sample_ID)
colnames(coldata) <- c("Sample_Name", "MouseType", "condition") # change name of one of the columns
rownames(coldata) <- coldata[,1] # move the sample names into row names
# The metadata can be found in a df called coldata!
head(coldata)

## the elements from Sample_Name from coldata must the the colnames of countsmatrix
colnames(countsmatrix) <- coldata$Sample_Name


# Sanity Check for DDS
all(rownames(coldata) %in% colnames(countsmatrix))
ncol(countsmatrix) == nrow(coldata)
dim(countsmatrix)

## Creating the DESeq Data set Object
dds_infected <- DESeqDataSetFromMatrix(
  countData = countsmatrix,
  colData = coldata,
  design = ~condition
)
nrow(dds_infected)


keep <- rowSums(counts(dds_infected)) > 10
dds_infected <- dds_infected[keep, ]
nrow(dds_infected)


## Applying VST transformation
vsd <- vst(dds_infected, blind = FALSE)
# head(assay(vsd), 3)
colData(vsd)
vsd_coldata <- colData(vsd)
dds_infected <- estimateSizeFactors(dds_infected)

### Euclidean Distance between samples


sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$Sample_Name
colnames(sampleDistMatrix) <- vsd$Sample_Name
colors <- colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(255)
(EuclideanDistanceHeatmap <- pheatmap(sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  main = "Sample-to-Sample Euclidean Distance of Samples",
  col = colors
))

### Poisson Distance between Samples
poisd <- PoissonDistance(t(counts(dds_infected))) # raw counts or unnormalised data
samplePoisDistMatrix <- as.matrix(poisd$dd)
rownames(samplePoisDistMatrix) <- dds_infected$Sample_Name
colnames(samplePoisDistMatrix) <- dds_infected$Sample_Name
colors <- colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(255)
(poisson_dist_plot <- pheatmap(samplePoisDistMatrix,
                               clustering_distance_rows = poisd$dd,
                               clustering_distance_cols = poisd$dd,
                               main = "Sample-to-Sample Poisson Distance of Samples",
                               col = colors
))

# PCA Plot

### Functions for Plot aethetics and saving PCA Plots
color_values <- c(
  "red", "red", "red", "red", "black", "black", "red", "red", "red",
  "red", "red", "red", "red", "red", "red", "blue", "red", "red", "red", "blue"
)
# the basic set of common aesthetic settings for PCA plots,
theme.my.own <- list(
  theme_bw(),
  geom_point(size = 3),
  coord_fixed(),
  scale_y_continuous(
    breaks = seq(-20, 20, 5),
    sec.axis = sec_axis(~ . * 1,
                        labels = NULL,
                        breaks = NULL
    )
  ),
  scale_x_continuous(
    breaks = seq(-20, 20, 5),
    sec.axis = sec_axis(~ . * 1,
                        labels = NULL,
                        breaks = NULL
    )
  ),
  theme_classic(),
  geom_hline(yintercept = 0, color = "gray", size = 1),
  geom_vline(xintercept = 0, color = "gray", size = 1),
  theme(
    text = element_text(size = 15),
    axis.text = element_text(size = 15),
    legend.position = "bottom",
    aspect.ratio = 1
  ),
  # geom_text(size = 4, hjust = 0, vjust = 0)
  geom_text_repel(size = 5, min.segment.length = 0.5)
)


## Calculating all PCA Values


plotPCA_local <- function(object,
                          intgroup = "condition",
                          ntop = 500,
                          returnData = TRUE,
                          nPC = 4) {
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
  intgroup.df <-
    as.data.frame(colData(object)[, intgroup, drop = FALSE])
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

### Function for calculating percentvar

percentvar_calculation <- function(pcaData_variable) {
  # function to calculate percentvar for different variables
  percentvar_variable <- round(100 * attr(pcaData_variable, "percentVar"), digits = 3)
  return(percentvar_variable)
}
savingFunction <- function(plotname, metadatacolumn) {
  # Function to save the PCA plots
  ggsave(
    filename =
      glue("~/Documents/Eva_RNAseq/PCAplot_{metadatacolumn}.png"),
    plot = plotname,
    dpi = 300,
    width = 10,
    height = 10,
    units = "in"
  )
}

pcaData_infected <- plotPCA_local(vsd, intgroup = c("condition", "Sample_Name"), returnData = T)
pcaData_infected
percentVar_infected <- percentvar_calculation(pcaData_infected)
percentVar_infected

(PCAplot_vst <- ggplot(
  pcaData_infected,
  aes(
    x = PC1,
    y = PC2,
    color = Sample_Name,
    label = Sample_Name
  )
) +
    xlab(paste0("PC1: ", percentVar_infected[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar_infected[2], "% variance")) +
    ggtitle("PCA") +
    scale_colour_manual(values = color_values) +
    theme.my.own)
savingFunction(PCAplot_vst, "condition")


# calculate the variance for top 500 gene
rv <- rowVars(assay(vsd))
ntop <- 500
# select the ntop genes by variance
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
df1 <- t(assay(vsd)[select, ])

res.pca <- PCA(df1, graph = FALSE, scale.unit = FALSE)
summary.PCA(res.pca)

# Visualize eigenvalues/variances
fviz_screeplot(res.pca, addlabels = TRUE)

library("factoextra")
eig.val <- get_eigenvalue(res.pca)
eig.val

# var <- get_pca_var(res.pca)
# fviz_pca_var(res.pca, repel = TRUE)

# library("corrplot")
# corrplot(var$cos2, is.corr=T)


## Genes + PCA Biplots

fviz_pca_biplot(res.pca,
                repel = TRUE,
                gradient.cols = c("pink", "blue", "yellow", "green", "red", "black")
)

heat.colors <- brewer.pal(6, "RdYlBu")
fviz_pca_var(res.pca,
             col.var = "contrib", repel = TRUE,
             gradient.cols = c("Gray", "blue", "yellow", "orange", "green", "red", "black"),
)

# Contributions of variables to PC2
fviz_pca_contrib(res.pca, choice = "var", axes = 2, top = 25)
#
# # Contributions of variables to PC1
fviz_contrib(res.pca, choice = "var", axes = 1, top = 25)

## Hierarchical Clustering

### applying rlog Transformation


rld <- rlog(dds_infected, blind = FALSE)
head(assay(rld), 3)

### Extract the rlog matrix from the object
rld_mat <- assay(rld) # assay() is function from the "SummarizedExperiment" package that was loaded when you loaded DESeq2
### Compute pairwise correlation values
rld_cor <- cor(rld_mat) ## cor() is a base R function
head(rld_cor) ## check the output of cor(), make note of the rownames and colnames

### Plot heatmap
heat.colors <- brewer.pal(6, "RdYlBu")
(Hclust_plot <- pheatmap(rld_cor,
                         color = heat.colors,
                         main = "Heirarchical Clustering of Samples - Correlation Matrix"
                         # filename = '/home/keshavprasadgubbi/Documents/AlinaRnaSeq/IvC/Hclust_plot.tiff'
))
# Hclust_plot


# DGE Results

### Running the differential expression pipeline
dds1_infected <- DESeq(dds_infected)
# str(dds1)
### Building the results table
res_infected <- results(dds1_infected,
                        contrast = c("condition", "adult", "d7")
)
head(res_infected, 30)
summary(res_infected)

resdf_infected <- as.data.frame(res_infected) # convert the results table to a df

head(resdf_infected, 20)