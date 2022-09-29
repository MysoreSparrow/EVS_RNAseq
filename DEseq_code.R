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