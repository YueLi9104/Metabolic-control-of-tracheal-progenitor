## set up a work directory
setwd("E:/L3vs0hTr5/L3no1vs0hTr5no3 flybaseID")

data <- read.csv("L3no1vs0hTr5no3_Geneid for heatmap.csv")
data
datatrix <- as.matrix(data[, 2:7])
rownames(datatrix) <- data$Geneid
head(datatrix)

grouping_file <- data.frame(row.names <- colnames(datatrix),
  Condition = c(
    "L3",
    "L3",
    "L3",
    "APF0h",
    "APF0h",
    "APF0h"
  )
)

## install Deseq2##
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("DESeq2")
## it will ask if wants update, put 'a' to update all##
a
library(DESeq2)

## option: prefilter low count genes. two ways to do this filter 
## 1:delete genes which have total counts less than N (e.g 10) ##
keepdata <- datatrix[rowSums(datatrix) >= 21, ]
keepdata1 <- as.matrix(keepdata)
## 2:delete genes which have count less than N in more than Y samples##
keepdata2 <- keepdata1[rowSums(keepdata1 < 1) < 4, ]

## the data needs in matrix format rather than table
## so convert the table to matrix##

dataset2 <- as.matrix(keepdata2)

## make a data matrix for Deseq2, countdata=the matrix we just made, colData=sample information file, design normally=~Condition##

dmatrix2 <- DESeqDataSetFromMatrix(
  countData = dataset2, colData = grouping_file,
  design = ~Condition
)

## analyse with DEseq##
analyseddata2 <- DESeq(dmatrix2)
res <- results(analyseddata2, contrast = c("Condition", "APF0h", "L3"))

## to make more graphs, need transform data, three ways to transform counts,I usedvstrans##

ntdtrans <- normTransform(analyseddata2)
assay(ntdtrans)

selected <- which((results(analyseddata2)$padj < 0.05))

## make heatmaps##
library("pheatmap")


test <- assay(ntdtrans)[selected, ]

p <- pheatmap(test,
  fontsize = 10,
  scale = "row",
  border_color = NA,
  na_col = "grey",
  cluster_rows = T, cluster_cols = T,
  show_rownames = F, show_colnames = T,
  treeheight_row = 30, treeheight_col = 10,
  cellheight = , cellwidth = 20,
  cutree_row = , cutree_col = ,
  display_numbers = F, legend = T
)
p

# Adjust the order of clusters
callback <- function(hc, mat) {
  sv <- svd(t(mat))$v[, 1]
  dend <- reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

pheatmap(test,
  clustering_callback = callback,
  fontsize = 10,
  scale = "row",
  border_color = NA,
  na_col = "grey",
  cluster_rows = T, cluster_cols = T,
  show_rownames = F, show_colnames = T,
  treeheight_row = 30, treeheight_col = 20,
  cellheight = , cellwidth = 20,
  cutree_row = , cutree_col = ,
  display_numbers = F, legend = T
)