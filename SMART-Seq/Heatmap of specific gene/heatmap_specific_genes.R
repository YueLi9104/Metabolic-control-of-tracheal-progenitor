## set up a work directory##
setwd("~")

data <- read.csv("L3vs0hr_Geneid for heatmap.csv")
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

library(DESeq2)

## option: prefilter low count genes. two ways to do this filter ##
## 1:delete genes which have total counts less than N (e.g 10) ##
keepdata <- datatrix[rowSums(datatrix) >= 21, ]
keepdata1 <- as.matrix(keepdata)
## 2:delete genes which have count less than N in more than Y samples##
keepdata2 <- keepdata1[rowSums(keepdata1 < 1) < 4, ]

## the data needs in matrix format rather than table##
## so convert the table to matrix##

dataset2 <- as.matrix(keepdata2)

## make a data matrix for Deseq2, countdata=the matrix we just made, colData=sample information file, design normally=~Condition##

dmatrix2 <- DESeqDataSetFromMatrix(
  countData = dataset2, colData = grouping_file,
  design = ~Condition
)

## analyse with DEseq##
analyseddata2 <- DESeq(dmatrix2)


## to make more graphs, need transform data, three ways to transform counts,I usedvstrans###

ntdtrans <- normTransform(analyseddata2)
p <- assay(ntdtrans)

## make heatmaps##
library("pheatmap")

candidate_genes <- c(
  "par-6",
  "sav",
  "dally",
  "upd3",
  "Patj",
  "lft",
  "Rassf",
  "Act42A",
  "Mer",
  "p53",
  "hth",
  "Act57B",
  "kibra",
  "ed",
  "Ldh",
  "Mipp1",
  "Ald1",
  "Pglym78",
  "Tpi",
  "Gapdh1",
  "AcCoAS",
  "Pgk",
  "Eno"
)


selected <- which(rownames(ntdtrans) %in% candidate_genes)

test <- assay(ntdtrans)[selected, ]

p <- pheatmap(test,
  fontsize = 10,
  scale = "row",
  border_color = NA,
  na_col = "grey",
  cluster_rows = T, cluster_cols = T,
  show_rownames = T, show_colnames = T,
  treeheight_row = 30, treeheight_col = 30,
  cellheight = 15, cellwidth = 10,
  cutree_row = , cutree_col = ,
  display_numbers = F, legend = T
)
p

callback <- function(hc, mat) {
  sv <- svd(t(mat))$v[, 1]
  dend <- reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

p1 <- pheatmap(test,
  clustering_callback = callback,
  fontsize = 10,
  scale = "row",
  border_color = NA,
  na_col = "grey",
  cluster_rows = T, cluster_cols = T,
  show_rownames = F, show_colnames = T,
  treeheight_row = 30, treeheight_col = 20,
  cellheight = 15, cellwidth = 10,
  cutree_row = , cutree_col = ,
  display_numbers = F, legend = T
)



gn <- rownames(test)[p1$tree_row[["order"]]]
gn
sn <- colnames(test)[p1$tree_col[["order"]]]
sn
new_test <- test[gn, sn]


pheatmap(new_test,
  fontsize = 10,
  scale = "row",
  border_color = NA,
  na_col = "grey",
  cluster_rows = F, cluster_cols = F,
  show_rownames = T, show_colnames = T,
  treeheight_row = 30, treeheight_col = 30,
  cellheight = 10, cellwidth = 10,
  cutree_row = , cutree_col = ,
  display_numbers = F, legend = T
)