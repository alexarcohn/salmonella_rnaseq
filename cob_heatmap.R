# cob_heatmap.R
# Alexa R. Cohn
# October 23, 2019
# Last revised January 21, 2020

cobeutpdu_spi1 <- read.delim("~/cobeutpdu_spi1_avg.txt", row.names = "Gene")

cobeutpdu_spi1_matrix <- as.matrix(cobeutpdu_spi1)


library(gplots)
mycol <- colorpanel(1000,"turquoise","white","lightpink1")
heatmap.2(cobeutpdu_spi1_matrix, scale = "row", Rowv = NULL, Colv = NULL,
          dendrogram = "none", 
          col = mycol,
          labRow = rownames(cobeutpdu_spi1),
          labCol = colnames(cobeutpdu_spi1),
          key = TRUE, keysize = 1.0,
          trace = "none", density.info  = "none")
