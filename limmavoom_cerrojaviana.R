### limmavoom_cerrojaviana.R
### Alexa R. Cohn
### July 3, 2019
### Last edited March 23, 2020
### input file: counts file from rna_seq_featureCounter.R

# install limma+voom and edgeR
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")
BiocManager::install("edgeR")
BiocManager::install("Glimma") 
BiocManager::install("GO.db")
install.packages("gplots")

# load limma+voom, edgeR, glimma, and gplots
library("limma")
library("edgeR")
library("Glimma")
library("gplots")
library("GO.db")

# use cjt matrix from rna_seq_featureCounter.R
CJ_fc <- read.delim("~/Documents/RNA_seq_Salmonella/Cerro vs Javiana/cerrojaviana_featurecounts_core.txt", stringsAsFactors = FALSE, header = TRUE, row.names = 1)

# convert raw counts to log-CPM values
CJ_lcpm <- cpm(CJ_fc, log=TRUE)

# create a DGEList using edgeR
dge <- DGEList(counts = CJ_fc)

# apply TMM normalization method to counts
dge <- calcNormFactors(dge)
dge$samples$norm.factors

# remove rows that consistently have zero or very low counts
keep <- filterByExpr(dge)
dge <- dge[keep,,keep.lib.sizes=FALSE]
dim(dge)

# derive experiment information from sample names
snames <- colnames(CJ_fc)
snames
serovar <- substr(snames, 1, nchar(snames) - 5)

# plot data in a multidimensional scaling plot
plotMDS(dge, col=as.numeric(serovar))

# specify the model to be fitted
mm <- model.matrix(~0 + serovar)

# adjust for sample-level variation using voom
v <- voom(dge, mm, plot=T)

# fit a linear model with limma
fit <- lmFit(v, mm)
head(coef(fit))

# specify which groups to compare
contr <- makeContrasts(serovarCerro - serovarJaviana, levels = colnames(coef(fit)))
contr

# estimate contrasts for each gene
tmp <- contrasts.fit(fit, contr)

# empirical Bayes smoothing of standard errors (shrinks standard errors that are much larger or smaller than those from other genes towards the average standard error)
tmp <- eBayes(tmp)

# determine which genes are most differentially expressed
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)

# determine how many differentially expressed genes there are
length(which(top.table$adj.P.Val < 0.05))

# write top.table to a file
top.table$Gene <- rownames(top.table)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]
write.table(top.table, file = "CvJ.txt", row.names = F, sep = "\t", quote = F)

# adjust logFC cut off
tfit <- treat(tmp, lfc = 2)
dt <- decideTests(tfit)
summary(dt)

# plot results as a mean-difference plot (displays log-FCs from linear model fit against the log CPM values)
glMDPlot(tfit, coef=1, status=dt, main=colnames(tfit)[1], side.main="Gene", counts=CJ_fc, groups=serovar, launch=TRUE)

# create a heatmap of most differentially expressed genes
CvJ <- topTreat(tfit, coef=1, n=Inf)
CvJ.topgenes <- rownames(CvJ)[1:50]
i <- which(rownames(CJ_lcpm) %in% CvJ.topgenes)
mycol <- colorpanel(1000,"turquoise","white","lightpink1")
dist <- dist(CJ_lcpm[i,], method = "euclidean")
c_dist <- dist(t(CJ_lcpm[i,]), method = "euclidean")
clust <- hclust(dist, method = "ward.D2")
c_clust <- hclust(c_dist, method = "ward.D2")
dev.off()
heatmap.2(CJ_lcpm[i,], scale="row",
          labRow=rownames(CJ_lcpm)[i], labCol=serovar, 
          col=mycol, trace="none", density.info="none", 
          margin=c(8,6), lhei=c(2,10), Rowv = as.dendrogram(clust),
          Colv = as.dendrogram(c_clust),
          key = TRUE, keysize = 3)

coolmap(CJ_lcpm[i,], col = mycol)
