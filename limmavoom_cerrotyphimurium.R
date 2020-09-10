### limmavoom_cerrotyphimurium.R
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
CT_fc <- read.delim("~/Documents/RNA_seq_Salmonella/Cerro vs. Typhimurium/cerrotyphimurium_featurecounts_core.txt", stringsAsFactors = FALSE, header = TRUE, row.names = 1)

# convert raw counts to log-CPM values
CT_lcpm <- cpm(CT_fc, log=TRUE)

# create a DGEList using edgeR
CT_dge <- DGEList(counts = CT_fc)

# apply TMM normalization method to counts
CT_dge <- calcNormFactors(CT_dge)
CT_dge$samples$norm.factors

# remove rows that consistently have zero or very low counts
keep <- filterByExpr(CT_dge)
CT_dge <- CT_dge[keep,,keep.lib.sizes=FALSE]
dim(CT_dge)

# derive experiment information from sample names
snames <- colnames(CT_fc)
snames
serovar <- substr(snames, 1, nchar(snames) - 5)

# plot data in a multidimensional scaling plot
plotMDS(CT_dge, col=as.numeric(serovar))

# specify the model to be fitted
mm <- model.matrix(~0 + serovar)

# adjust for sample-level variation using voom
v <- voom(CT_dge, mm, plot=T)

# fit a linear model with limma
CT_fit <- lmFit(v, mm)
head(coef(CT_fit))

# specify which groups to compare
contr <- makeContrasts(serovarCerro - serovarTyphimurium, levels = colnames(coef(CT_fit)))
contr

# estimate contrasts for each gene
tmp <- contrasts.fit(CT_fit, contr)

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
write.table(top.table, file = "CvT.txt", row.names = F, sep = "\t", quote = F)

# adjust logFC cut off
tfit <- treat(tmp, lfc = 2)
dt <- decideTests(tfit)
summary(dt)

# plot results as a mean-difference plot (displays log-FCs from linear model fit against the log CPM values)
glMDPlot(tfit, coef=1, status=dt, main=colnames(tfit)[1], side.main="Gene", counts=CT_fc, groups=serovar, launch=TRUE)

# create a heatmap of most differentially expressed genes
CvT <- topTreat(tfit, coef=1, n=Inf)
CvT.topgenes <- rownames(CvT)[1:50]
i <- which(rownames(CT_lcpm) %in% CvT.topgenes)
mycol <- colorpanel(1000,"turquoise","white","lightpink1")
dist <- dist(CT_lcpm[i,], method = "euclidean")
c_dist <- dist(t(CT_lcpm[i,]), method = "euclidean")
clust <- hclust(dist, method = "ward.D2")
c_clust <- hclust(c_dist, method = "ward.D2")
pdf(file="heatmap_CvT.pdf", width = 7, height = 10)
heatmap.2(CT_lcpm[i,], scale="row",
          labRow=rownames(CT_lcpm)[i], labCol=serovar, 
          col=mycol, trace="none", density.info="none", 
          margin=c(8,6), lhei=c(2,10), Rowv=as.dendrogram(clust), 
          Colv=as.dendrogram(c_clust))
dev.off()

coolmap(CT_lcpm[i,], col = mycol)

