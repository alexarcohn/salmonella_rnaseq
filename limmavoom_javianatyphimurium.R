### limmavoom_javianatyphimurium.R
### Alexa R. Cohn
### June 30, 2019
### Last edited July 2, 2019
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
fc <- read.delim("~/Documents/RNA_seq_Salmonella/Javiana vs. Typhimurium/javianatyphimurium_featurecounts_core.txt", stringsAsFactors = FALSE, header = TRUE, row.names = 1)

# convert raw counts to log-CPM values
lcpm <- cpm(fc, log=TRUE)

# create a DGEList using edgeR
dge <- DGEList(counts = fc)

# apply TMM normalization method to counts
dge <- calcNormFactors(dge)
dge$samples$norm.factors

# remove rows that consistently have zero or very low counts
keep <- filterByExpr(dge)
dge <- dge[keep,,keep.lib.sizes=FALSE]
dim(dge)

# derive experiment information from sample names
snames <- colnames(fc)
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
contr <- makeContrasts(serovarJaviana - serovarTyphimurium, levels = colnames(coef(fit)))
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
write.table(top.table, file = "JvT.txt", row.names = F, sep = "\t", quote = F)

# adjust logFC cut off
tfit <- treat(tmp, lfc = 2)
dt <- decideTests(tfit)
summary(dt)

# plot results as a mean-difference plot (displays log-FCs from linear model fit against the log CPM values)
anno <- rownames(fc)
glMDPlot(tfit, coef=1, status=dt, main=colnames(tfit)[1], side.main="Gene", counts=fc, groups=serovar, launch=TRUE, anno=NULL)

# create a heatmap of most differentially expressed genes
JvT <- topTreat(tfit, coef=1, n=Inf)
JvT.topgenes <- rownames(JvT)[1:50]
i <- which(rownames(lcpm) %in% JvT.topgenes)
mycol <- colorpanel(1000,"turquoise","white","lightpink1")
dist <- dist(lcpm[i], method = "euclidean")
c_dist <- dist(t(lcpm), method = "euclidean")
clust <- hclust(dist, method = "average")
c_clust <- hclust(c_dist, method = "average")
heatmap.2(lcpm[i,], scale="row",
          labRow=rownames(lcpm)[i], labCol=serovar, 
          col=mycol, trace="none", density.info="none", 
          margin=c(8,6), lhei=c(2,10), Rowv=as.dendrogram(clust), Colv = as.dendrogram(c_clust))

coolmap(lcpm[i,], col = mycol)

JvT.DE <- read.delim("~/Documents/RNA_seq_Salmonella/JvT_DE.txt")
goana(de=tfit)

genes <- as.integer(top.table$adj.P.Val<0.05)
names(genes) = rownames(top.table)
genes

roary <- read.csv("~/Documents/RNA_seq_Salmonella/gene_presence_absence_48.csv")
typhimurium <- read.delim("~/Documents/RNA_seq_Salmonella/typhimurium/14028s_NCBI.txt")
typhimurium$length <- typhimurium$end - typhimurium$start + 1
lengths <- merge(roary, typhimurium, by = "X14028s_NCBI", all = T)

length_data <- as.integer(lengths$length)
names(length_data) = lengths$Gene
length_data <- length_data[!is.na(length_data)]

GOterms <- read.delim("~/Documents/RNA_seq_Salmonella/typhimurium/gene_GOterms.txt", stringsAsFactors = FALSE, header = TRUE)

pwf <- nullp(genes, bias.data = length_data, plot.fit = TRUE)
