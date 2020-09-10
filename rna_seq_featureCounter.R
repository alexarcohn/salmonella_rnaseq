### rna_seq_featureCounter.R
### Laura M. Carroll
### April 20, 2019
### input: bam files of mapped reads, saf file output by gff2saf.py script

# uncomment to install Rsubread if necessary
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("Rsubread", version = "3.8")


# read.saf function to read saf file from gff2saf.py
read.saf <- function(saf){
  saf.og <- read.delim(file = saf, header = T, sep = "\t")
  saf <- saf.og[,1:5]
  # returns saf file for featureCounts, as well as the original saf file (saf_full) from gff2saf.py
  # see featureCounts for description of saf file
  # featureCounts: https://www.rdocumentation.org/packages/Rsubread/versions/1.22.2/topics/featureCounts
  return(list(saf = saf, saf_full = saf.og))
}


# load Rsubread
library(Rsubread)

# read Cerro saf file
saf.cerro <- read.saf(saf = "~/Documents/RNA_seq_Salmonella/cerro/FSL-R8-2349_S17_L001_contigs_long.saf")

# run featureCounts for Cerro, supplying sorted bam files and Cerro saf file
# keep strandSpecific to 2 (reverse-stranded RNA-Seq)
# keep isPairedEnd to F (this is single-end data)
# change other options if desired; see featureCounts for details
fc.cerro <- featureCounts(files = c("~/Documents/RNA_seq_Salmonella/cerro/cerro_AC4.sorted.bam", 
                              "~/Documents/RNA_seq_Salmonella/cerro/cerro_AC5.sorted.bam",
                              "~/Documents/RNA_seq_Salmonella/cerro/cerro_AC6.sorted.bam"), 
                    annot.ext = saf.cerro$saf, strandSpecific = 2, isPairedEnd = F)
# stats to see how many reads mapped to features, didn't map, etc.
fc.cerro$stat
# see the top part of the annotation file just because
 head(fc.cerro$annotation)
 write.table(x = data.frame(fc.cerro$annotation[,c("GeneID", "Length")], fc.cerro$counts, stringsAsFactors = FALSE), file = "cerro_featurecounts.txt", append = F, quote = F, sep = "\t", row.names = F)

# read Javiana saf file
saf.javiana <- read.saf(saf = "~/Documents/RNA_seq_Salmonella/javiana/S5_0395_long.saf")

# run featureCounts for Javiana, supplying sorted bam files and Javiana saf file
# keep strandSpecific to 2 (reverse-stranded RNA-Seq)
# keep isPairedEnd to F (this is single-end data)
# change other options if desired; see featureCounts for details
fc.javiana <- featureCounts(files = c("~/Documents/RNA_seq_Salmonella/javiana/javiana_AC1.sorted.bam", 
                                    "~/Documents/RNA_seq_Salmonella/javiana/javiana_AC2.sorted.bam",
                                    "~/Documents/RNA_seq_Salmonella/javiana/javiana_AC3.sorted.bam"), 
                          annot.ext = saf.javiana$saf, strandSpecific = 2, isPairedEnd = F)
# stats to see how many reads mapped to features, didn't map, etc.
fc.javiana$stat
# see the top part of the annotation file just because
head(fc.javiana$annotation)
write.table(x = data.frame(fc.javiana$annotation[,c("GeneID", "Length")], fc.javiana$counts, stringsAsFactors = FALSE), file = "javiana_featurecounts.txt", append = F, quote = F, sep = "\t", row.names = F)

# read Typhimurium saf file
saf.typhimurium <- read.saf(saf = "~/Documents/RNA_seq_Salmonella/typhimurium/14028s_NCBI.saf")

# run featureCounts for Typhimurium, supplying sorted bam files and Typhimurium saf file
# keep strandSpecific to 2 (reverse-stranded RNA-Seq)
# keep isPairedEnd to F (this is single-end data)
# change other options if desired; see featureCounts for details
fc.typhimurium <- featureCounts(files = c("~/Documents/RNA_seq_Salmonella/typhimurium/typhimurium_AC7.sorted.bam", 
                                      "~/Documents/RNA_seq_Salmonella/typhimurium/typhimurium_AC8.sorted.bam",
                                      "~/Documents/RNA_seq_Salmonella/typhimurium/typhimurium_AC9.sorted.bam"), 
                            annot.ext = saf.typhimurium$saf, strandSpecific = 2, isPairedEnd = F)
# stats to see how many reads mapped to features, didn't map, etc.
fc.typhimurium$stat
# see the top part of the annotation file just because
head(fc.typhimurium$annotation)
write.table(x = data.frame(fc.typhimurium$annotation[,c("GeneID", "Length")], fc.typhimurium$counts, stringsAsFactors = FALSE), file = "typhimurium_featurecounts.txt", append = F, quote = F, sep = "\t", row.names = F)

# merge Cerro and Javiana counts, including pan genome
# genes present in one reference strain but not the other will be denoted by NA
cj <- merge(fc.cerro$counts, fc.javiana$counts, by="row.names", all.x=T, all.y=T)
rownames(cj) <- cj$Row.names
cj <- cj[,2:ncol(cj)]
cj[is.na(cj)] <- 0
table(is.na(as.vector(cj)))
#merge Cerro and Typhimurium counts, including pan genome
ct <- merge(fc.cerro$counts, fc.typhimurium$counts, by="row.names", all.x=T, all.y=T)
rownames(ct) <- cj$Row.names
ct <- ct[,2:ncol(ct)]
ct[is.na(ct)] <- 0
table(is.na(as.vector(ct)))

# merge Cerro and Javiana annotation, including pan genome
# genes present in one reference strain but not the other will be denoted by NA
cj <- merge(fc.cerro$annotation, fc.javiana$annotation, by="row.names", all.x=T, all.y=T)
rownames(cj) <- cj$Row.names
cj <- cj[,2:ncol(cj)]
cj[is.na(cj)] <- 0
table(is.na(as.vector(cj)))

#merge Cerro and Typhimurium annotation, including pan genome
ct <- merge(fc.cerro$annotation, fc.typhimurium$annotation, by="row.names", all.x=T, all.y=T)
rownames(ct) <- cj$Row.names
ct <- ct[,2:ncol(ct)]
ct[is.na(ct)] <- 0
table(is.na(as.vector(ct)))

# merge CerroJaviana and Typhimurium counts, including pan genome
# genes absent from a reference strain will be denoted by NA
cjt <- merge(cj, fc.typhimurium$counts, by="row.names", all.x=T, all.y=T)
# sanity check; how many NAs do we have?
table(is.na(as.vector(cjt)))
# replace NA with 0 counts
cjt[is.na(cjt)] <- 0
# sanity check; are our NAs gone?
table(is.na(as.vector(cjt)))
# change the column names of the count matrix to make it more readable
colnames(cjt) <- c("Feature", "Cerro-rep1", "Cerro-rep2", "Cerro-rep3", 
                   "Javiana-rep1", "Javiana-rep2", "Javiana-rep3", 
                   "Typhimurium-rep1", "Typhimurium-rep2", "Typhimurium-rep3")
# write this zero-inflated pan-genome count matrix to a tab-delimited file
# file is named CerroJavianaTyphimurium_featureCounts_panGenome_as_zeros.txt
write.table(x = cjt, file = "CerroJavianaTyphimurium_featureCounts_panGenome_as_zeros.txt", 
            append = F, quote = F, sep = "\t", row.names = F, col.names = T)

# merge CerroJaviana and Typhimurium annotations, including pan genome
# genes absent from a reference strain will be denoted by NA
cjt <- merge(cj, fc.typhimurium$annotation, by="row.names", all.x=T, all.y=T)
# sanity check; how many NAs do we have?
table(is.na(as.vector(cjt)))
# replace NA with 0 counts
cjt[is.na(cjt)] <- 0
# sanity check; are our NAs gone?
table(is.na(as.vector(cjt)))
# change the column names of the count matrix to make it more readable
colnames(cjt) <- c("Feature", "Cerro-rep1", "Cerro-rep2", "Cerro-rep3", 
                   "Javiana-rep1", "Javiana-rep2", "Javiana-rep3", 
                   "Typhimurium-rep1", "Typhimurium-rep2", "Typhimurium-rep3")
# write this zero-inflated pan-genome count matrix to a tab-delimited file
# file is named CerroJavianaTyphimurium_featureCounts_panGenome_as_zeros.txt
write.table(x = cjt, file = "CerroJavianaTyphimurium_featureCounts_panGenome_as_zeros.txt", 
            append = F, quote = F, sep = "\t", row.names = F, col.names = T)

# merge Cerro and Javiana counts, excluding pan genome
# only genes present in both reference strains will be maintained
cj.core <- merge(fc.cerro$counts, fc.javiana$counts, by="row.names", all.x=F, all.y=F)
rownames(cj.core) <- cj.core$Row.names
table(is.na(as.vector(cj.core)))
colnames(cj.core) <- c("Feature", "Cerro-rep1", "Cerro-rep2", "Cerro-rep3", "Javiana-rep1", "Javiana-rep2", "Javiana-rep3")
write.table(x = cj.core, file="cerrojaviana_featurecounts_core.txt", append = F, quote = F, sep = "\t", row.names = F, col.names = T)

# merge Cerro and Javiana annotation, excluding pan genome
# only genes present in both reference strains will be maintained
cj.core <- merge(fc.cerro$annotation, fc.javiana$annotation, by="row.names", all.x=F, all.y=F)
rownames(cj.core) <- cj.core$Row.names
table(is.na(as.vector(cj.core)))
colnames(cj.core) <- c("Feature", "Cerro-rep1", "Cerro-rep2", "Cerro-rep3", "Javiana-rep1", "Javiana-rep2", "Javiana-rep3")
write.table(x = cj.core, file="cerrojaviana_annotation_core.txt", append = F, quote = F, sep = "\t", row.names = F, col.names = T)

# merge Cerro and Typhimurium counts, excluding pan genome
# only genes present in both reference strains will be maintained
ct.core <- merge(fc.cerro$counts, fc.typhimurium$counts, by="row.names", all.x=F, all.y=F)
rownames(ct.core) <- ct.core$Row.names
table(is.na(as.vector(ct.core)))
colnames(ct.core) <- c("Feature", "Cerro-rep1", "Cerro-rep2", "Cerro-rep3", "Typhimurium-rep1", "Typhimurium-rep2", "Typhimurium-rep3")
write.table(x = ct.core, file="cerrotyphimurium_featurecounts_core.txt", append = F, quote = F, sep = "\t", row.names = F, col.names = T)

# merge Javiana and Typhimurium counts, excluding pan genome
# only genes present in both reference strains will be maintained
jt.core <- merge(fc.javiana$counts, fc.typhimurium$counts, by="row.names", all.x=F, all.y=F)
rownames(jt.core) <- jt.core$Row.names
table(is.na(as.vector(jt.core)))
colnames(jt.core) <- c("Feature", "Javiana-rep1", "Javiana-rep2", "Javiana-rep3", "Typhimurium-rep1", "Typhimurium-rep2", "Typhimurium-rep3")
write.table(x = jt.core, file="javianatyphimurium_featurecounts_core.txt", append = F, quote = F, sep = "\t", row.names = F, col.names = T)

# merge CerroJaviana and Typhimurium counts, excluding pan genome
# only genes present in all three reference strains will be maintained
cjt.core <- merge(cj.core, fc.typhimurium$counts, by="row.names", all.x=F, all.y=F)
# sanity check; there should be no NAs (this is just the core genome of the 3 strains)
table(is.na(as.vector(cjt.core)))
# change the column names of the count matrix to make it more readable
colnames(cjt.core) <- c("Feature", "Cerro-rep1", "Cerro-rep2", "Cerro-rep3", 
                   "Javiana-rep1", "Javiana-rep2", "Javiana-rep3", 
                   "Typhimurium-rep1", "Typhimurium-rep2", "Typhimurium-rep3")
# write this core-genome count matrix to a tab-delimited file
# file is named CerroJavianaTyphimurium_featureCounts_drop_panGenome.txt
write.table(x = cjt.core, file = "CerroJavianaTyphimurium_featureCounts_drop_panGenome.txt", 
            append = F, quote = F, sep = "\t", row.names = F, col.names = T)

# merge CerroJaviana and Typhimurium annotation, excluding pan genome
# only genes present in all three reference strains will be maintained
cjt.core <- merge(cj.core, fc.typhimurium$annotation, by="row.names", all.x=F, all.y=F)
# sanity check; there should be no NAs (this is just the core genome of the 3 strains)
table(is.na(as.vector(cjt.core)))
# change the column names of the count matrix to make it more readable
colnames(cjt.core) <- c("Feature", "Cerro-rep1", "Cerro-rep2", "Cerro-rep3", 
                        "Javiana-rep1", "Javiana-rep2", "Javiana-rep3", 
                        "Typhimurium-rep1", "Typhimurium-rep2", "Typhimurium-rep3")
# write this core-genome count matrix to a tab-delimited file
# file is named CerroJavianaTyphimurium_featureCounts_drop_panGenome.txt
write.table(x = cjt.core, file = "cjt_annotation.txt", 
            append = F, quote = F, sep = "\t", row.names = F, col.names = T)
