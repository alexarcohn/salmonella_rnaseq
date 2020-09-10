#!/usr/bin/env Rscript
# ABSTRACT: Create R plots
# PODNAME: create_plots.R
# Take the output files from the pan genome pipeline and create nice plots.

library(ggplot2)

newgenes = read.table("~/Documents/RNA_seq_Salmonella/javiana/roary/number_of_new_genes_javiana.Rtab")
boxplot(newgenes, data=newgenes, main="Number of new genes",
        xlab="No. of genomes", ylab="No. of genes",varwidth=TRUE, ylim=c(0,max(newgenes)), outline=FALSE)

conservedgenes = read.table("~/Documents/RNA_seq_Salmonella/javiana/roary/number_of_conserved_genes_javiana.Rtab")
boxplot(conservedgenes, data=conservedgenes, main="Number of conserved genes",
        xlab="No. of genomes", ylab="No. of genes",varwidth=TRUE, ylim=c(0,max(conservedgenes)), outline=FALSE)

pangenome_genes = read.table("~/Documents/RNA_seq_Salmonella/javiana/roary/number_of_genes_in_pan_genome_javiana.Rtab")
boxplot(pangenome_genes, data=pangenome_genes, main="No. of genes in the pan-genome",
        xlab="No. of genomes", ylab="No. of genes",varwidth=TRUE, ylim=c(0,max(pangenome_genes)), outline=FALSE)

uniquegenes = read.table("~/Documents/RNA_seq_Salmonella/javiana/roary/number_of_unique_genes_javiana.Rtab")
boxplot(uniquegenes, data=uniquegenes, main="Number of unique genes",
        xlab="No. of genomes", ylab="No. of genes",varwidth=TRUE, ylim=c(0,max(uniquegenes)), outline=FALSE)

blast = read.table("~/Documents/RNA_seq_Salmonella/javiana/roary/blast_identity_frequency_javiana.Rtab")
plot(blast,main="Number of blastp hits with different percentage identity",  xlab="Blast percentage identity", ylab="No. blast results")


library(ggplot2)
conserved = colMeans(read.table("~/Documents/RNA_seq_Salmonella/javiana/roary/number_of_conserved_genes_javiana.Rtab"))
total = colMeans(read.table("~/Documents/RNA_seq_Salmonella/javiana/roary/number_of_genes_in_pan_genome_javiana.Rtab"))

genes = data.frame( genes_to_genomes = c(conserved,total),
                    genomes = c(c(1:length(conserved)),c(1:length(conserved))),
                    Key = c(rep("Conserved genes",length(conserved)), rep("Total genes",length(total))) )

ggplot(data = genes, aes(x = genomes, y = genes_to_genomes, group = Key, linetype=Key)) +geom_line()+
  theme_classic() +
  ylim(c(3500,max(total)))+
  xlim(c(1,length(total)))+
  xlab("No. of genomes") +
  ylab("No. of genes")+ theme_bw(base_size = 16)
######################

unique_genes = colMeans(read.table("~/Documents/RNA_seq_Salmonella/Cerro/Roary/number_of_unique_genes_cerro.Rtab"))
new_genes = colMeans(read.table("~/Documents/RNA_seq_Salmonella/Cerro/Roary/number_of_new_genes_cerro.Rtab"))

genes = data.frame( genes_to_genomes = c(unique_genes,new_genes),
                    genomes = c(c(1:length(unique_genes)),c(1:length(unique_genes))),
                    Key = c(rep("Unique genes",length(unique_genes)), rep("New genes",length(new_genes))) )

ggplot(data = genes, aes(x = genomes, y = genes_to_genomes, group = Key, linetype=Key)) +geom_line()+
  theme_classic() +
  ylim(c(1,max(unique_genes)))+
  xlim(c(1,length(unique_genes)))+
  xlab("No. of genomes") +
  ylab("No. of genes")+ theme_bw(base_size = 16)

conserved = colMeans(read.table("~/number_of_conserved_genes_cerro.Rtab"))
total = colMeans(read.table("~/number_of_genes_in_pan_genome_cerro.Rtab"))

genes = data.frame( genes_to_genomes = c(conserved,total),
                    genomes = c(c(1:length(conserved)),c(1:length(conserved))),
                    Key = c(rep("Conserved genes",length(conserved)), rep("Total genes",length(total))) )

ggplot(data = genes, aes(x = genomes, y = genes_to_genomes, group = Key, linetype=Key)) +geom_line()+
  theme_classic() +
  ylim(c(3500,max(total)))+
  xlim(c(1,16))+
  xlab("No. of genomes") +
  ylab("No. of genes")+ theme_bw(base_size = 16)
#### javiana ####

unique_genes_javiana = colMeans(read.table("~/Documents/RNA_seq_Salmonella/javiana/number_of_unique_genes.Rtab"))
new_genes_javiana = colMeans(read.table("~/Documents/RNA_seq_Salmonella/javiana/number_of_new_genes.Rtab"))

genes = data.frame( genes_to_genomes = c(unique_genes_javiana,new_genes_javiana),
                    genomes = c(c(1:length(unique_genes_javiana)),c(1:length(unique_genes_javiana))),
                    Key = c(rep("Unique genes",length(unique_genes_javiana)), rep("New genes",length(new_genes_javiana))) )

library(ggplot2)
ggplot(data = genes, aes(x = genomes, y = genes_to_genomes, group = Key, linetype=Key)) +geom_line()+
  theme_classic() +
  ylim(c(1,max(unique_genes_javiana)))+
  xlim(c(1,length(unique_genes_javiana)))+
  xlab("No. of genomes") +
  ylab("No. of genes")+ theme_bw(base_size = 16) +  theme(legend.justification=c(1,1),legend.position=c(1,1))

conserved_javiana = colMeans(read.table("~/Documents/RNA_seq_Salmonella/javiana/number_of_conserved_genes.Rtab"))
total_javiana = colMeans(read.table("~/Documents/RNA_seq_Salmonella/javiana/number_of_genes_in_pan_genome.Rtab"))

genes_javiana = data.frame( genes_to_genomes = c(conserved_javiana,total_javiana),
                    genomes = c(c(1:length(conserved_javiana)),c(1:length(conserved_javiana))),
                    Key = c(rep("Conserved genes",length(conserved_javiana)), rep("Total genes",length(total_javiana))) )

ggplot(data = genes_javiana, aes(x = genomes, y = genes_to_genomes, group = Key, linetype=Key)) +geom_line()+
  theme_classic() +
  ylim(c(3500,max(total_javiana)))+
  xlim(c(1,length(total_javiana)))+
  xlab("No. of genomes") +
  ylab("No. of genes")+ theme_bw(base_size = 16) +  theme(legend.justification=c(0,1),legend.position=c(0,1))

###typhimurium###

conserved_typhimurium = colMeans(read.table("~/Documents/RNA_seq_Salmonella/Typhimurium/Roary/number_of_conserved_genes_typhimurium.Rtab"))
total_typhimurium = colMeans(read.table("~/Documents/RNA_seq_Salmonella/Typhimurium/Roary/number_of_genes_in_pan_genome_typhimurium.Rtab"))

genes_typhimurium = data.frame( genes_to_genomes = c(conserved_typhimurium,total_typhimurium),
                            genomes = c(c(1:length(conserved_typhimurium)),c(1:length(conserved_typhimurium))),
                            Key = c(rep("Conserved genes",length(conserved_typhimurium)), rep("Total genes",length(total_typhimurium))) )

ggplot(data = genes_typhimurium, aes(x = genomes, y = genes_to_genomes, group = Key, linetype=Key)) +geom_line()+
  theme_classic() +
  ylim(c(3500,max(total_typhimurium)))+
  xlim(c(1,length(total_typhimurium)))+
  xlab("No. of genomes") +
  ylab("No. of genes")+ theme_bw(base_size = 16)
########## all ##########

conserved_all = colMeans(read.table("~/Documents/RNA_seq_Salmonella/Roary/number_of_conserved_genes.Rtab"))
total_all = colMeans(read.table("~/Documents/RNA_seq_Salmonella/Roary/number_of_genes_in_pan_genome.Rtab"))

genes_all = data.frame( genes_to_genomes = c(conserved_all,total_all),
                                genomes = c(c(1:length(conserved_all)),c(1:length(conserved_all))),
                                Key = c(rep("Conserved genes",length(conserved_all)), rep("Total genes",length(total_all))) )

ggplot(data = genes_all, aes(x = genomes, y = genes_to_genomes, group = Key, linetype=Key)) +geom_line()+
  theme_classic() +
  ylim(c(3500,max(total_all)))+
  xlim(c(1,length(total_all)))+
  xlab("No. of genomes") +
  ylab("No. of genes")+ theme_bw(base_size = 16)
