#!/bin/bash
# Lory henderson
# 01_03_17
# RNA seq analysis
#Edited by Renato Orsi to comply with BioHPC
#sh rna_seq.sh <directory with .fastq.gz files> <reference genome>

cd $1

#Remote access to Zeus: ssh loh9@zeus.foodscience.cornell.edu
#Run bwa to create index for reference genome
#which bwa gives you path (/usr/bin/bwa) to bwa program

bwa index $2 

#Map/align reads to reference genome using BWA MEM
#Add extension: creates a new file type (.sam)

for f in *.fastq.gz
do
bwa mem $2 $f > ${f%.fastq.gz}.sam
 

#Use Samtools to sort files by left most coordinate or by read name (maybe it aligns identical reads?)
#Samtools requires path to the program


samtools view -b -S -u -q 1 ${f%.fastq.gz}.sam | samtools sort -o ${f%.fastq.gz}.sorted.bam 


#Create a file with reads that map to the anti-sense strand from reference genome; b = bam file
#Why 16? 16 stands for reverse strand. -f stands for match this and -F stands for do not match this. -f 16 means match reverse strand
#How does it distinguish sense and antisense

samtools view -f 16 ${f%.fastq.gz}.sorted.bam -b -o ${f%.fastq.gz}_as.sorted.bam 

samtools view -F 16 ${f%.fastq.gz}.sorted.bam -b -o ${f%.fastq.gz}_s.sorted.bam

#Index the bam file:
#Creates index for sorted sense and antisense strand BAM file
#You have to do everything twice, for sense and antisense file

samtools index ${f%.fastq.gz}_as.sorted.bam

samtools index ${f%.fastq.gz}_s.sorted.bam

#Create coverage files using BEDTOOLS
#creates a 3-column coverage file; first column is a genome, second is a nucleotide position and third is number of reads matching; 
#This is a coverage file: -d gets coverage for each nucleotide position; -ibam it is reading from bam file -g is the reference genome file
#Can be written as full path to bedtools genomecov or genomeCoverageBed 
#Can be found by googling genomecoveragebed

bedtools genomecov -d -ibam ${f%.fastq.gz}_as.sorted.bam -g $2 > ${f%.fastq.gz}_antisense_cov.txt

bedtools genomecov -d -ibam ${f%.fastq.gz}_s.sorted.bam -g $2 > ${f%.fastq.gz}_sense_cov.txt

#UNIX command (because Artemis cannot open that file)
#Creates a one-column coverage file we will open in Artemis
#Need one column as an input for Bayseq

awk '{ print $3 }' ${f%.fastq.gz}_antisense_cov.txt > ${f%.fastq.gz}_antisense.cov 

awk '{ print $3 }' ${f%.fastq.gz}_sense_cov.txt > ${f%.fastq.gz}_sense.cov 

done

#Notes: Things to check if something goes wrong
#1. Check for correct path lengths
#2.Check for correct reference file name
#3.Check for updated program version, which would lead to a different program name
#4. Check where in the process the script failed




