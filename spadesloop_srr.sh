#!/bin/bash

#Matt Stasiewicz 7-1-14 modified by LC Carroll 12-18-14 by S Harrand 10-12-16 by jk2739 10-26-16
#Use SPades to assemble the genome
#nohup sh spades.sh <inpath>

#run spades
python /programs/spades/bin/spades.py -s SRR653931_1.trimmedP.fastq.gz -o SRR653931_1.trimmedP.fastq.gz;
done
#check the created log file for any issues

#collect contigs files and rename them
mkdir contigs  
cd SRR653931_1.trimmedP.fastq.gz}
cat contigs.fasta > SRR653931_1.trimmedP.fastq.gz_contigs.fasta
cp SRR653931_1.trimmedP.fastq.gz_contigs.fasta ../contigs
cd ..

mkdir scaffolds  
cd SRR653931_1.trimmedP.fastq.gz
cat scaffolds.fasta > SRR653931_1.trimmedP.fastq.gz_scaffolds.fasta
cp SRR653931_1.trimmedP.fastq.gz_scaffolds.fasta ../scaffolds
cd ..