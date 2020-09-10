#!/bin/bash
#QUAST - assembly quality control
#jk2739
#112415

mkdir quast_results

for f in *.fasta
do
python /programs/quast-4.0/quast.py -o ./quast_results/quast_${f%_contigs.fasta} --min-contig 1 $f
done

#collect report txt files

mkdir quast_reports
for f in *.fasta
do cd quast_results/quast_${f%_contigs.fasta}
cat report.txt > ${f%_contigs.fasta}_report.txt
cp ${f%_contigs.fasta}_report.txt ../../quast_reports
cd ../..
done


