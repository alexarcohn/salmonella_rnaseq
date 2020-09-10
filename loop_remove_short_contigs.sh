#!/bin/bash
#jk2739
#Loop through the contig files and remove contigs shorter than 200 bp
#Usage: Run the script from the directory with the loop_remove_short_contigs.sh, remove_short_contigs.py, contigs fasta files
#Usage: sh loop_remove_short_contigs.sh .


for f in *_contigs.fasta
do
python remove_short_contigs.py $f ${f%.fasta}_long.fasta
done
