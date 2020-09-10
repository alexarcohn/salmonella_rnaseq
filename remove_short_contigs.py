#!/bin/python
#Remove contigs shorter than 200 bp
#jk2739
#092716
#Usage: run from a directory with contigs file and a script
#Usage: python remove_short_contigs.py <infile.fasta> <outfile.fasta>

import sys
from Bio import SeqIO

infile = sys.argv[1]
parsed_infile = SeqIO.parse(open(infile,"rU"), "fasta")
remove_short = (contig for contig in parsed_infile if len(contig.seq) > 500)

outfile= sys.argv[2]
output = open(outfile, "w")
SeqIO.write(remove_short, output, "fasta")
output.close()
