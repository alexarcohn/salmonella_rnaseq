#!/bin/bash
# sh sistrloop.sh 
# run from your directory with fasta files; deposits results in sistr_out directory
# Laura Carroll lmc297@cornell.edu
# Jan 25, 2019

export PYTHONPATH=/programs/sistr_cmd/lib/python2.7/site-packages/

export PATH=/programs/mafft/bin:$PATH

mkdir sistr_out

for f in *.fasta
do
/programs/sistr_cmd/bin/sistr -i ./$f $f -o sistr_out/${f%.fasta} -f csv --qc -t 4 --verbose 
done
