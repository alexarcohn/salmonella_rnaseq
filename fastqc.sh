#!/bin/bash
#Matt Stasiewicz 3-7-14
#Use fastq to perform intial qc and save results
#nohup sh <inpath> <outpath> <suff>

#$1=/media/drive2/NYSDOH_ENV_SRA1/MJS
#$2=/media/drive2/NYSDOH_ENV_SRA1/MJS/fastqc_res
#$3=[1-2].fastq.gz

#O: fastqc html in $gp/fastqc_res

cd $1
for file in *$3
    do 
	cd $2	
	if [ -d "${file%.fastq.gz}_fastqc" ]
	then
	echo 'skip'${file}
	continue
	fi
	cd $1
    echo 'fastqc 1'${file}
    fastqc $file -o $2;
done
