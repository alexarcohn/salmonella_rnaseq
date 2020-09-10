#!/bin/bash
#Matt Stasiewicz 7-1-14
#Use Trimmomatic to trim the raw reads
#sh trimmomatic2.sh <inpath>

#$1=/media/drive2/NYSDOH_ENV_SRA1/MJS

#O: *.trimmed[S/P].fastq.gz files in $gp
#deletes .fastq.gz

#loops through the output from fastq, the _[1/2].fastq.gz and does read trimming
#Raw read trimming with Trimmomatic, ref below.  MJS comments each step of the loop with documentation from:
#http://www.usadellab.org/cms/?page=trimmomatic
#It appears Henk used the default parameters settings for all steps, adjusting file paths and names as appropriate

cd $1
echo | pwd
for f in *1.fastq.gz
	do 
		if [ -f "${f%_1.fastq.gz}_1.trimmedP.fastq.gz" ]
		then
		echo 'skip'${f}
		continue
		fi
	echo 'trim' ${f}
	java -jar /programs/trimmomatic/trimmomatic-0.36.jar SE -threads 30 -phred33 -trimlog log SRR653931_1.fastq.gz SRR653931_1.trimmedP.fastq.gz ILLUMINACLIP:/programs/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:26;
done;

#rm *_1.fastq.gz
#rm *_2.fastq.gz

