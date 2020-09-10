#!/bin/bash
# average_coverage.sh <path to directory with contigs and reads>
# written November 2, 2015 by LC Carroll
# Shout out to Matt S. for giving me the bam_coverage.sh script

cd $1
# BBMap to determine coverage
for f in *_contigs_long.fasta
do
echo "Indexing $f with BBMap..."
/programs/bbmap-37.50/bbmap.sh ref=$f
echo "Mapping reads to $f with BBMap..." #the suffix if not dealing with trimmedP files
/programs/bbmap-37.50/bbmap.sh in=${f%_contigs_long.fasta}_1.trimmedP.fastq.gz in2=${f%_contigs_long.fasta}_2.trimmedP.fastq.gz out=${f%_contigs_long.fasta}.sam
echo "SAM file created.  BBMap finished."
mv ref/ ${f%_contigs_long.fasta}_ref/

# Now let's use samtools to covert, sort, and index
echo "Converting SAM to BAM with samtools..."
samtools view -Sb ${f%_contigs_long.fasta}.sam > ${f%_contigs_long.fasta}.bam
echo "BAM file created."
echo "Removing sam file..."
rm -r *.sam
echo "Sorting BAM file with samtools..."
samtools sort ${f%_contigs_long.fasta}.bam -o ${f%_contigs_long.fasta}_sorted.bam
echo "Finished sorting."
echo "Indexing sorted BAM file..."
samtools  index ${f%_contigs_long.fasta}_sorted.bam
echo "Index complete."
echo "Using samtools depth to obtain average genome coverage..."
X=$(samtools depth ${f%_contigs_long.fasta}_sorted.bam | awk '{sum+=$3} END { print sum/NR}');
echo "${f%_contigs_long.fasta}_sorted.bam";
echo "$X";
echo "${f%_contigs_long.fasta}_sorted.bam $X">> average_coverage.txt;
done

