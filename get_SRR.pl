#!usr/bin/perl
#This is a script to go over several files with SRA IDs and then download SRA files usign fastq-dump

#Get all files with SRA Ids from directory
my @txt = glob "srr_cantaloupe*.txt";

#read each file at a time
foreach $filename (@txt){

#open file using file handler
open INPUT, "<", $filename;

#check each line at a time                         
while(defined($line=<INPUT>)){

#chomp($line);

#run fastq-dump

open FASTQ, "|fastq-dump --gzip --split-files $line";
close FASTQ;
}
close INPUT;
}

