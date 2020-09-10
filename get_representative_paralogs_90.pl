#!usr/bin/perl
 

print ("Cluster\tNumber_of_Genes\tNumber_of_Taxa\tRepresentative\n");                 
while(defined($line=<>)){

if ($line=~m/\S+,\d+,\d+,+\w+_\d+,+/){

$line=~ m/^(\S+),(\d+),(\d+),+(\w+_\d+),+/;
$cluster= $1;
$number_genes= $2;
$number_taxa= $3;
$representative= $4;

print ("" . $cluster . "\t" . $number_genes . "\t" . $number_taxa . "\t" . $representative . "\n");
}
}


