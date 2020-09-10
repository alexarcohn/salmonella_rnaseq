#!/usr/bin/python
# Laura M. Carroll
# April 22, 2019
# using roary's clustered_proteins assignments and gff files, create a saf file which can be used with featureCounts R package
# note that this only uses genes in the 3 reference strains (only protein IDs associated with the 3 reference genomes are queried)
# move the script to a directory with clustered_proteins file, as well as all gff files
# py gff2saf.py

# import required packages
import sys, os, glob

# dictionary relating each protein ID to its orthologous cluter (gene) name assigned in Roary
d = {}
# dictionary relating each roary gene name to the serotypes in which it is found (only queries the three reference strains, i.e., the genomes with gff files)
serotype = {}
# open clustered_proteins
infile = open("clustered_proteins", "r")
# loop through each roary orthologous cluster
for line in infile:
	# get the name of the gene
	gene = line.split(":")[0].strip()
	# get the protein IDs for that gene
	ids = line.split(":")[1].strip()
	add each protein ID to the d dictionary, with the roary gene treated as a value
	for i in ids.split():
		d[i.strip()] = gene
	# add the roary gene to the serotype dictionary
	serotype[gene] = []
	# if the Typhimurium reference gff protein ID prefix is listed in the protein IDs, include Typhimurium in the dictionary
	if "EEOABMNO" in ids:
		serotype[gene].append("Typhimurium")
	# if the Cerro reference gff protein ID prefix is listed in the protein IDs, include Cerro in the dictionary
	if "AJGAFMDA" in ids:
		serotype[gene].append("Cerro")
	# if the Javiana reference gff protein ID prefix is listed in the protein IDs, include Javiana in the dictionary
	if "IBJCNMFN" in ids:
		serotype[gene].append("Javiana")
# close clustered_proteins
infile.close()


# loop through each reference gff file
for f in glob.glob("*.gff"):
	# make a saf file for each gff file, adding column headers
	outfile = open(f.replace(".gff", ".saf"), "a")
	print >> outfile, "\t".join(["GeneID", "Chr", "Start", "End", "Strand", "Annotation", "Phenotype"])
	outfile.close()
	# open the gff file
	infile = open(f, "r")
	# loop through the lines of the gff file
	for line in infile:
		# if the line is not commented out
		if "##" not in line:
			# get the chromosome (or contig ID)
			Chr = line.split("\t")[0].strip()
			# get the start of the gene
			Start = line.split("\t")[3].strip()
			# get the end of the gene
			End = line.split("\t")[4].strip()
			# get the strand of the gene
			Strand = line.split("\t")[6].strip()
			# get the gene annotation, just because
			Annotation = line.split("\t")[-1].strip()
			# get the protein ID of the gene for the reference strain
			testGeneID = Annotation.split(";")[0].strip()
			testGeneID = testGeneID.replace("ID=", "")
			# check if the protein ID has a roary gene name associated with it
			try:
				GeneID = d[testGeneID]
				# assign a "phenotype" to the gene
				# if the gene is present in all 3 reference strains, the phenotype is "CerroJavianaTyphimurium"
				# if the gene is present in 2 reference strains, the phenotype is "CerroJaviana", "CerroTyphimurium", or "JavianaTyphimurium"
				# if the gene is present in 1 reference strain, the phenotype is "Cerro", "Javiana", or "Typhimurium"
				phenotype = sorted(serotype[GeneID])
				phenotype = "".join(phenotype)
				# open the saf file for the reference genome and print the roary gene ID, chromosome ID, start, end, strand, annotation, and phenotype of each gene
				outfile = open(f.replace(".gff", ".saf"), "a")
				print >> outfile, "\t".join([GeneID, Chr, Start, End, Strand, Annotation, phenotype])
				outfile.close()
			# if there is no roary gene name associated with a protein ID (e.g., the ID is associated with a tRNA or other RNA), write to a file ending in .errors
			except KeyError:
				outfile = open(f.replace(".gff", ".errors"), "a")
				print >> outfile, line.strip()
				outfile.close()
		# ignore the fasta part of the gff file	
		else:
			if "##FASTA" in line:
				break
	infile.close()	
