#!/usr/bin/env python2

from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import translate
from Bio.Seq import reverse_complement
import os
import string
import random
import sys


''' given all the gff files, summarise them to create a big CSV file
with all the details of these genomes
the details:
1. Number of CDSs
2. Number of pseudogene
3. Number of other elements in GFF files
4. Genome size
After summary -> add ST, pathotype and metadata and see how stratifying the data changes anything
and will give a general description of the pan-genome'''


def read_gff(gff_file):
	categories = ["hypothetical protein", "transposase",
	    "pseudogene", "conjuga", "phage", "fimbrial", "plasmid",
	     "crispr", "resistance", "virulence", "secretion system"]
	counts = {}
	
	for cat in categories:
		counts[cat] = 0
	try:
		f = open(gff_file)
	except IOError:
		print("Could not read file:", gff_file)
		return counts

	for line in f:
		line = line.lower()
		if line.startswith("##fasta"):
			break
		if line.startswith("#"):
			continue
		toks = line.strip().split()
		product = toks[2]
		if product not in counts:
			counts[product] = 0
		counts[product] += 1
		for cat in categories:
			if cat in line:
				counts[cat] += 1
	
	f.close()

	return counts

header = ["cds", "trna", "hypothetical protein", "transposase",
    "pseudogene", "conjuga", "phage", "fimbrial", "plasmid",
     "crispr", "resistance", "virulence", "secretion system"]

out = open("gff_summary.csv", "w")
out.write("ID, file_name," + ",".join(header) + "\n")

cnt = 0
with open(sys.argv[1]) as f:
	for line in f:
		toks = line.strip().split("\t")
		if line.startswith("ID"):
			annot_loc_index = toks.index("Annotation_Location")
			continue
		ID = toks[0]
		files = toks[annot_loc_index].split(",")
		for f1 in files:
			print(f1)
			counts = read_gff(f1)
			out.write(ID + "," + f1)
			for cat in header:
				out.write("," + str(counts[cat]))
			out.write("\n")
		cnt += 1


out.close()