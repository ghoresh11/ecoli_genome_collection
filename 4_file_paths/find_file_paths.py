from dirs import *
import os.path


## on the way check that all the files that have a "NO" for artificial reads actually exist on the farm
## if I can't find them, then understand why I can't (maybe deleted them in the past because of not being similar enough)
## then, all the ones that I do download, should be downloaded using fastq-dump, some should point back to accessions that
## I already have their assemblies, others need to be assembled and annotated!

# out = open("require_artificial.csv", "w")

out2 = open("metadata_fixed_with_loc.csv", "w")

fasta_suffixes =  [".fa", ".fasta", "_trimmed.contigs_velvet.fa", "_contigs_pacbio.fa",  ".contigs_velvet.fa","_contigs_canu_1_6.fa", "_contigs_hgap_4_0.fa"]
gff_suffixes = [".gff", ".fasta.gff", ".contigs_velvet.fa.gff", "_trimmed.contigs_velvet.fa.gff", "_trimmed.gff", "_trimmed.velvet.gff",".fa.gff", "_contigs_pacbio.fa.gff", "_contigs_hgap_4_0.fa.gff", "_contigs_canu_1_6.fa.gff"]
reads_suffixes = [".fastq.gz", "_1.fastq", "_1.fastq.gz","_trimmed_1.fastq.gz","_2.fastq", "_2.fastq.gz","_trimmed_2.fastq.gz"]

### 

def find_path(name, suffixes, dirs):
	all_files = set()
	for d in dirs:
		check =  os.path.join(d,name)
		for s in suffixes:
			if os.path.isfile(check + s):
				all_files.add(check + s) ## keep track of all files meeting the requirements

	all_files = list(all_files)
	if len(all_files) == 0:
		return None

	return all_files

with open("metadata_fixed_with_reads.csv") as f:
	for line in f:
		if line.startswith("ID"):
			out2.write(line.strip() + "\tAssembly_Location\tAnnotation_Location\tReads_Location\n")
			continue

		flag = False
		toks = line.strip().split("\t")
		
		name1 = toks[0]  ## file could have either identifier (confusing...)
		name2 = toks[1]
		name3 = toks[-2].split(",")
		


		all_names = [name1.upper(), name2.upper(), name1, name2]
		for n in name3:
			all_names.append(n.upper())
			all_names.append(n.lower())
		all_names = list(set(all_names))

		# if "7067_5#4" not in all_names:  ## if debugging change here to focus on one example
		# 	continue

		### Search for the assembly file
		loc = ["Not found"]
		for name in all_names:
			curr_loc = find_path(name, fasta_suffixes, all_dirs_fasta)
			if curr_loc != None and loc == ["Not found"]:
				loc = set(curr_loc) ## enable more than a single locations
			elif curr_loc != None:
				for item in curr_loc:
					loc.add(item)

		if toks[-1] == "No" and loc == ["Not found"]: ## we do have the read files for this genome
			loc = ["No assembly"]		 ## but couldn't find an assembly


		## Search for annotation file
		loc_annot = ["No annotation"]
		for name in all_names:
			curr_loc = find_path(name, gff_suffixes, all_dirs_gff)
			if curr_loc != None and loc_annot == ["No annotation"]:
				loc_annot = set(curr_loc)
			elif curr_loc != None:
				for item in curr_loc:
					loc_annot.add(item)

		## Search for the read files
		loc_reads = ["No reads"]

		reads = toks[-2].split(",") ## get all the reads
		reads = reads + [toks[0]] ## sometimes the read files have a different name from the SRR 
		
		names = set()
		for r in reads:
			names.add(r.lower())
			names.add(r.upper())	
		names = list(names)
		for name in names:
			curr_loc = find_path(name, reads_suffixes, all_dirs_reads)
			if curr_loc != None and loc_reads == ["No reads"]:
				loc_reads = set(curr_loc)
			elif curr_loc != None:
				for item in curr_loc:
					loc_reads.add(item)


		## write all the locations to a file
		out2.write(line.strip() + "\t" + ",".join(list(loc)) + "\t" + \
					",".join(list(loc_annot)) + "\t" + ",".join(list(loc_reads)) + "\n")


out2.close()