from Bio.SeqIO.FastaIO import SimpleFastaParser
import os
import re
## change the names of files downloaded from NCBI so they match the accessions I have

dr = "ExPEC/"

fasta_files = os.listdir(dr + "others0000/")
gff_files = os.listdir(dr + "others0000/")


# # first remove the older version if there is more than one
# versions = {}
# for f in fasta_files:  ## switch between GFF and FASTA
# 	if not f.endswith(".fna"):
# 		continue
# 	toks = f.split("v")
# 	name = toks[0].split(".1")[0]
# 	name = name.split(".2")[0]
# 	name = name.split(".3")[0]
# 	print(toks)
# 	v = toks[1][0]
# 	if name not in versions:
# 		versions[name] = []
# 	versions[name].append(f)


# for v in versions:
# 	if len(versions[v]) > 1:
# 		print(versions[v])
# 		# for f in versions[v]:
# 		# 	if "v1" in f:
# 		# 		os.remove(dr + "fastas-2018-07-27/" + f)
# 		# 	elif "v2" in f and len(versions[v]) > 2:
# 		# 		os.remove(dr + "fastas-2018-07-27/" + f)
# #print(versions)
# quit()

## list files again after duplocates have been removed
# fasta_files = os.listdir(dr + "fastas-2018-07-27/")
# gff_files = os.listdir(dr + "gffs-2018-07-27/")

def change_identifier(identifier):
	print("** " + identifier)
	if identifier.startswith("NZ"):
		identifier = identifier.split("_")[-1]
	identifier = identifier.replace(".1","")
	identifier = identifier.replace(".2","")
	identifier = identifier.replace(".3","")
	if identifier.startswith("A"):
		identifier = re.sub('\d', '0', identifier)
	return identifier


for f in gff_files:
	if not f.endswith(".gff"):
		continue
	with open(dr + "others0000/" + f) as handle:
		for line in handle:
			if line.startswith("#"):
				continue
			identifier = line.strip().split()[0]
			identifier = change_identifier(identifier)
			break
		print(identifier)
		os.rename(dr + "others0000/" + f, dr + "others0000/" + identifier + ".gff")

for f in fasta_files:
	if not f.endswith(".fna"):
		continue
	with open(dr + "others0000/" + f) as handle:
		for values in SimpleFastaParser(handle):
			identifier = values[0].split()[0]
			identifier = change_identifier(identifier)
			break
	print(identifier)
	os.rename(dr + "others0000/" + f, dr + "others0000/" + identifier + ".fa")
