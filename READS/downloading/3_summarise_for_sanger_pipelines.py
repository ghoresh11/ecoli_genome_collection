import os
import gzip
import sys



loc =  os.path.abspath(sys.argv[1]) ## location of the trimmed reads to import to farm
orig_loc = os.path.abspath(sys.argv[2]) ## locations of the original reads (before trimming)
batch = sys.argv[3]



out = open("HoreshGal_ecoli_mobilome_" + batch + ".csv", "w")
out.write("Supplier Name,NCBI\nSupplier Organisation,NCBI\nSanger Contact Name,Gal Horesh\n\
	Sequencing Technology,Illumina\n")
out.write("Study Name, ecoli_mobilome_" + str(batch) + "\n")
out.write("Study Accession number,\nTotal size of files in GBytes,\nData to be kept until,31.12.21\n")
out.write("Path: " + loc + "\n")
out.write("Filename,Mate File,Sample Name,Sample Accession number,Taxon ID,Library Name,Fragment Size,Read Count,Base Count,Comments\n")


fastq_files = os.listdir(loc)


failed = [] ## add here numbers of download files that failed

## move and delete files that failed (you can save files that failed trimming if you want to work on them later)
failed_loc = "failed/"
for name in failed:
	pair1_trimmed =  os.path.join(loc,name + "_trimmed_1.fastq.gz")
	pair2_trimmed =  os.path.join(loc, name + "_trimmed_2.fastq.gz")
	pair1 = os.path.join(orig_loc, name + "_1.fastq.gz")
	pair2 = os.path.join(orig_loc, name + "_2.fastq.gz")
	
	## delete the trimmed files 
	if os.path.isfile(pair1_trimmed):
		os.remove(pair1_trimmed)
	if os.path.isfile(pair2_trimmed):
		os.remove(pair2_trimmed)

	## move original files to other directory	os.rename("../" + pair1, long_loc + pair1)
	if os.path.isfile(pair1):
		os.rename(pair1, os.path.join(failed_loc, name + "_1.fastq.gz"))
	if os.path.isfile(pair2):
		os.rename(pair2, os.path.join(failed_loc, name + "_2.fastq.gz"))


fastq_files = os.listdir(loc) ## relook up all the files (without the failed files)

for f in fastq_files:
	if not f.endswith(".fastq.gz"):
		continue
	toks = f.split("_")
	if toks[2].startswith("2"):
		continue ## the other pair
	name = toks[0]

	pair1 =  name + "_trimmed_1.fastq.gz"
	pair2 =  name + "_trimmed_2.fastq.gz"

	if not pair2 in fastq_files:
		pair2 = "" 


	out.write(",".join(map(str, 
		[pair1, pair2, name, name, 562, "", "-", "-", 
		"-", ""])) + "\n")

out.close()