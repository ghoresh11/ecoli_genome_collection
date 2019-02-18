
out = open("run_shredder_missing2.sh", "w")

with open("/Users/gh11/e_colis/genomes/4_file_paths/missing/artificial_reads.txt") as f:
	for line in f:
		toks = line.strip().split("/")
		filename = line.strip()
		outfile = toks[-1].split(".")[0]
		out.write("\t".join(["/nfs/pathogen/sh16_scripts/fasta2fastq_shredder.py",
			filename, outfile, "100", "3" , "l", "350"]) + "\n")
out.close()