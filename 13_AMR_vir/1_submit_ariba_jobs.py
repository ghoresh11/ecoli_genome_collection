import subprocess
import os
import shutil

command_prefix = ["/nfs/users/nfs_g/gh11/.local/bin/ariba", "run", "ariba_db"] ## + two files plus output file
mem = "500"
threads = "1"

failed = os.listdir("failed/")
for i in range(0,len(failed)):
	failed[i] = failed[i].replace(".o","")

i = 0


with open("../e_coli_collections/FILTERED_MD_FINAL_ALL.tab", 'rU', encoding="ISO-8859-1") as f:
	for line in f:
		toks = line.strip().split("\t")
		if line.startswith("ID"):
			reads_loc = toks.index("Reads_Location")
			annot_loc = toks.index("New_annot_loc")
			continue
		if i == 0:
			i+=1
			continue

		reads_files = toks[reads_loc].split(",")
		gff_name = toks[annot_loc].replace(".gff","").split("/")[-1]
		if gff_name not in failed:
			continue
		## directory must be deleted for ariba to work!
		shutil.rmtree(gff_name)

		if len(reads_files) < 2: ## long read, can't be used with ariba
			continue
		if len(reads_files) == 2:
			if reads_files[0].endswith("_1.fastq.gz"):
				file1 = reads_files[0]
				file2 = reads_files[1]
			else:
				file1 = reads_files[1]
				file2 = reads_files[0]
		else:
			files = []
			for curr_file in reads_files:
				if gff_name in curr_file:
					files.append(curr_file)
			if len(files) >= 2:
				for curr_file in files:
					if curr_file.endswith("_1.fastq.gz"):
						file1 = curr_file
					else:
						file2 = curr_file
		commad_suffix = [file1, file2, gff_name]
		job_name = gff_name
		lsf_prefix = ["bsub", "-q", "normal", "-J", job_name, "-G", "team216", "-o", job_name + ".o",
		"-e", job_name + ".e", '-R"select[mem>' + mem + '] rusage[mem=' + mem + '] span[hosts=1]"', '-M' + mem, "-n" + threads]
		subprocess.call(lsf_prefix + command_prefix + commad_suffix)


