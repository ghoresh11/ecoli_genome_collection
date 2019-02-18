import os

def get_val_from_line(line):
	val = line.strip().split(" ")[-1]
	val = val.replace("(","")
	val = val.replace(").","")
	val = val.replace(")","")
	val = val.replace("%","")
	val = val.split("/")
	if len(val) == 2:
		if val[1] == "0":
			val = "0"
		else:
			val = str(float(val[0]) / float(val[1]))
	else:
		val = str(val[0])

	return val

def get_qc_vals(f, run_id):
	
	d = {}
	header = ["kraken", "mapped_bases", "error_rate", "indel_ratio", "insert_size", "insert_size_rev", "heterozygous_snps"]
	for h in header:
		d[h] = "NA"

	kraken_file = f
	report_file = "/".join(kraken_file.split("/")[:-1] + ["qc-sample/auto_qc.txt"])
	hg_snps_file = "/".join(kraken_file.split("/")[:-1] + ["qc-sample/" + run_id + "_heterozygous_snps_report.txt"])

	### get the kraken report output
	with open(kraken_file) as kf:
		for line in kf:
			if line.startswith("#"):
				continue
			toks = line.strip().split("\t")
			if toks[3] == "S":
				if toks[-1].strip() != "Escherichia coli":
					d["kraken"] = "NA"
				else: 
					d["kraken"] = toks[0]
				break
	### get the output from the auto report
	with open(report_file) as rf:
		for line in rf:
			if line.startswith("Mapped"):
				d["mapped_bases"] = get_val_from_line(line)
				continue
			if line.startswith("Error"):
				d["error_rate"] = get_val_from_line(line)
				continue
			if line.startswith("InDel"):
				d["indel_ratio"] = get_val_from_line(line)
				continue
			if line.startswith("Insert size:"):
				d["insert_size"] = get_val_from_line(line)
				continue
			if line.startswith("Insert size (rev)"):
				d["insert_size_rev"] = get_val_from_line(line)
				continue

	with open(hg_snps_file) as hgf:
		for line in hgf:
			if line.startswith("No."):
				continue
			d["heterozygous_snps"] = line.strip().split("\t")[-1]
	return d


if __name__ == '__main__':
	
	qc_report_files = os.listdir(".")
	all_runs = {}
	for qc_file in qc_report_files:
		
		if not qc_file.endswith("_qc.txt"):
			continue

		print("Obtaining QC values for %s..." %qc_file)

		with open(qc_file) as f:
			for line in f:
				if not line.startswith("/"): # not a location line
					continue
				run_id = line.split("/")[-2]
				run_id = run_id.replace("_trimmed","")
				run_id = run_id.replace("_gal","")
				all_runs[run_id] = get_qc_vals(line.strip(), line.split("/")[-2])

	print("Generating output file...")
	header = ["kraken", "mapped_bases", "error_rate", "indel_ratio", "insert_size", "insert_size_rev", "heterozygous_snps"]
	with open("report.csv", "w") as out:
		out.write(",".join(["ID"] + header) + "\n")
		for run_id in all_runs:
			out.write(run_id)
			for h in header:
				out.write(", ")
				out.write(all_runs[run_id][h])
			out.write("\n")




