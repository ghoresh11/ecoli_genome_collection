## after converting between names and run IDs
## add that information to the large metadata file
## also mark genomes that need to be broken into artifical reads

## this will be inefficient, but i'd rather it be slow and correct

## 1. read in the conversions
conversions = {}
conversions_opposite = {}
with open("/Users/gh11/e_colis/genomes/READS/conversion_all.csv") as f:
	for line in f:
		if line.startswith("identifier"):
			continue
		toks = line.strip().split(",")
		if "?" in toks[2]:
			continue
		
		if not toks[0].lower() in conversions: ## mapping the names to the read IDs
			conversions[toks[0].lower()] = set()
		conversions[toks[0].lower()].add(toks[1])


		if not toks[1] in conversions_opposite: ## opposite mapping read IDs to their corresponding names
			conversions_opposite[toks[1]] = set()
		conversions_opposite[toks[1]].add(toks[0].lower())
		

## make a list of all the SRR ids that have reads available to them
available_reads = []  ## IDs in this list are available to download from the SRA

with open("/Users/gh11/e_colis/genomes/READS/all_run_ids.txt") as f:
	for line in f:
		if line.startswith("identifier"): ## header
			continue
		toks = line.strip().split()
		## toks[0] == the identifier I have as the read
		## toks[1] == the actual read ID, i need to swap them in the original "conversion" DF
		
		if toks[0] in conversions_opposite:

			for item in conversions_opposite[toks[0]]:
				if item not in conversions:
					conversions[item] = set()
				conversions[item].add(toks[1])

		available_reads.append(toks[1])

## rewrite metadata_fixed with two new columns
## first: SRR ID (multiple at times)
## second: whether it requires to create artifial reads

out = open("metadata_fixed_with_reads.csv", "w")


def get_run_ids(run_ids_set, available, all_unique_run_ids):
	out_vals = []
	for id_var in run_ids_set:
		if (id_var.startswith("SRR") or id_var.startswith("ERR")) and id_var in available:
			out_vals.append(id_var)
			all_unique_run_ids.add(id_var)
	if len(out_vals) == 0:
		return ("ND","Yes")
	return (",".join(out_vals), "No")

all_unique_run_ids = set()
SRR_to_line_nums = {}
lines = {}
line_num = 0
with open("metadata_fixed.csv") as f:
	for line in f:
		if line.startswith("ID"):
			line = line.strip()
			line += "\tRun_ID\tMake_artificial\n"
			out.write(line)
			continue

		toks = line.strip().split()
		name1 = toks[0]
		name2 = toks[1]

		line_num += 1
		run_id = "ND"
		artificial = "Yes"
		

		if name1.lower() in conversions:
			run_id, artificial = get_run_ids(conversions[name1], available_reads, all_unique_run_ids)
		if name2.lower() in conversions and run_id == "ND": ## couldn't find runs using first name
			run_id, artificial = get_run_ids(conversions[name2], available_reads, all_unique_run_ids)

		## sanger sequences have "#" in them (other than NCTC files, but there is a problem with them!)
		if "#" in name1 or toks[-3] == "nctc":
			run_id, artificial = name1, "No"
		elif "#" in name2:
			run_id, artificial = name2, "No"
		
		srrs = run_id.split(",")
		for srr in srrs:
			if srr not in SRR_to_line_nums:
				SRR_to_line_nums[srr] = []
			SRR_to_line_nums[srr].append(line_num)
		
		## instead of "write" ass to lines
		lines[line_num] = line.strip() + "\t" + run_id + "\t" + artificial + "\n"



### over here -> change it to merge lines that have the same SRR ids, this means that they are in fact the same sequence 
## and will significantly reduce the size of my table
written = set() ## keep track of the lines that had been written

for srr in SRR_to_line_nums:

	if len(SRR_to_line_nums[srr]) == 1 or srr == "ND":
		for i in SRR_to_line_nums[srr]:
			if not i in written:
				out.write(lines[i]) ## nothing to choose from
				written.add(i)
		continue

	chosen = SRR_to_line_nums[srr][0]
	chosen_line = 0
	for i in SRR_to_line_nums[srr]:
	 	if not lines[i].startswith("srs"): ## choose whatever isn't SRS at random
			chosen = lines[i]
			chosen_line = i


	if not chosen_line in written:
		written.add(chosen_line)
		out.write(chosen)


out.close()


## regardless write to a file all the SRR files that need to be downlaoded
with open("READS/check_existence.txt","w") as out:
	for id_var in all_unique_run_ids:
		out.write(id_var + "\n")

