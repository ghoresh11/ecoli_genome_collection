

with open("../final_metadata_with_loc.csv") as f:
	for line in f:
		toks = line.strip().split("\t")
		if line.startswith("ID"):
			annot_loc = toks.index("Assembly_Location")
			continue
		annotations = toks[annot_loc].split(",")
		for a in annotations:
			print(a)