import sys

names = []
# with open("/Users/gh11/poppunk_pangenome/4_pairwise_roary/071019_mode_rep/freqs.csv") as f:
# 	for line in f:
# 		toks = line.strip().split(",")
# 		if toks[0].startswith("yehB"):
# 			names.append(toks[0])

## to make the barcharts
with open("/Users/gh11/poppunk_pangenome/5_classify_genes/interpro_scan_results/Pertactin.csv") as f:
	for line in f:
		toks = line.strip().split("\t")
		names.append(toks[0].split(":")[1])


for name in names:
	## step one -> get the frequency of the gene in each cluster
	freq_per_cluster = {}
	with open("/Users/gh11/poppunk_pangenome/4_pairwise_roary/071019_mode_rep/freqs.csv") as f:
		for line in f:
			if line.startswith("Gene"):
				clusters = line.strip().split(",")[1:]
				continue
			toks = line.strip().split(",")
			if toks[0] != name:
				continue
			freqs = toks[1:]
			for i in range(0,len(freqs)):
				freq_per_cluster[clusters[i]] = freqs[i]
			break

	## step 2 -> create the ITOL file
	out = open("barcharts/" + name + ".txt", "w")
	out.write("DATASET_MULTIBAR\nSEPARATOR COMMA\nDATASET_LABEL," + name + "_bar\nCOLOR,#ff0000\nFIELD_COLORS,#ff0000\nFIELD_LABELS,f1\nWIDTH,50\nDATASET_SCALE,1\nDATA\n")

	with open("/Users/gh11/poppunk_pangenome/X_choose_reps/270619_chosen_treemer.csv") as f:
		for line in f:
			toks = line.strip().split(",")
			if line.startswith("Name"):
				phylogroup_index = toks.index("Phylogroup")
				assembly_index = toks.index("Annotation")
				cluster_index = toks.index("popppunk_cluster")
				continue

			curr_genome = toks[assembly_index].split("/")[-1].replace(".gff","")
			phylogroup = toks[phylogroup_index]
					
			cluster = toks[cluster_index]
			if cluster not in freq_per_cluster:
				continue
			if curr_genome == "SRR7187888_trimmed":
				continue
			out.write(curr_genome + "," + freq_per_cluster[cluster] + "\n")
	out.close()