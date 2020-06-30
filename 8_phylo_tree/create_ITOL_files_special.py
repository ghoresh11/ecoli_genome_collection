import seaborn as sns

genes_of_interest_file = "/Users/gh11/Documents/longer_26.txt"

with open(genes_of_interest_file) as f:
	for line in f:
		toks = line.strip().split(",")
		name = toks[0]
		poppunk_clusters_to_mark = []
		for c in toks[1].split():
			poppunk_clusters_to_mark.append(c.split("_")[0])
		out = open("genes/" + name + ".txt", "w")
		out.write("DATASET_BINARY\nSEPARATOR COMMA\nDATASET_LABEL," + name + "\nCOLOR,#ff0000\nFIELD_SHAPES,1\nFIELD_LABELS,f1\nDATA\n")


		with open("/Users/gh11/poppunk_pangenome/X_choose_reps/270619_chosen_treemer_corrected.csv") as f:
			for line in f:
				toks = line.strip().split(",")
				if line.startswith("Name"):
					phylogroup_index = toks.index("Phylogroup")
					assembly_index = toks.index("Annotation")
					cluster_index = toks.index("popppunk_cluster")
					continue

				name = toks[assembly_index].split("/")[-1].replace(".gff","")
				phylogroup = toks[phylogroup_index]
						
				cluster = toks[cluster_index]
				if cluster in poppunk_clusters_to_mark:
					out.write(name + ",1" + "\n")

		out.close()