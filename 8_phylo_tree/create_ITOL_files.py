import seaborn as sns

out_clusters = open("cluster_labels.txt", "w")
out_clusters.write("DATASET_SYMBOL\nSEPARATOR COMMA\nDATASET_LABEL,cluster\nCOLOR,#ff0000\nDATA\n")

out_labels = open("labels.txt", "w")
out_labels.write("DATASET_TEXT\nSEPARATOR COMMA\nDATASET_LABEL,label\nCOLOR,#ff0000\nDATA\n")


out_phylogroup = open("phylogroup_labels.txt", "w")
out_phylogroup.write("DATASET_COLORSTRIP\nSEPARATOR COMMA\nDATASET_LABEL,phylogroup\nCOLOR,#ff0000\nDATA\n")


cols = sns.color_palette("Set2").as_hex()
cols = cols + [cols[7], cols[7]]
phylogroup_cols = {}
index = 0
for p in ["B2", "B1", "A", "D", "E", "C", "F", "Shigella", "U"]:
	phylogroup_cols[p] = cols[index]
	index += 1

cluster_cols = {}
shapes = {"15": [1, 1], "12": [1,0], "7": [5,1],
			"2": [4,1], "17": [4,0],
			"16": [2,1], "10": [2,0],
			"9": [2,1], "5": [2,0],
			"8": [5,0], "4": [2,0] }

with open("/Users/gh11/Submissions/my_thesis/Chapter4/figures/cluster_graphics.csv") as f:
	for line in f:
		if line.startswith("Cluster"):
			continue
		toks = line.strip().split(",")
		cluster_cols[toks[0]] = {"color": toks[2], "shape":shapes[toks[3]]}


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
		
		out_phylogroup.write(",".join([name, phylogroup_cols[phylogroup], phylogroup]) + "\n")
		
		cluster = toks[cluster_index]
		curr_cluster_color = cluster_cols[cluster]["color"]
		curr_cluster_shape = cluster_cols[cluster]["shape"]		
		curr_out_cluster = [name, curr_cluster_shape[0], 3, curr_cluster_color, curr_cluster_shape[1], 1]
		out_clusters.write(",".join(map(str, curr_out_cluster)) + "\n")

		out_labels.write(",".join([name, cluster, "-1", curr_cluster_color, "bold", "3", "0"]) + "\n")

out_phylogroup.close()
out_clusters.close()
out_labels.close()