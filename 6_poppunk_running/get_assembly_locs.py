

mash = open("mash_clusters.csv", "w")
st = open("st_cluster.csv", "w")

mash.write("Taxon,Cluster\n")
st.write("Taxon,Cluster\n")

with open("/Users/gh11/e_colis/FINAL_METADATA_CLEANED.csv") as f:
	for line in f:
		toks = line.strip().split("\t")
		if line.startswith("ID"):
			assembly_index = toks.index("Assembly_Location")
			mash_index = toks.index("MASH")
			st_index = toks.index("ST")
			continue
		assemblies = toks[assembly_index].split(",")
		for a in assemblies:
			mash.write(a + "," + toks[mash_index] + "\n")
			st.write(a + "," + toks[st_index] + "\n")
			print(a)


mash.close()
st.close()