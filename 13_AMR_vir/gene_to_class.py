
d = {}
with open("resfinder_genes_by_class.txt") as f:
    for line in f:
        if line.startswith("#"):
            curr_class = line.strip().replace(":","").replace("#","")

            continue
        curr_gene = line.strip().split(":")[0]
        for c in ["'","(",")","-"]:
            curr_gene = curr_gene.replace(c,"_").lower()
        curr_gene = curr_gene.split(",")
        for g in curr_gene:
            if g not in d:
                d[g] = []
            d[g].append(curr_class)

out = open("Gene_to_class.csv", "w")
out.write("Gene,Altered_name,Class\n")
with open("/Users/gh11/poppunk_pangenome/7_AMR_vir_plasmid/results/resfinder.csv") as f:
    for line in f:
        toks = line.strip().split(",")
        if line.startswith("Genome"):
            genes = toks[1:]
            for g in genes:
                name = "_".join(g.split("_")[:-1]).lower()
                if name not in d:
                    name = "_".join(name.split("_")[:-1]).lower()
                if name not in d:
                    out.write(g + "," + name + ",?\n")
                else:
                    out.write(g + ","+ name +","+ "/".join(d[name]) + "\n")

out.close()
