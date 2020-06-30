import networkx as nx
from Bio.SeqIO.FastaIO import SimpleFastaParser
import os
import numpy as np


pairwise_dir = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/new_roary/pairwise/analysis/071019_mode_rep/"

G = nx.read_gml("wrong.gml")
# Step 4: use the connected components of the WRONG graph to correct
## any genes that should be marked as the same gene
cc = sorted(nx.connected_components(G), key=len, reverse=True) ##
gene_presence_absence = {}
with open(os.path.join(pairwise_dir,"complete_presence_absence.csv")) as f:
    for line in f:
        toks = line.strip().split(",")
        if line.startswith("Strain"):
            strains = line
            continue
        if line.startswith("Cluster"):
            clusters = line
            continue
        curr_counts = map(int,toks[1:])
        gene_presence_absence[toks[0]] = curr_counts
        ## gene gene_presence_absence has the old name

# This section is for writing the new presence absence file:
print("Generating new complete presence absence file...")
out = open("corrected_complete_presence_absence.csv","w")
out.write(strains + clusters)
strains = strains.split(",")[1:]
used_names = set()
old_gene_to_new_gene = {}
## here I can also apply a dbscan step-> similar to what I did for the panaroo outputs...
for members in cc:  # each members is one gene with all its members
    ## first check how many different clusters of each cluster this has
    curr_presence_absence = [0] * len(strains)
    names = []
    for m in members:
        names.append(m.replace("*",""))
        curr_presence_absence = [min(sum(x),1) for x in zip(curr_presence_absence, gene_presence_absence[m])] # merge with existing (nothing if the first time)
    gene_name = max(set(names), key=names.count) # get the most common gene name
    while gene_name in used_names: ## add stars to prevent duplicate names
        gene_name = gene_name + "*"
    used_names.add(gene_name)
    for m in members:
        old_gene_to_new_gene[m] = gene_name
    out.write(gene_name + "," + ",".join(list(map(str, curr_presence_absence))) + "\n")
out.close()

with open("old_gene_to_new_gene.csv", "w") as out:
    out.write("Prev_name,New_name\n")
    for m in old_gene_to_new_gene:
        out.write(m + "," + old_gene_to_new_gene[m] + "\n")

## create the new members.csv file by reading the old one and combining the lines
new_gene_to_all_members = {}
print("Generating new members file...")
with open(os.path.join(pairwise_dir, "members.csv")) as f:
    for line in f:
        if line.startswith("Gene"):
            continue
        toks = line.strip().split(",")
        old_gene = toks[0]
        new_gene = old_gene_to_new_gene[old_gene]
        if new_gene not in new_gene_to_all_members:
            new_gene_to_all_members[new_gene] = []
        new_gene_to_all_members[new_gene] += toks[1].split()

with open("corrected_members.csv","w") as out:
    out.write("Gene,Members\n")
    for gene in new_gene_to_all_members:
        out.write(gene + "," + "\t".join(new_gene_to_all_members[gene]) + "\n")

### Generate R outputs is taken directly from 2_build_pairwise_connections
''' reopen the complete presence absence file and write files
that can easily be plotted in R'''
print("Generating new freqs file...")
clusters_to_remove = [21,43,49]
melted_freqs = open("corrected_melted_gene_freqs.csv","w")
freqs = open("corrected_freqs.csv","w")
melted_freqs.write("Gene,Cluster,Freq,Class\n")
freqs.write("Gene")
for c in range(1,52):
    if c in clusters_to_remove:
        continue
    freqs.write("," + str(c))
freqs.write("\n")
with open("corrected_complete_presence_absence.csv") as f:
    for line in f:
        if line.startswith("Strain"):
            continue
        toks = line.strip().split(",")
        if line.startswith("Cluster"):
            indices = {}
            for cluster in range(1,52):
                if cluster in clusters_to_remove:
                    continue
                indices[cluster] = [i for i, x in enumerate(toks) if x == str(cluster)]
            continue
        gene_name = toks[0]
        curr_freqs = {}
        freqs.write(gene_name)
        for cluster in range(1,52):
            if cluster in clusters_to_remove:
                continue
            if len(indices[cluster]) == 0:
                freq = 0
            else:
                total = 0
                for i in indices[cluster]:
                    total += int(toks[i])
                freq = total/float(len(indices[cluster]))
            if freq < 0.15:
                gene_class = "rare"
            elif freq < 0.95:
                gene_class = "inter"
            elif freq < 0.99:
                gene_class = "soft_core"
            else:
                gene_class = "core"
            melted_freqs.write(",".join([gene_name, str(cluster), str(freq), gene_class]) + "\n")
            freqs.write(',' + str(freq))
        freqs.write("\n")
freqs.close()
melted_freqs.close()


## The last bit needs a lot of memory and I can't really understand why...
## bsub -J wrong -R"select[mem>20000] rusage[mem=20000]" -M20000  -G team216 -o wrong.o -e wrong.e python 3_correct_wrong_splits.py
