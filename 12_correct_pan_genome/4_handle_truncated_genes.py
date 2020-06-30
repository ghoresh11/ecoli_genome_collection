import networkx as nx
import os
from Bio.SeqIO.FastaIO import SimpleFastaParser
import subprocess

## read the truncated gene file and get the connected components
## of that file to generate a CSV file of all the genes that are truncated
## versions of other genes
## TODO: Genetic context to confirm? -> MAYBE, would add a lot but would be a lot more work

pairwise_dir = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/new_roary/pairwise/analysis/071019_mode_rep/"

clusters = []
for i in range(1,52):
    if i not in [21,43,49,50]:
        clusters.append(str(i))

print("Getting counts per gene...")
gene_to_counts = {}
with open("corrected_members.csv") as f:
    for line in f:
        if line.startswith("Gene"):
            continue
        toks = line.strip().split(",")
        curr_counts = [0] * len(clusters)
        for m in toks[1].split():
            curr_cluster = m.split("_")[0]
            curr_counts[clusters.index(curr_cluster)] += 1
        gene_to_counts[toks[0]] = curr_counts

print("Getting new members...") ## get the corrected_members.csv file and group ones that are the same
G_wrong = nx.read_gml("wrong.gml")
# Step 4: use the connected components of the WRONG graph to correct
## any genes that should be marked as the same gene
cc_wrong = sorted(nx.connected_components(G_wrong), key=len, reverse=True) ##
member_to_rep = {}
#rep_to_test = "?"
for cc in cc_wrong:  # each members is one gene with all its members
    members = list(cc)
    rep = members[0]
    for m in members:
        member_to_rep[m] = rep

old_gene_to_new_gene = {}
with open("old_gene_to_new_gene.csv") as f:
    for line in f:
        toks = line.strip().split(",")
        old_gene_to_new_gene[toks[0]] = toks[1]

print("Getting components from graph...")
G = nx.read_gml("trunc.gml")
ccs = sorted(nx.connected_components(G), key=len, reverse=True)
cluster_count = 0 ## id for the connected components
visited = [] ## keep track of visited nodes

out2 = open("trunc_relationships.csv", "w")
out2.write("ID,GeneA,GeneB,LengthA,LengthB,StdA,StdB,Alignment,Location,Ratio,Start,Stop\n")
with open("trunc_components.csv","w") as out:
    out.write("ID,Gene,Length,Std,Num_clusters," + ",".join(clusters) + "\n")
    for cc in ccs:
        if len(cc) == 1:
            break
        cluster_count += 1
        print("On connected component number %d..." %cluster_count)
        nodes_as_list = list(cc)
        for node in nodes_as_list:
            if member_to_rep[node] != node: ## duplicate and was removed
                continue
            counts = gene_to_counts[old_gene_to_new_gene[node]]
            num_clusters = len([x for x in counts if x > 0])
            out.write(",".join(list(map(str,[cluster_count, old_gene_to_new_gene[node], G.nodes()[node]["length"], G.nodes()[node]["sd"], num_clusters] + counts))) + "\n")
        for i in range(0,len(nodes_as_list) - 1):
            for j in range(i+1, len(nodes_as_list)):
                nodeA = nodes_as_list[i]
                nodeB = nodes_as_list[j]

                if G.nodes()[nodeA]["length"] < G.nodes()[nodeB]["length"]:
                    tmp = nodeA
                    nodeA = nodeB
                    nodeB = tmp
                if not G.has_edge(nodeA, nodeB):
                    continue

                if member_to_rep[nodeA] != nodeA or member_to_rep[nodeB] != nodeB: ## they're duplicated and they have a different rep
                    continue

                out2.write(",".join(list(map(str,[cluster_count, old_gene_to_new_gene[nodeA], old_gene_to_new_gene[nodeB],
                G.nodes()[nodeA]["length"], G.nodes()[nodeB]["length"],
                G.nodes()[nodeA]["sd"], G.nodes()[nodeB]["sd"],
                G[nodeA][nodeB]["alignment"], G[nodeA][nodeB]["location"],
                 G[nodeA][nodeB]["ratio"], G[nodeA][nodeB]["start"], G[nodeA][nodeB]["stop"]]))) + "\n")
out2.close()
