import networkx as nx
from Bio.SeqIO.FastaIO import SimpleFastaParser
import os

''' build a graph where the relationship is "truncate" or wrong -> then I can
colour the edges by the relationship. Dont draw CCs of size 1
For each gene also save the std value and length of that node'''

print("Initiating graph...")
Gs = {"wrong":nx.Graph(), "trunc":nx.Graph()}
## step 1: build a graph from the combined pan genomes reference
with open("combined_pan_genome_reference_prot.fa") as handle:
    for values in SimpleFastaParser(handle):
        length = len(values[1])
        sd = values[0].split()[-1].replace("-","")
        for key in Gs:
            Gs[key].add_node(values[0].split()[0], length = length, sd = sd, sequence = values[1])

print("Reading blast results...")
## step 2: add edges from blast results
with open("blast_results_2") as f:
    for line in f:
        toks = line.strip().split()
        if toks[0] == toks[1]: ## ignore self loops
            continue
        ## definition: if they differ in length by less than 80% and the %identity is over 95 they should be the same
        ## otherwise: they're truncated
        coverage = max(int(toks[3])/float(toks[4]),int(toks[3])/float(toks[5]))
        length_ratio  = min( float(toks[5])/float(toks[4]), float(toks[4])/float(toks[5]))
        if coverage < 0.95 or float(toks[2]) < 95: ## always look for 95% identity and almost full coverage of shorter sequence
            continue

        if length_ratio > 0.8: ## that's how I originally defined it
            Gs["wrong"].add_edge(toks[0], toks[1], identity = toks[2], alignment = toks[3], ratio = length_ratio)
        else:
            ## check where the truncation is
            protein_lengths = list(map(int, [toks[4], toks[5]]))
            longer_prot_index = protein_lengths.index(max(protein_lengths))
            longer_length = protein_lengths[longer_prot_index]
            start = int(toks[7 + (longer_prot_index*2)])
            stop = int(toks[8 + (longer_prot_index*2)])
            if start < 0.1 * longer_length:
                loc = "N-Terminus"
            elif stop < 0.9 * longer_length:
                loc = "Middle"
            else:
                loc = "C-Terminus"
            Gs["trunc"].add_edge(toks[0], toks[1], identity = toks[2], alignment = toks[3], location = loc, ratio = length_ratio, start = start, stop = stop)

print("Removing singletons...")
## step 3: remove nodes of size 1

for key in Gs:
    # singletons = list(nx.isolates(Gs[key]))
    # with open("singletons.csv","w") as out:
    #     out.write("Name,Length,SD\n")
    #     for n in singletons:
    #         out.write(",".join([n, str(Gs[key].nodes()[n]["length"]), Gs[key].nodes()[n]["sd"]]) + "\n")
    # Gs[key].remove_nodes_from(singletons)
    print("Generating output..")
    nx.write_gml(Gs[key], key + ".gml")
