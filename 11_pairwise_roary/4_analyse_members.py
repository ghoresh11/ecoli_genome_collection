from numpy import mean
from sys import 
'''
- Draw a distribution of how many gene cluster have how many from each genomic group, see if some groups are causing problems
- count how many times a gene has been merged with a another gene from the same roary output (i.e. the same cluster appears twice in a row)
-> i.e. these connection of genes that weren't merged before and suddenly they are - that means they were marked in the original roary as different but now are treated as the same
- When that happens, have a look at the gene properties table and see that these genes are usually like (length, position etc.)
I can extract the networks from the edges and nodes file for clusters that are interesting
'''
clusters = range(1,21) + range(22,40)
out = open(argv[2], 'w')
out.write("Gene, " + ",".join(map(str,clusters)) + ", mean" + "\n")

with open(argv[1]) as f:
    for line in f:
        if line.startswith("Name"):
            continue
        toks = line.strip().split(",")
        name = toks[0]
        members = toks[1].split("\t")
        counts = {}
        flag = False
        for m in members:
            curr_cluster = int(m.split("_")[0])
            if curr_cluster not in counts:
                counts[curr_cluster] = 1
            else: ## this cluster was seen already for this gene
                counts[curr_cluster] += 1
                flag = True

        if not flag:
            continue

        out.write(name)
        for cluster in clusters:
            if cluster not in counts:
                out.write(",0")
            else:
                out.write("," + str(counts[cluster]) )
        out.write("," + str(mean(counts.values())) + "\n")
out.close()
