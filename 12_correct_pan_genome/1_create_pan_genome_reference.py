import os
from Bio.SeqIO.FastaIO import SimpleFastaParser
from numpy import std
from Bio.Seq import translate

''' create the pan-genome reference file for the combined genomes
Because I had made the correction with the SNP distances I think
it does make sense to take the LONGEST rep within each cluster,
but I will also compare the difference in length between the
longest and the mode within each cluster'''


combine_dir = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/new_roary/pairwise/analysis/071019_mode_rep/"
roary_dir = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/new_roary/"

## 1. Read the members.csv file to get the member from each popppunk cluster1
print("Reading members file...")
cluster_to_rep_to_gene = {}
with open(os.path.join(combine_dir, "members.csv")) as f:
    for line in f:
        if line.startswith("Gene"):
            continue
        toks = line.strip().split(",")
        for rep in toks[1].split():
            cluster = rep.split("_")[0]
            rep = "_".join(rep.split("_")[1:])
            if cluster not in cluster_to_rep_to_gene:
                cluster_to_rep_to_gene[cluster] = {}
            cluster_to_rep_to_gene[cluster][rep] = toks[0]

## 2. Read the pan-genome-reference of each poppunk cluster and save the rep for each gene
gene_to_rep_to_sequence = {}
print("Reading pan genome files...")
for cluster in cluster_to_rep_to_gene:
    print("Cluster: %s..." %cluster)
    with open(os.path.join(roary_dir, cluster, "mode_pan_genome_reference.fa")) as handle:
        for values in SimpleFastaParser(handle):
            rep = values[0].split()[-1]
            if rep not in cluster_to_rep_to_gene[cluster]:
                continue
            curr_gene = cluster_to_rep_to_gene[cluster][rep]
            if curr_gene not in gene_to_rep_to_sequence:
                gene_to_rep_to_sequence[curr_gene] = {}
            gene_to_rep_to_sequence[curr_gene][cluster + "_" + rep] = values[1]

print(len(gene_to_rep_to_sequence))
## 3. Choose either the longest or the mode and generate output
out2 = open("mode_v_longest.csv","w")
out2.write("Gene,Mode,Longest\n")
out3 = open("combined_pan_genome_reference_prot.fa","w")
print("Generating output per gene...")
cnt = 0
with open("combined_pan_genome_reference.fa","w") as out:
    for gene in gene_to_rep_to_sequence:
        cnt += 1
        print("Gene num %d" %cnt)
        curr_reps = gene_to_rep_to_sequence[gene]
        lengths = list(map(len, curr_reps.values()))
        sd = str(std(lengths))
        max_length = max(lengths)
        try:
            mode_length = mode(lengths)
        except:
            ## this is a problem actually.. this means that roary clustered these badly:
            lengths.sort()
            mid_value = int(len(lengths)/2)-1
            mode_length = lengths[mid_value]

        out2.write(",".join([gene, str(mode_length), str(max_length)]) + "\n")
        for rep in curr_reps:
            if len(curr_reps[rep]) == mode_length:
                out.write(">" + gene + "\t" + rep + "\t" + sd + "\n" + curr_reps[rep] + "\n")
                out3.write(">" + gene + "\t" + rep + "\t" + sd + "\n" + translate(curr_reps[rep]) + "\n")
                break
out2.close()
out3.close()


### after this run blast (only used 3GB and took 4 minutes)
## makeblastdb -in combined_pan_genome_reference_prot.fa -dbtype prot
## runjob5 blastp -db combined_pan_genome_reference_prot.fa -query combined_pan_genome_reference_prot.fa -out blast_results_2 -outfmt  "6 qseqid sseqid pident length qlen slen evalue qstart qend sstart send"  -evalue 0.01 -num_threads 4
