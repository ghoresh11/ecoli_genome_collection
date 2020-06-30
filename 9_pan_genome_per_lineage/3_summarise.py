
import argparse
import os
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import reverse_complement
from Bio.Seq import translate


def get_input_dirs(input_dir):
    ''' check all directories in the input dir
    return a list of directories that have all the required files'''
    print("Getting the input directories...")
    input_dir = os.path.abspath(input_dir)
    directories = [d for d in os.listdir(input_dir)]
    dirs_to_return = []
    for d in directories:
        d = os.path.join(input_dir, d)
        cluster = d.split("/")[-1].split("_")[0]
        if not os.path.isdir(d):
            continue
        if os.path.isfile(os.path.join(d, "summary_statistics.txt")):
            dirs_to_return.append(os.path.join(input_dir, d))
    return dirs_to_return

def classify_genes(d):
    ''' go over the Rtab file for each roary cluster,
    calculate the frequency of each gene in the cluster
    return: dict classification of genes into core, soft_core, intermediate, rare'''
    print("Calculating gene frequencies....")
    gene_class = {}
    with open(os.path.join(d, "gene_presence_absence.Rtab")) as f:
        for line in f:
            if line.startswith("Gene"):
                continue
            toks = line.strip().split("\t")
            freq = sum(map(int, toks[1:])) / float(len(toks)-1)
            gene_class[toks[0]] = freq
    output_genes_per_class(d, gene_class)
    return


def output_genes_per_class(d, gene_freqs):
    ''' use the dictionary, and the pan_genome_reference fasta file
    to create four different files for each cluster with the
    sequences of the genes in each category
    These will be used to postprocess the rare genes, and later
    to compare between different clusters.'''
    print("Getting the gene sequences....")
    counts = {}
    outputs = {}
    types = ["core", "soft_core", "inter", "rare"]
    for type in types:
        outputs[type] = open(os.path.join(d, type + "_genes.fa"), "w")
        counts[type] = 0
    with open(os.path.join(d, "pan_genome_reference.fa")) as handle:
        for values in SimpleFastaParser(handle):
            gene_name = values[0].split()[-1]
            protein_sequence = translate(values[1])
            if gene_freqs[gene_name] < 0.15:
                outputs["rare"].write(">" + values[0] + "\n" + protein_sequence + "\n")
                counts["rare"] += 1
            elif gene_freqs[gene_name] < 0.95:
                outputs["inter"].write(">" + values[0] + "\n" + protein_sequence + "\n")
                counts["inter"] += 1
            elif gene_freqs[gene_name] < 0.99:
                outputs["soft_core"].write(">" + values[0] + "\n" + protein_sequence + "\n")
                counts["soft_core"] += 1
            else:
                outputs["core"].write(">" + values[0] + "\n" + protein_sequence + "\n")
                counts["core"] += 1
    for o in outputs:
        outputs[o].close()
    return


def get_cluster_sizes():
    sizes = {}
    with open("cluster_sizes_updated.csv") as f:
        for line in f:
            if line.startswith("Cluster"):
                continue
            toks = line.strip().split(",")
            sizes[toks[0]] = toks[1]
    return sizes

def run(args):
    dirs = get_input_dirs(args.d)
    # for input_dir in dirs:
    #     classify_genes(input_dir)
    sizes = get_cluster_sizes()

    out = open("FINAL_summary_per_cluster.csv", "w")
    out.write("cluster,variable,count,size\n")

    for input_dir in dirs:
        cluster = input_dir.split("/")[-1].split("_")[0]
        with open(os.path.join(input_dir, "summary_statistics.txt")) as f:
            for line in f:
                toks = line.strip().split("\t")
                if toks[0] == "Total genes":
                    continue
                out.write(",".join([cluster,toks[0],toks[2],sizes[cluster]]) + "\n")
    out.close()
    return

def get_options():
    parser = argparse.ArgumentParser(description='Extract the gene sequences from roary outputs, and merge rare genes')
    # input options
    parser.add_argument('--d', required=False,
                        type=str, default =  "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/new_roary/",
                        help='path to inpqut directory [%(default)s]')
    return parser.parse_args()


if __name__ == "__main__":
    # get arguments from user
    options = get_options()
    run(options)
