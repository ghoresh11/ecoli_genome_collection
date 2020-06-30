import argparse
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import translate
import subprocess
import os


''' given that some genes stand out as having many members, run this
analysis on these genes which will enable to open them in cytoscape,
and also to view a multiple sequence alignment from all the members'''

def get_options():
    parser = argparse.ArgumentParser(description='Get the subnetwork and edges and nodes relevant to one bigger gene cluster')
    ## for running on debug mode
    ## python 5_get_network_for_gene.py --gene group_11315 --members_file members.csv --edges_file edges.txt --nodes_file nodes.txt
    # input options
    parser.add_argument('--gene', required=True,
                        type=str,
                        help='name of gene to retrieve')
    parser.add_argument('--members_file', required=False,
                        type=str, default = "030519_members.csv",
                        help='name of members file')
    parser.add_argument('--edges_file', required=False,
                        type=str, default = "030519_edges.txt",
                        help='name of edges file')
    parser.add_argument('--nodes_file', required=False,
                        type=str, default = "030519_nodes.txt",
                        help='name of nodes file')

    return parser.parse_args()


def get_members_for_gene(gene, members_file):
    with open(members_file) as f:
        for line in f:
            toks = line.strip().split(",")
            if toks[0] == gene:
                return toks[1].split("\t")
    print("Gene %s not found!!" %gene)
    quit()
    return

def define_members_per_cluster(member):
    cluster_members = {}
    for m in members:
        cluster = m.split("_")[0]
        name = "_".join(m.split("_")[1:])
        if cluster not in cluster_members:
            cluster_members[cluster] = []
        cluster_members[cluster].append(name)
    return cluster_members

def create_fasta_of_members(cluster_members, gene, outdir):
    for curr_cluster in cluster_members:

        curr_fasta = os.path.join(outdir, curr_cluster + "_" + gene + "_members.fa")
        out = open(curr_fasta, "w")
        ref_dir = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk//dists_analysis/roary_outputs/pan_genome_references/"
        # ref_dir = "/Users/gh11/poppunk_pangenome/2_dists_roary_analysis/pan_genome_references/"
        with open(os.path.join(ref_dir ,curr_cluster + "_reference.fa")) as handle:
            for values in SimpleFastaParser(handle):
                if curr_cluster == "1":
                    name = values[0].split()[0]
                else:
                    name = values[0].split()[1]
                if name in cluster_members[curr_cluster]:
                    out.write(">" + curr_cluster + "_" + name + "\n" + values[1] + "\n")
        out.close()
        # create a multiple sequence alignment of current file
        out_msa = open(os.path.join(outdir,  curr_cluster + "_" + gene + "_msa.fa"), "w")
        subprocess.Popen(["mafft", "--leavegappyregion", curr_fasta], stdout=out_msa, stderr=subprocess.PIPE)
        out.close()
    return


def filter_file(filename, gene, members, filetype, outdir):
    out = open(os.path.join(outdir, gene + "_" + filetype + ".txt"), "w")
    with open(filename) as f:
        for line in f:
            if line.startswith("Name"):
                out.write(line)
                continue
            toks = line.strip().split("\t")
            if toks[0] in members:
                out.write(line)
    out.close()
    return

if __name__ == "__main__":
    options = get_options()

    outdir = os.path.join("specific_genes", options.gene)
    try:
        os.makedirs(outdir)
    except OSError as e:
        pass

    members = get_members_for_gene(options.gene, options.members_file)
    cluster_members = define_members_per_cluster(members)
    create_fasta_of_members(cluster_members, options.gene, outdir)
    #filter_file(options.nodes_file, options.gene, members, "nodes", outdir)
    #filter_file(options.edges_file, options.gene, members, "edges", outdir)


## THe best thing to do at this point is probably an MSA of all these genes...
