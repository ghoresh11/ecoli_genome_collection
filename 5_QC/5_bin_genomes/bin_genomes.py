#!/usr/bin/env python2

import os
import sys
import subprocess
import glob
import shutil
from collections import Counter
import networkx as nx
import argparse
import random
import string
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import translate
from Bio.Seq import reverse_complement


class Error (Exception):
    pass


def qc_fasta(fasta_file, qc, out, name=""):
    ''' go over FASTA file and see if it meets QC requirements
    Output: LOG of number of contigs and length.
    Return: True, if it does
    False: If it doesn't, and remove the file '''
    if name == "":
        name = fasta_file

    num_contigs = 0
    length = 0
    with open(fasta_file) as handle:
        for values in SimpleFastaParser(handle):
            num_contigs += 1
            length += len(values[1].strip())
    length = length / 1000000.0
    if qc["max_contigs"] < num_contigs or qc["min_length"] > length or qc["max_length"] < length:
        out.write(name + "," + str(num_contigs) + "," + str(length) + "\n")
        return False
    return True


def qc_genomes(job_id, genomes_file, out_dir, args, gff, verbose):
    ''' Go over genomes in genome file
    Calc number of contigs and length, remove if they don't pass quality control
    If in GFF format = convert to FASTA
    If in FASTA format = create a symlink
    Genomes which don't pass QC criteria are removed
    Return: a mapping between old location and new location
    Output: genomes which were removed in QC'''

    tmp_out = os.path.join(out_dir, "fasta_files")
    log = open(job_id + "_removed_in_qc.csv", "w")
    log.write("File, Num_contigs, Length (MBP)\n")

    try:
        os.mkdir(tmp_out)
    except OSError:
        sys.stderr.write(
            "Output directory for this job already exists (%s). Please remove directory\n" % tmp_out)

    if gff:
        new_loc_to_old_loc = _convert_gffs_to_fasta(genomes_file, tmp_out,
                                                    args, log, verbose)
    else:
        new_loc_to_old_loc = _symlink_fastas(
            genomes_file, tmp_out, args, log, verbose)
    log.close()
    return new_loc_to_old_loc


def gff_to_fasta(gff_file, fasta_file, protein_coding=False, qc=None, log=None):
    '''Convert a gff file with the appended FASTA to protein/all fasta file
    gff_file = input gff file
    fasta_file = output file
    output: a protein coding FASTA file OR nucleotide FASTA file'''
    out_tmp = ''.join(random.choice(string.ascii_lowercase +
                                    string.ascii_uppercase + string.digits) for _ in range(7)) + "_fasta.fa"
    out = open(out_tmp, "w")
    contigs = {}
    with open(gff_file) as f:
        fasta = False
        for line in f:
            if fasta:
                out.write(line)
                continue
            if line.startswith("##FASTA"):
                fasta = True
                continue
            if line.startswith("#"):
                continue
            toks = line.strip().split()
            if toks[2] != "CDS":
                continue

            name = toks[-1].split("|")[-1]
            if toks[0] not in contigs:
                contigs[toks[0]] = []

            contigs[toks[0]].append({"name": name, "start": int(toks[3]) - 1,
                                     "stop": int(toks[4]), "strand": toks[6]})
    out.close()

    # not protein coding, in this case we apply QC measures
    if not protein_coding:
        if qc_fasta(out_tmp, qc, log, name=gff_file):  # if passed quality control
            # rename the temp to the final and return
            os.rename(out_tmp, fasta_file)
            return fasta_file
        else:
            os.remove(out_tmp)  # remove the temp file if didnt pass QC
            return None

    # read the contigs and save the final fasta file
    out = open(fasta_file, "w")
    with open(out_tmp) as handle:
        for values in SimpleFastaParser(handle):
            curr_contig = values[0]
            if curr_contig not in contigs:  # no CDSs in this contig
                continue
            for cds in contigs[curr_contig]:
                out.write(">" + cds["name"] + "\n")
                seq = values[1][cds["start"]:cds["stop"]]
                if cds["strand"] == "-":
                    seq = reverse_complement(seq)
                out.write(translate(seq) + "\n")
    out.close()
    os.remove(out_tmp)
    return



def _convert_gffs_to_fasta(genomes_file, tmp_out, args, log, verbose):
    ''' Convert a GFF files into FASTA file format so that they can be used for specific methods
    like MASH, orthofinder etc.
    Genomes which don't pass QC criteria are removed
    Return: a mapping between the GFF location and FASTA location.'''
    fasta_loc_to_gff_loc = {}
    if verbose:
        print("Converting the GFF files to FASTA format and QC")

    with open(genomes_file) as f:
        for line in f:
            gff_file = line.strip()
            new_loc = os.path.join(
                tmp_out, os.path.basename(gff_file) + ".fa")

            if not args["debug"]:
                gff_to_fasta(gff_file, new_loc, qc=args, log=log)
            fasta_loc_to_gff_loc[new_loc] = gff_file

    return fasta_loc_to_gff_loc


def _symlink_fastas(genomes_file, out_dir, args, log, verbose):
    ''' if genomes are provided in FASTA format, simply create
    a symlink to a single directory
    Genomes which don't pass QC criteria are removed
    Return: a mapping between the original location and new location.'''
    new_loc_to_old_loc = {}
    with open(genomes_file) as f:
        for line in f:
            fasta_file = os.path.abspath(line.strip())
            # create symilnk and add to dict
            new_loc = os.path.join(
                out_dir, os.path.basename(fasta_file) + ".fa")

            new_loc_to_old_loc[new_loc] = fasta_file

            if not args["debug"] and qc_fasta(fasta_file, args, log):
                if not os.path.isfile(new_loc):
                    os.symlink(fasta_file, new_loc)

    return new_loc_to_old_loc


def run_binning_method(out_dir, fasta_dir, method, cpu, verbose):
    ''' Depending on the binning method used (MASH, popPunk, MLST),
    run the specified method to retrieve an output that can be read for binning the genomes'''

    if verbose:
        print("Running the binning method %s" % method)
    if method.lower() == "mash":
        _run_mash(out_dir, fasta_dir, cpu)
    else:  # if adding any new method, add postprocessing step to match the MASH outputs
        sys.stderr.write(
            "binning method %s not supported. Supported methods: [mash]" % method)
        sys.exit()
    return


def _run_mash(out_dir, fasta_dir, cpu):
    ''' Run MASH on all the genomes in fasta_dir
    Output: MASH output of all against all MASH distances
    Return: Nothing '''
    # Step 1: sketch
    files = glob.glob(os.path.join(fasta_dir, "*.fa"))
    subprocess.call(["mash", "sketch", "-o", os.path.join(out_dir, "reference"),
                     "-p", str(cpu)] + files)

    # Step 2: dist
    out = subprocess.check_output(["mash", "dist", "-p", str(cpu),
                                   os.path.join(out_dir, "reference.msh")] + files)

    with open(os.path.join(out_dir, "distances.tab"), "w") as f:
        f.write(out)
    return


def _bin_nodes(line, G, threshold):
    ''' read a line from the MASH output and bin the nodes accordingly
    to threshold (species and cluster) '''
    toks = line.strip().split("\t")
    node1, node2, dist = toks[0], toks[1], float(toks[2])
    G.add_node(node1)
    G.add_node(node2)

    if node1 == node2:
        return
    if dist < threshold:
        G.add_edge(node1, node2)
    return


def bin(job_id, out_dir, distance, species_threshold, fasta_loc_to_gff_loc, verbose):
    ''' Using the output of "method", bin the genomes according to the threshold given.
    Ouput: A file with genomes and bin ID in the second column that can be loaded into "randomise"
               A file with all the genomes which are not the same species
    Return: Nothing '''
    species_network = nx.Graph()
    cluster_network = nx.Graph()

    if verbose:
        print("Reading distance file...")

    with open(os.path.join(out_dir, "distances.tab")) as f:
        for line in f:
            species_bins_id = _bin_nodes(
                line, species_network, species_threshold)  # pick apart contaminents
            # bin into clusters according to threshold
            clusters_bins_id = _bin_nodes(
                line, cluster_network, distance)

    # get the largest connected component
    species_cc = sorted(nx.connected_components(
        species_network), key=len, reverse=True)[0]
    clusters_ccs = sorted(nx.connected_components(
        cluster_network), key=len, reverse=True)

    # write the output file, ignore isolates which don't
    # have a species cluster of "species_cluster"
    if verbose:
        print("Creating output files...")
    bin_id = 1
    with open(job_id + "_binned.txt", "w") as out:
        with open(job_id + "_contaminents.txt", "w") as out2:
            for c in clusters_ccs:
                for genome in c:
                    if genome not in species_cc:  # contaminent
                        out2.write(fasta_loc_to_gff_loc[genome] + "\n")
                        continue
                    out.write(fasta_loc_to_gff_loc[genome] + ","
                              + str(bin_id) + "\n")
                bin_id += 1
    return


def remove_temp(job_id, keep_temp, out_dir, new_loc_to_old_loc, verbose):
    ''' Remove the temporary files generated in the binning process'''
    if verbose:
        print("Removing temporary files...")

    ## rewrite the distance file
    with open(job_id + "_distances.tab", "w") as out:
        with open(os.path.join(out_dir, "distances.tab")) as f:
            for line in f:
                toks = line.split("\t")
                toks[0] = new_loc_to_old_loc[toks[0]]
                toks[1] = new_loc_to_old_loc[toks[1]]
                out.write("\t".join(toks))

    if keep_temp:
        return

    shutil.rmtree(out_dir)
    return


def run(args):

    out_dir = args["job_id"] + "_bin"
    try:
        os.mkdir(out_dir)
    except OSError:
        sys.stderr.write(
            "Output directory for this job already exists (%s). Please remove directory\n" % out_dir)
        # sys.exit()

    new_loc_to_old_loc = qc_genomes(args["job_id"], args["genomes_file"], out_dir, args, args["gff"],
                                    args["verbose"])
    if not args["debug"]:
        run_binning_method(out_dir, os.path.join(out_dir, "fasta_files"),
                           args["method"], args["cpu"], args["verbose"])

    bin(args["job_id"], out_dir, args["distance"], args["species_cutoff"],
        new_loc_to_old_loc, args["verbose"])

    remove_temp(args["job_id"], args["keep_temp"], out_dir, new_loc_to_old_loc, verbose=args["verbose"])

    return


def init():
    parser = argparse.ArgumentParser(
        description='Bin the input genomes according to sequence identity',
        usage='python bin_genomes.py [options] <job_id> <genomes_file>')
    parser.add_argument('--gff',  action='store_true',
                        help='Set if input in GFF format [Default: FASTA]', default=False)
    parser.add_argument('--species_cutoff', type=float,
                        help='Maximum distance between species to remove contaminents [%(default)s]', default=0.04, metavar='FLOAT')
    parser.add_argument('--distance', type=float,
                        help='Maximum distance between two genomes to be considered in same bin [%(default)s]', default=0.005, metavar='FLOAT')
    parser.add_argument('--max_contigs', type=int, metavar='INT',
                        help='Skip genomes with more than num_contigs [%(default)s]', default=600)
    parser.add_argument('--min_length', type=float, metavar='FLOAT',
                        help='Skip genomes shorter than this length, in MBP [%(default)s]', default=4)
    parser.add_argument('--max_length', type=float, metavar='FLOAT',
                        help='Skip genomes longer than this length, in MBP [%(default)s]', default=6)
    parser.add_argument('--cpu', type=int,
                        help='Number of CPUs to use [%(default)s]', default=16, metavar='INT')
    parser.add_argument('--keep_temp', action='store_true',
                        help='Keep temporary files', default=False)
    parser.add_argument('--verbose', action='store_true',
                        help='Verbose output while run', default=False)
    parser.add_argument('--debug', action='store_true',
                        help='Set for Debug mode (doesnt run MASH)', default=False)
    parser.add_argument(
        '--method', type=str, help='Method to bin the genomes. Options: (MASH), [%(default)s]', default="MASH", metavar='STR')
    parser.add_argument(
        'job_id', help="Name used to describe this run", metavar='STR')
    parser.add_argument(
        'genomes_file', help="File with list of genome files (Default: FASTA, set --gff for GFF)", metavar='FILE')

    args = vars(parser.parse_args())
    run(args)
    return


if __name__ == "__main__":
    init()
