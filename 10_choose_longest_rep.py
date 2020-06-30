import os
from Bio.SeqIO.FastaIO import SimpleFastaParser
from random import sample


''' choose the longest rep for the pan-genome reference file
I already found mistakes in the combined roary output because of
a short rep being chosen, I'm going to rerun eveyrthing with
this correction and see how it changes things.'''

if __name__ == "__main__":

    for cluster in range(1,52):
        cluster = str(cluster)
        roary_dir = os.path.join("/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/new_roary", cluster)
        cds_dir = os.path.join("/nfs/pathogen004/gh11/", cluster)

        if not os.path.isdir(roary_dir):
            continue

        cds_to_gene = {} ## dictionary pointing between a CDS and the gene cluster
        with open(os.path.join(roary_dir, "clustered_proteins")) as f:
            for line in f:
                toks = line.strip().split()
                gene = toks[0].replace(":","")
                for m in toks[1:]:
                    cds_to_gene[m] = gene

        gene_to_rep = {}
        all_cdss_files = os.listdir(cds_dir)
        for fasta_file in all_cdss_files:
            with open(os.path.join(cds_dir,fasta_file)) as handle:
                for values in SimpleFastaParser(handle):
                    if values[0] not in cds_to_gene:
                        continue
                    gene = cds_to_gene[values[0]]
                    if gene not in gene_to_rep:
                        gene_to_rep[gene] = []
                    gene_to_rep[gene].append(len(values[1])) ## keep track of lengths
        for gene in gene_to_rep:
            try:
                gene_to_rep[gene] = mode(gene_to_rep[gene])
            except:
                ## this is a problem actually.. this means that roary clustered these badly:
                gene_to_rep[gene].sort()
                mid_value = int(len(gene_to_rep[gene])/2)-1
                gene_to_rep[gene] = gene_to_rep[gene][mid_value]

        for fasta_file in all_cdss_files:
            with open(os.path.join(cds_dir,fasta_file)) as handle:
                for values in SimpleFastaParser(handle):
                    if values[0] not in cds_to_gene:
                        continue
                    gene = cds_to_gene[values[0]]
                    if gene_to_rep[gene] == len(values[1]):
                        gene_to_rep[gene] = {"member": values[0], "seq":values[1]}

        with open(os.path.join(roary_dir, "mode_pan_genome_reference.fa"),"w") as out:
            for gene in gene_to_rep:
                try:
                    out.write(">" + gene_to_rep[gene]["member"] + " " + gene + "\n" + gene_to_rep[gene]["seq"] + "\n")
                except:
                    pass
