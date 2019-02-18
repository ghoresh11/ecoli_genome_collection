import os

'''Locations of all the different E. coli genomes that have been downloaded from all the different sources'''


### FASTA dir ####

brodrick2017_fasta = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/genomes/brodrick2017/"
chen2013_fasta = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/genomes/chen2013/complete/"
human_fasta = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/genomes/enterobase/human/complete/"
ecor_fasta = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/genomes/enterobase/ecor/complete/"
hazen2013_fasta = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/genomes/hazen2013/complete/"
nctc_fasta = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/genomes/nctc/complete/"
salipante2014_fasta = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/genomes/salipante2014/complete/"
others0000_fasta = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/genomes/others0000/complete/"
bsac_fasta = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/genomes/bsac/"
etec_fasta = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/genomes/etec/complete/"
hazen2016_fasta = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/genomes/hazen2016/complete/"
ingle2016_fasta = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/genomes/ingle2016/complete/"
murray_fasta = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/genomes/murray/"
other_fasta = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/genomes/other/complete/"
missing_assemblies = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/genomes/missing_and_added/assemblies/"

## Downloaded sequences
downloaded = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/genomes/downloaded/"
study_names = ["ecoli_mobilome", "ecoli_mobilome_1_2", "ecoli_mobilome_3_4", "ecoli_mobilome_5_6", "ecoli_mobilome_7_long", "exist"]

downloaded_reads = []
downloaded_assemblies = []
downloaded_annotations = []
for study in study_names:
	downloaded_reads.append(os.path.join(downloaded, study, "reads"))
	downloaded_assemblies.append(os.path.join(downloaded, study, "assembly"))
	downloaded_annotations.append(os.path.join(downloaded, study, "annot"))

all_dirs_fasta =  downloaded_assemblies + [missing_assemblies, salipante2014_fasta, brodrick2017_fasta, chen2013_fasta, human_fasta , ecor_fasta, hazen2013_fasta, nctc_fasta, others0000_fasta, bsac_fasta, etec_fasta, hazen2016_fasta, ingle2016_fasta, murray_fasta, other_fasta]

###### GFF dirs  ####

brodrick2017_gff = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/genomes/brodrick2017/"
chen2013_gff = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/genomes/chen2013/annotation/"
human_gff = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/genomes/enterobase/human/annotation/"
ecor_gff = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/genomes/enterobase/ecor/annotation/"
hazen2013_gff = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/genomes/hazen2013/annotation/"
nctc_gff = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/genomes/nctc/annotation/"
others0000_gff = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/genomes/others0000/annotation/"
bsac_gff = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/genomes/bsac/"
etec_gff = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/genomes/etec/annotation/"
hazen2016_gff = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/genomes/hazen2016/annotation/"
ingle2016_gff = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/genomes/ingle2016/annotation/"
murray_gff = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/genomes/murray/"
other_gff = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/genomes/other/annotation/"
salipante2014_gff = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/genomes/salipante2014/annotation/"
missing_annots = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/genomes/missing_and_added/annots/"

all_dirs_gff =  downloaded_annotations + [missing_annots, salipante2014_gff, brodrick2017_gff, chen2013_gff, human_gff, ecor_gff, hazen2013_gff, nctc_gff, others0000_gff, bsac_gff, etec_gff, hazen2016_gff, ingle2016_gff, murray_gff, other_gff]


### READS dirs  ####

brodrick2017_reads = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/genomes/brodrick2017/reads/"
nctc_reads = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/genomes/nctc/reads/"
bsac_reads = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/genomes/bsac/reads/"
etec_reads = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/genomes/etec/reads/"
ingle2016_reads = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/genomes/ingle2016/reads/"
murray_reads = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/genomes/murray/reads/"
artificial_reads = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/genomes/reads/artificial_reads/"

all_dirs_reads =  downloaded_reads + [artificial_reads, brodrick2017_reads, nctc_reads, bsac_reads, etec_reads, ingle2016_reads, murray_reads]
