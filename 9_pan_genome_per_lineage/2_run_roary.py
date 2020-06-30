import os
import subprocess

## Roary usage:
# Usage:   roary [options] *.gff
#
# Options: -p INT    number of threads [1]
#          -o STR    clusters output filename [clustered_proteins]
#          -f STR    output directory [.]
#          -e        create a multiFASTA alignment of core genes using PRANK
#          -n        fast core gene alignment with MAFFT, use with -e
#          -i        minimum percentage identity for blastp [95] ## leaving on 95
#          -cd FLOAT percentage of isolates a gene must be in to be core [99]
#          -b STR    blastp executable [blastp]
#          -m STR    makeblastdb executable [makeblastdb]
#          -s        dont split paralogs


jobs_dir= "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/new_roary/new_jobs/"
roary_out_dir= "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/new_roary/"
jobs = os.listdir(jobs_dir)
threads = "16"
blastp = "/software/pubseq/bin/ncbi_blast+/blastp"
makeblastdb = "/software/pubseq/bin/ncbi_blast+/makeblastdb"

for jobs_file in jobs:
    if not jobs_file.endswith(".txt"):
        continue
    cluster_id = jobs_file.split("_")[1].replace(".txt","")
    gffs = []
    with open(os.path.join(jobs_dir, jobs_file)) as f:
        for line in f:
            gffs.append(line.strip())

    ## if some failed: apply here:
    # if cluster_id not in ["10", "11", "13", "16", "17","19", "31", "33", "37", "7"]:
    #     continue
    ## calculate how much memory is required (more or less)
    ## according to my calcs each genomes needs about 30MBs
    # I might be off, but I'll start with lower memory request
    mem = 30 * len(gffs)
    mem = max(mem, 10000)
    queue = "parallel"

    if mem > 200000:
        queue = "hugemem"

    mem = str(mem)
    ## create output directory for this cluster
    outdir = os.path.join(roary_out_dir, cluster_id)
    try:
        os.makedirs(outdir)
    except Exception:
        pass

    job_name = "cluster_"+cluster_id


    lsf_prefix = ["bsub", "-q", queue, "-J", job_name, "-G", "team216","-o", job_name + ".o",
     "-e", job_name + ".e", '-R"select[mem>' + mem + '] rusage[mem='+ mem + '] span[hosts=1]"', '-M' + mem, "-n" + threads]

    command = map(str,["roary", "-p", threads, "-f", outdir, "-e", "-mafft", "-b", blastp, "-m", makeblastdb, "-s"]+ gffs)
    ## submit the job
    subprocess.call(lsf_prefix + command)
