#!/usr/bin/env python3

import argparse
import os
import subprocess
import sys

class Error (Exception): pass


def syscall(command):
    completed_process = subprocess.run(command, shell=True, stderr=subprocess.STDOUT, stdout=subprocess.PIPE, universal_newlines=True)
    if completed_process.returncode != 0:
        print('Error running this command:', command, file=sys.stderr)
        print('Return code:', completed_process.returncode, file=sys.stderr)
        print('Output from stdout and stderr:', completed_process.stdout, sep='\n', file=sys.stderr)
        raise Error('Error in system call. Cannot continue')

    return completed_process


def run_trimmomatic(
        reads1,
        reads2,
        out1,
        out2,
        trimmo_root=None,
        adapters='TruSeq3-PE-2.fa', ## what Kate also used, must be the most standard illumina adaptors
        minlen=36,
        verbose=1,
        threads=5,
        qual_trim='LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15',
        adapters_included=True,
        quality_encoding='phred33'):

    trimmo_root = '/software/pathogen/external/apps/usr/bin/' ## where the trimmomatic files sit
    #jar_files = [x for x in os.listdir(trimmo_root) if x.endswith('.jar')]
    # if len(jar_files) != 1:
    #     raise Error('Error finding Trimmoatic jar file in directory "' + trimmo_root + '". Found ' + str(len(jar_files)) + ' jar files. Cannot continue')
    jar_file = os.path.join(trimmo_root, 'trimmomatic-0.33.jar') ## this is the name of the jar file found in root

    # if adapters_included:
    #     adapters = os.path.join(trimmo_root, 'adapters', adapters) ## which adaptors to ise
    adapters = "TruSeq3-PE-2.fa" ## change here for adaptors file

    if not os.path.exists(adapters):
        raise Error('Cannot find adapters file "' + adapters + '".')

    if not os.path.isfile(reads2):
        ## theres only a single read file
            ## the actual command it runs
        cmd = ' '.join([
            'java -Xmx1000m -jar',
            jar_file, 
            'SE', ## quality encoding is automatic
            '-threads', str(threads),
            reads1,
            out1,
            'ILLUMINACLIP:' + os.path.abspath("TruSeq3-SE.fa") + ':2:30:10',
            qual_trim,
            'MINLEN:' + str(minlen) 
        ])
    else:
        ## the actual command it runs
        cmd = ' '.join([
            'java -Xmx1000m -jar',
            jar_file, 
            'PE', ## quality encoding is automatic
            '-threads', str(threads),
            reads1,
            reads2,
            out1,
            '/dev/null',
            out2,
            '/dev/null',
            'ILLUMINACLIP:' + os.path.abspath(adapters) + ':2:30:10',
            qual_trim,
            'MINLEN:' + str(minlen) 
        ])

    if verbose:
        print('Run trimmomatic:', cmd)
    syscall(cmd)



reads_dir = sys.argv[1]
trimmed_dir = sys.argv[2]
read_names_file = os.path.join("long_jobs", "download_" + sys.argv[3] + ".txt") ## the file name with the reads to be trimmed

dont_trim = [] ## some SRRs I've found are not illumina sequenced, so I won't trim them
with open("dont_trim.txt") as f:
    for line in f:
        dont_trim.append(line.strip())

with open(read_names_file) as f:
    for line in f:
        name = line.strip()        
        if name in dont_trim:
            print("TAKE CARE!!! %s has been marked not to be trimmed!") ## just copy the files?
            continue 


        file1 = os.path.join(reads_dir, name + "_1.fastq.gz" )
        file2 = os.path.join(reads_dir, name +  "_2.fastq.gz")

        if not os.path.isfile(file1):
            continue 

        out1 = os.path.join(trimmed_dir,  name + "_trimmed_1.fastq.gz")  
        out2 = os.path.join(trimmed_dir,  name + "_trimmed_2.fastq.gz")

        run_trimmomatic(
            file1, file2,
            out1, out2)

## to run this job:
## arg1=downloaded/ ## here: the reads directory
## arg2=1
## job_name=trimming${arg2}
## bsub -J ${job_name} -G team216 -o ${job_name}.o -e ${job_name}.e -R"select[mem>2000] rusage[mem=2000]" -M2000 -n5 -R"span[hosts=1]" python3.6 2_trimming.py ${arg1} ${arg2}


