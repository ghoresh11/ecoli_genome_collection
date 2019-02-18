#!/bin/bash

## create directories for genome files
# mkdir -p ${1}
#mkdir -p ${1}/reads/

# # get all the relevant files for the strains
# for f in $(<${1}/${1}_lanes.txt)
# do
# 	echo $f
# 	pf assembly --type  lane --id $f --symlink ${1}/ 
# 	pf annotation --type  lane --id $f --filetype gff --symlink ${1}/ 
# done 

for f in $(<${1}/${1}_samples.txt)
do
	#echo $f
	#pf assembly --type  lane --id $f --symlink ${1}/ 
	#pf annotation --type  lane --id $f --filetype gff --symlink ${1}/ 
	#pf data --type lane --id $f --filetype fastq --symlink ${1}/reads/
	pf qc -t sample -i $f 
done 

exit
## to get all the files for a complete study
study=ecoli_mobilome_1_2
study=ecoli_mobilome_3_4
study=ecoli_mobilome_5_6_new
study=ecoli_mobilome_7
pf status -t study -i $study > status_${study}

pf annotation -t study -i $study --filetype gff --symlink .
pf assembly -t study -i $study --symlink .
pf data -t study -i $study --filetype fastq --symlink .


## to run
## job_name=bsac
## bsub -J ${job_name} -G team216 -o ${job_name}.o -e ${job_name}.e -R"select[mem>500] rusage[mem=500]" -M500 bash get_sanger_sequences.sh ${job_name}