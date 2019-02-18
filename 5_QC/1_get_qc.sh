


### to run
job_name=tutorial1
bsub -J ${job_name} -G team216 -o ${job_name}.o -e ${job_name}.e -R"select[mem>500] rusage[mem=500]" -M500 -n8 -R"span[hosts=1]" bash get_qc.sh ${job_name}

for f in $(<${1}/${1}_samples.txt)
do
	#echo $f
	#pf assembly --type  lane --id $f --symlink ${1}/ 
	#pf annotation --type  lane --id $f --filetype gff --symlink ${1}/ 
	#pf data --type lane --id $f --filetype fastq --symlink ${1}/reads/
	pf qc -t sample -i $f 
done 


