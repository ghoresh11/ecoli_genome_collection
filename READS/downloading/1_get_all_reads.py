import subprocess
import sys

CREATE_FILES = False ## set to True if creating the downloading jobs
RUN_JOB = False ## set to True if sumbitting a batch

if CREATE_FILES:
	with open("DOWNLOAD.txt") as f:  ## DOWNLOAD.txt has all the SRS ids that you need to download
		i = 0
		job_num = 1
		ids = []
		for line in f:
			line = line.strip()
			if i == 11: ## download 12 files per job, I had one job failed after downloading only 15
				ids.append(line)
				with open("download_jobs/download_" + str(job_num) + ".txt","w" ) as out: ## all the jobs will be saved in a directory "download_jobs"
					for j in ids:
						out.write(j + "\n")
				job_num += 1
				ids = []
				i = 0
			else:
				ids.append(line)
				i += 1

	## last one out of the loop
	if len(ids) > 0:
		with open("download_jobs/download_" + str(job_num) + ".txt","w" ) as out:
			for j in ids:
				out.write(j + "\n")

## now there are 588 files -> I want to do this 7 times, getting 1000 genomes at a time
## 588 divides by 84 7 times
## I'll run the job seven times in duplicates of 84

num_jobs_per_batch = 84 ## I had 7 batches of 84 files each, change depending on your number of files

if RUN_JOB:
	batch = int(sys.argv[1]) ## change here depending what batch I'm on (start with 1)
	
	for i in range((num_jobs_per_batch * (batch - 1)) + 1,(num_jobs_per_batch * (batch - 1)) + num_jobs_per_batch + 1):  ## 
		
		filename = "download_jobs/download_" + str(i) + ".txt" ## open the file with the SRR ids
		jobname = str(i)
		
		subprocess.call(["bsub", "-J", jobname,\
			"-o", jobname + ".o", "-e", jobname + ".e", \
			"bash", "X_download_fastqs.sh", filename])
