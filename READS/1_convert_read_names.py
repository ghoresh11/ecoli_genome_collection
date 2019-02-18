import os
import subprocess
import sys


f = sys.argv[1] 

identifiers = set()
with open(f) as f_open:
	for line in f_open:
		line = line.strip()
		if line == "":
			continue
		curr_identifiers = line.split(",")
		for i in curr_identifiers:
			identifiers.add(i)


identifiers = list(identifiers)

for i in identifiers:
	subprocess.call(["python", "getSraRunsFromAccIds.py",\
		"--identifiers", i, "--outStub", "out/" + i, "--overWrite"])


# ## to run
for filename in *_reads.txt; do
    bsub -J ${filename} -G team216 -o ${filename}.o -e ${filename}.e python convert_read_names.py ${filename}
done