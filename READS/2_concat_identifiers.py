import os

files = os.listdir("out/")

out = open("all_run_ids.txt", "w")
out.write("identifier\trun\turl\n")
# i = 0  
for f in files:
	if not f.endswith(".txt"):
		continue
	with open("out/" + f) as f_open:
		for line in f_open:
			if line.startswith("identifier"):
				continue
			out.write(line)
			# i += 1
			# if i == 500:
			# 	break
out.close()

