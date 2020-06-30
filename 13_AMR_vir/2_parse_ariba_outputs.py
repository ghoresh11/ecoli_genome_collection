import os

''' because I'm unsure about the ARIBA outputs, I will go over all the report files to 
1. split them into the different databases
2. Make sure the length coverage is correct and see when a gene is truncated
3. Look how many reads there are - maybe it will indicate copies?
'''
class Report:

	def __init__(self, name):
		self.name = name
		self.dbs = {
		"ARGannot" : set(),
		"card" : set(),
		"resfinder" : set(),
		"stx" : set(),
		"VFDB" : set(),
		"virulence" : set(),
		"plasmid" : set()
		}
		return


	def print_database(self, database_name, order):
		''' print a vector of 0 and 1 in the order the the required db'''
		out = self.name
		for item in order:
			if item in self.dbs[database_name]:
				out += ",1"
			elif item + "*" in self.dbs[database_name]:
				out += ",1*"
			else:
				out += ",0"
		return out

ordered_outputs = {"ARGannot" : set(),
		"card" : set(),
		"resfinder" : set(),
		"stx" : set(),
		"VFDB" : set(),
		"virulence" : set(),
		"plasmid" : set()}

directories = os.listdir(".")
cnt = 0
all_reports = []
for d in directories:
	if not os.path.isfile(os.path.join(d, "report.tsv")):
		continue
	cnt += 1
	print("Parsing file: %d..." %cnt)
	curr = Report(d)
	all_reports.append(curr)
	with open(os.path.join(d, "report.tsv")) as f:
		for line in f:
			toks = line.strip().split("\t")
			if line.startswith("#"):
				read_index = toks.index("reads")
				ref_len = toks.index("ref_len")
				aligned_len = toks.index("ref_base_assembled")
				continue
			db = toks[0].split("_")[0]
			hit = "_".join(toks[0].split("_")[1:])
			alignment_length = float(toks[aligned_len])
			reference_length = float(toks[ref_len])

			if alignment_length/reference_length < 0.8:
				curr.dbs[db].add(hit + "*") ## mark as a truncated version
			else:
				curr.dbs[db].add(hit)
				ordered_outputs[db].add(hit)

for key in ordered_outputs:
	out = open(key + ".csv", "w")
	order = list(ordered_outputs[key])
	out.write("Genome," + ",".join(order) + "\n")
	for r in all_reports:
		out.write(r.print_database(key, order) + "\n")
	out.close()

				
