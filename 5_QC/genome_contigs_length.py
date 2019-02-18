from Bio.SeqIO.FastaIO import SimpleFastaParser

def get_vals(assembly_file):
	num_contigs = 0
	length = 0
	with open(assembly_file) as handle:
		for values in SimpleFastaParser(handle):
			num_contigs += 1
			length += len(values[1].strip())
	return ",".join(map(str, [length, num_contigs]))

out = open("lengths.csv", "w")
out.write("ID, Assembly_file, length, num_contigs\n")

cnt = 0
with open("final_metadata_with_loc.csv") as f:
	for line in f:
		toks = line.strip().split("\t")
		if line.startswith("ID"):
			assembly_loc = toks.index("Assembly_Location")
			continue
		assemblies = toks[assembly_loc].split(",")
		for a in assemblies:
			out.write(toks[0] + "," + a + "," + get_vals(a) + "\n")

out.close()