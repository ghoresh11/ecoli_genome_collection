
import sys

## for each genome, see if the pipelines have completed

tot = 0
tot_imp = 0
tot_ass = 0
failed_ass = 0
tot_annot = 0

with open(sys.argv[1]) as f:
	for line in f:
		if line.startswith("Name"):
			continue

		toks = line.strip().split()
		annot = toks[-1]
		assem = toks[-2]
		imp = toks[1]

		tot += 1
		if annot == "Done":
			tot_annot += 1
		if imp == "Done":
			tot_imp += 1
		if assem == "Done":
			tot_ass += 1
		if "Failed" in assem:
			failed_ass += 1


print("Total sequences : %d" %tot)
print("Total imported: %d" %tot_imp)
print("Total done assembled: %d" %tot_ass)
print("Total failed assembled: %d" %failed_ass)
print("Total assembled: %d" %(tot_ass + failed_ass))
print("Total annotated: %d" %tot_annot)
