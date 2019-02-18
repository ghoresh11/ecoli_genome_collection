
out = open("conversion.csv", "w")
out.write("identifier, reads, source\n")

projects = {}

with open("project_to_name_ingle2016.txt") as f:
	for line in f:
		if line.startswith("Organism"):

			if "str." in line:
				name = line.strip().split("str. ")[-1]
				name = name.split()[0]
			else:
				name = line.strip().split("coli ")[-1]
				name = name.split("(Taxonomy")[0]
				if len(name) == 0:
					continue
				
				name = name.split()[-1]

				if name == "AG'":
					name = "BL21"

			
			continue

		if line.startswith("BioProject"):
			project = line.strip().split()[-1]
			projects[name.replace(" ","")] = project.replace(" ","")
			continue
		


with open("metadata2.csv") as f:
	for line in f:
		if line.startswith("Strain"):
			continue
		toks = line.strip().split(",")
		if toks[-1].startswith("ERS"): ## sequenced
			continue
		name = toks[0].replace(" ","")
		
		if name not in projects:
			out.write(",".join([name, toks[-1], "ingle2016?"]) + "\n")
		else:
			out.write(",".join([toks[-1].lower(), projects[name], "ingle2016"]) + "\n")


out.close()