
line_num = 0

techs = set()


out = open("enterobase_fixed.txt", "w")
with open("enterobase.txt") as f:
	for line in f:
		toks = line.strip().split("\t")

		## removed a single sequence that doesn't parse correctly
		if len(toks)!=44:
			print("PROBLEM!!!")
			print(line_num)
			print(len(toks))
			print(line)
			continue




		line_num += 1

		for i in range(0,len(toks)):
			toks[i] = toks[i].lower()
			if toks[i] == "" or toks[i] == "unknown" or toks[i] == "na":
				toks[i] = "nd"
		if toks[6] == "1000" or toks[6] == "0":
			toks[6] = "nd"
		

		## remove sequencing using
		tech = toks[2].split(";")
		
		if tech[1] == "ls454" or tech[1] == "abi_solid" or tech[1] == "ion_torrent":
			continue

		techs.add(tech[1])

		out.write("\t".join(toks)  + "\n")

out.close()
print(techs)