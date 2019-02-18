from country_info import countries
import re


out = open("metadata_fixed.csv", "w")

## check duplicated accessions
accessions = {}

states = [
         'Alabama','Alaska','Arizona','Arkansas','California','Colorado',
         'Connecticut','Delaware','Florida','Georgia','Hawaii','Idaho', 
         'Illinois','Indiana','Iowa','Kansas','Kentucky','Louisiana',
         'Maine', 'Maryland','Massachusetts','Michigan','Minnesota',
         'Mississippi', 'Missouri','Montana','Nebraska','Nevada',
         'New Hampshire','New Jersey','New Mexico','New York',
         'North Carolina','North Dakota','Ohio',    
         'Oklahoma','Oregon','Pennsylvania','Rhode Island',
         'South  Carolina','South Dakota','Tennessee','Texas','Utah',
         'Vermont','Virginia','Washington','West Virginia',"south carolina",
         'Wisconsin','Wyoming', 'united states', "washington d.c.", "baltimore","missourri", "seattle", "usa"]
states = [s.lower() for s in states]
uk = ["england", "scotland", "wales","great britian","great britain", "united kingdom", "uk"]


with open("complete_metadata_to_fix.csv") as f:
	for line in f:
		if line.startswith("ID"):
			out.write(line)
			continue
		line = line.lower()
		toks = line.strip().split("\t")

		##### fix the pathotype
		if toks[2] == "na" or toks[2] == "other" or toks[2] == "none" or "esbl" in toks[2] or toks[2] == "ns" or "ndm" in toks[2] or "fecal" in toks[2]:
			toks[2] = "nd"
		elif "epec" in toks[2]:
			toks[2] =  "epec"
		elif "lab" in toks[2]:
			toks[2] = "lab"
		elif "aeec" in toks[2]:
			toks[2] = "aeec"
		elif "shigella" in toks[2] or 'e. albertii' in toks[2]:
			toks[2] = "not e_coli"
		elif 'no157' in toks[2]:
			toks[2] = "non-o157 ehec"
		elif "expec" in toks[2]:
			toks[2] = "expec"
		toks[2] = toks[2].replace(" ","")


		### fix country
		if toks[3] == "na" or toks[3] == "other" or toks[3] == "none" or toks[3] == "car" or toks[3] == "ns" or toks[3] == "nk" or toks[3] == "unresolved":
			toks[3] = "nd"
		elif toks[3] in states or "usa" in toks[3]:
			toks[3] = "united states"
		elif toks[3] in uk or "england" in toks[3]:
			toks[3] = "united kingdom"
		elif toks[3] == "bali":
			toks[3] = "indonesia"
		toks[3] = toks[3].strip()
		


		### fix continent
		if toks[3] == "nd":
			toks[4] = "nd"
		else:
			flag = False
			for item in countries:
				if item["name"].lower() == toks[3]:
					toks[4] = item["continent"]
					flag = True
					break
			if not flag:
				toks[4] = "nd"

		## fix year
		if toks[5] in ["na", "ns", "no data"]:
			toks[5] = "nd"
		elif "2000s" in toks[5]:
			toks[5] = "2000s"

		## toks[6] is the pateint ID

		### Source
		if toks[7] in ["2001","2003","2002", "1987", "na","ns"]:
			toks[7] = "nd"
		elif toks[7] in ['pig', 'lion', 'ape', 'sheep', 'celebese ape', 'cougar', 'elephant', 'goat', 'leopard', 'giraffe', 'kangaroo rat',  'dog', 'marmoset', 'gorilla', 'bison', 'orangutan', 'steer']:
			toks[7] = "animal"
		elif "human" in toks[7]:
			toks[7] = "human"
		
		## Isolation
		if toks[8] in ["ns", "na", "", "not collected"]:
			toks[8] = "nd"
		elif "urine" in toks[8]:
			toks[8] = "urine"
		elif "blood" in toks[8]:
			toks[8] = "blood"
		elif "stool" in toks[8] or "feces" in toks[8] or "faeces" in toks[8]:
			toks[8] = "feces"
		elif "environ" in toks[8]:
			toks[8] = "environment"
		elif "homo" in toks[8]:
			toks[8] = "human"
		else:
			toks[8] = toks[8].split(";")[0]
		

		### phylogroup
		if toks[9] not in [ "e" , "f" , "nd" ,"a" , "b2" ,"d" , "b1" ]:
			toks[9] = "nd"
		

		### Age
		if "age" in toks[10]:
			toks[10] = toks[10].split("age")[-1].split(";")[0].split(":")[-1].replace(" ","")

		if re.search('[a-zA-Z]', toks[10]) or toks[10] == "":
			toks[10] = "nd"
		

		## sex
		if "female" in toks[11]:
			toks[11] = "f"
		elif "male" in toks[11]:
			toks[11] = "m"
		elif toks[11] != "m" and toks[11] != "f":
			toks[11] = "nd"
		
		## keep track of duplicates
		if toks[0] in accessions:
			accessions[toks[0]].append(toks)
		else:
			accessions[toks[0]] = [toks]


for a in accessions:
	if len(accessions[a]) > 1:
		final = []
		for i in range(0, len(accessions[a][0])):
			vals = set()			
			for j in accessions[a]:
				if j[i] == "nd":
					continue
				vals.add(j[i])
			if len(vals) == 0:
				final.append("nd")
			else: 
				final.append("/".join(list(vals)))
		out.write("\t".join(final) + "\n")
	else:
		out.write("\t".join(accessions[a][0]) + "\n") ## no duplicates



