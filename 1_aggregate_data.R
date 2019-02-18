setwd("/Users/gh11/e_colis/genomes/")

## go through each data file and add the relevant rows to the complete file
## if anything is missing complete it with "ND"
## In the end, add the run IDs from the conversion table

enterobase = read.table(file = "NCBI/Enterobase/enterobase_minimal.csv", 
                        sep = "\t", header = T, comment.char = "", 
                        quote = "", stringsAsFactors = F)

## remove from enterobase dataframe the sequences I was told not to use by the broad
broad_remove = read.table("READS/broad_dont_use.csv", sep = ",", header = T, stringsAsFactors = F)
remove = unique(match(tolower(broad_remove$Assembly.barcode), enterobase$assembly.barcode))
remove = remove[-which(is.na(remove))]
enterobase = enterobase[-remove, ]

add_from_enterobase <- function(niche){
  eb = enterobase[which(enterobase$source.niche == niche),]
  eb_metadata = data.frame(ID = eb$assembly.barcode,
                                 Name = eb$assembly.barcode,
                                 Pathotype = eb$simple.patho,
                                 Country = eb$country,
                                 Continent = eb$continent,
                                 Year = eb$collection.year,
                                 Patient_ID = rep("ND",dim(eb)[1]),
                                 Source = eb$source.niche,
                                 Isolation = eb$source.details,
                                 Phylgroup = eb$ecor.cluster,
                                 Patient_age = eb$source.details,
                                 Patient_sex = eb$source.details,
                                 Publication = eb$lab.contact,
                                 Notes = eb$comment,
                                 Removed = rep("no", dim(eb)[1]),
                           stringsAsFactors = F)
  return(eb_metadata)
}


### retrieve metadata that I chose to use for human and environmental isolates
human_eb = add_from_enterobase("human")

#### ECOR collection ####
ecor = enterobase[which(startsWith(enterobase$name, "ecor-")),]
ecor_data = read.table("NCBI/Enterobase//ecor.csv", sep = ",", quote = "",
                       comment.char = "", header = T, stringsAsFactors = F)
ecor_ID = ecor$assembly.barcode
ecor_name = ecor$name

## get metadata from other file
ecor_matched = ecor_data[match(ecor$name,ecor_data$Isolate), ]
ecor_metadata = data.frame(ID = ecor_ID,
                           Name = ecor_name,
                           Pathotype = ecor_matched$Pathotype,
                           Country = ecor_matched$Locale,
                           Continent = ecor_matched$Continent,
                           Year = ecor_matched$Date,
                           Patient_ID = rep("ND",dim(ecor)[1]),
                           Source = ecor_matched$Host,
                           Isolation = rep("ND", dim(ecor)[1]),
                           Phylgroup = ecor_matched$Phylogroup,
                           Patient_age = ecor_matched$Age,
                           Patient_sex = rep("ND", dim(ecor)[1]),
                           Publication = rep("Ochman84_ECOR", dim(ecor)[1]),
                           Notes = rep("", dim(ecor)[1]),
                           Removed = rep("no", dim(ecor)[1]),
                           stringsAsFactors = F)

### Hazen 2016
hazen2016 = read.table("EPEC/Hazen2016/metadata.csv", sep = ",",
                       header =  T, stringsAsFactors = F, comment.char = "")
hazen2016_metadata = data.frame(ID = hazen2016$Accession.no.,
                           Name = hazen2016$Specimen.ID,
                           Pathotype = rep("EPEC", dim(hazen2016)[1]),
                           Country = hazen2016$Location,
                           Continent = rep("Asia/Africa", dim(hazen2016)[1]),
                           Year = rep("2008-2012", dim(hazen2016)[1]),
                           Patient_ID = rep("ND",dim(hazen2016)[1]),
                           Source = rep("human", dim(hazen2016)[1]),
                           Isolation = rep("stool", dim(hazen2016)[1]),
                           Phylgroup = hazen2016$Phylogroup,
                           Patient_age = rep("<5", dim(hazen2016)[1]),
                           Patient_sex = rep("ND", dim(hazen2016)[1]),
                           Publication = rep("Hazen2016_GEMS", dim(hazen2016)[1]),
                           Notes = hazen2016$Clinical.outcome,
                           Removed = rep("no", dim(hazen2016)[1]),
                           stringsAsFactors = F)

### Ingle 2016
ingle2016 = read.table("EPEC/Ingle2016/metadata.csv", sep = ",",
                       header = T, stringsAsFactors = F, comment.char = "")
#ingle2016 = ingle2016[-which(ingle2016$Accession == ""),]
### get names for ERS_to_lane file
lane_ids = read.table("EPEC/Ingle2016/ERS_to_Lane.txt", sep = ",",
                      header = T, stringsAsFactors = F, comment.char = "")
lane_ids = lane_ids[-which(lane_ids$Lane.accession == "not found"),]
lane_ids = lane_ids[match(ingle2016$Accession, lane_ids$Sample.accession ),]

lane_ids$Sample.name[which(is.na(lane_ids$Lane.accession))] = ingle2016$Accession[which(is.na(lane_ids$Lane.name))]
lane_ids$Sample.accession[which(is.na(lane_ids$Lane.accession))] = ingle2016$Accession[which(is.na(lane_ids$Lane.name))]
lane_ids$Lane.name[which(is.na(lane_ids$Lane.accession))] = ingle2016$Accession[which(is.na(lane_ids$Lane.name))]
lane_ids$Lane.accession[which(is.na(lane_ids$Sample.name))] = ingle2016$Accession[which(is.na(lane_ids$Lane.name))]


ingle2016 = data.frame(lane_ids,ingle2016)

ingle2016_metadata = data.frame(ID = ingle2016$Lane.name,
                                Name = ingle2016$Accession,
                                Pathotype = ingle2016$Pathotype,
                                Country = ingle2016$Country,
                                Continent = rep("ND", dim(ingle2016)[1]),
                                Year = ingle2016$Year.of.Isolation,
                                Patient_ID = rep("ND",dim(ingle2016)[1]),
                                Source = rep("human", dim(ingle2016)[1]),
                                Isolation = rep("stool", dim(ingle2016)[1]),
                                Phylgroup =rep("ND",dim(ingle2016)[1]),
                                Patient_age = rep("<5", dim(ingle2016)[1]),
                                Patient_sex = rep("ND", dim(ingle2016)[1]),
                                Publication = rep("Ingle2016_GEMS", dim(ingle2016)[1]),
                                Notes = rep("-", dim(ingle2016)[1]),
                                Removed = rep("no", dim(ingle2016)[1]),
                                stringsAsFactors = F)

### BSAC
bsac = read.table("ExPEC/BSAC/BSAC_paper.csv", sep = ",", header = T,
                  stringsAsFactors = F, comment.char = "")
### get names for ERS_to_lane file
lane_ids = read.table("ExPEC/BSAC/ERS_to_Lane.txt", sep = ",",
                      header = T, stringsAsFactors = F, comment.char = "")
lane_ids = lane_ids$Lane.name[match(bsac$ERS.accession, lane_ids$Sample.accession )]
bsac_metadata = data.frame(ID = lane_ids,
                                Name = bsac$ERS.accession,
                                Pathotype = rep("ExPEC", dim(bsac)[1]),
                                Country = rep("UK", dim(bsac)[1]),
                                Continent = rep("europe", dim(bsac)[1]),
                                Year = bsac$Year_of_isolation,
                                Patient_ID = rep("ND",dim(bsac)[1]),
                                Source = rep("human", dim(bsac)[1]),
                                Isolation = rep("blood", dim(bsac)[1]),
                                Phylgroup = bsac$Phylogroup,
                                Patient_age = rep("ND", dim(bsac)[1]),
                                Patient_sex = rep("ND", dim(bsac)[1]),
                                Publication = rep("kallonen2017_BSAC", dim(bsac)[1]),
                                Notes = rep("-", dim(bsac)[1]),
                                Removed = rep("no", dim(bsac)[1]))


### Chen 2013
chen2013 = read.table("ExPEC/Chen2013.csv", sep = ",", header = T,
                             stringsAsFactors = F, comment.char = "")
chen2013_metadata = data.frame(ID = chen2013$WGS.accession,
                           Name = chen2013$WGS.accession,
                           Pathotype = rep("ExPEC", dim(chen2013)[1]),
                           Country = rep("USA", dim(chen2013)[1]),
                           Continent = rep("north america", dim(chen2013)[1]),
                           Year = rep("2003-2006", dim(chen2013)[1]),
                           Patient_ID = chen2013$Patient,
                           Source = rep("human", dim(chen2013)[1]),
                           Isolation = chen2013$Source,
                           Phylgroup = rep("ND", dim(chen2013)[1]),
                           Patient_age = rep("18-49", dim(chen2013)[1]),
                           Patient_sex = rep("F", dim(chen2013)[1]),
                           Publication = rep("Chen2013", dim(chen2013)[1]),
                           Notes = rep("-", dim(chen2013)[1]),
                           Removed = rep("no", dim(chen2013)[1]),
                           stringsAsFactors = F)


## Salipante 2014
salipante2014 = read.table("ExPEC/salipante.csv", sep = ",", header = T,
                      stringsAsFactors = F, comment.char = "")
### get names identifiers
ids = read.table("ExPEC/salipante2014_names.txt", 
                      header = F, stringsAsFactors = F, comment.char = "#")
ids = ids$V3[match(salipante2014$Sample.ID, ids$V5)]
ids[which(is.na(ids))] = salipante2014$Sample.ID[which(is.na(ids))]
salipante2014_metadata = data.frame(ID = ids,
                               Name = salipante2014$Sample.ID,
                               Pathotype = rep("ExPEC", dim(salipante2014)[1]),
                               Country = rep("USA", dim(salipante2014)[1]),
                               Continent = rep("north america", dim(salipante2014)[1]),
                               Year = rep("ND", dim(salipante2014)[1]),
                               Patient_ID = salipante2014$Coded.Patient.Number,
                               Source = rep("human", dim(salipante2014)[1]),
                               Isolation = salipante2014$Isolation.Source,
                               Phylgroup = rep("ND", dim(salipante2014)[1]),
                               Patient_age = salipante2014$Patient.age,
                               Patient_sex = salipante2014$Sex,
                               Publication = rep("Salipante2014", dim(salipante2014)[1]),
                               Notes = rep("-", dim(salipante2014)[1]),
                               Removed = rep("no", dim(salipante2014)[1]),
                               stringsAsFactors = F)


### Others 
others = read.table("ExPEC/other.csv", sep = ",", header = T,
                           stringsAsFactors = F, comment.char = "")
others_metadata = data.frame(ID = others$Accession,
                                    Name = others$Name,
                                    Pathotype = others$pathotype,
                                    Country = others$Country,
                                    Continent = others$Continent,
                                    Year = others$Year,
                                    Patient_ID = rep("ND", dim(others)[1]),
                                    Source = rep("human", dim(others)[1]),
                                    Isolation = rep("urine", dim(others)[1]),
                                    Phylgroup = rep("ND", dim(others)[1]),
                                    Patient_age = others$Patient.age,
                                    Patient_sex = others$Sex,
                                    Publication = others$Paper,
                                    Notes = rep("-", dim(others)[1]),
                                    Removed = rep("no", dim(others)[1]),
                             stringsAsFactors = F)

### hazen 2013
hazen2013 = read.table("Hazen2013_EHEC_EPEC/metadata.csv", sep = ",", header = T,
                    stringsAsFactors = F, comment.char = "")
hazen2013_metadata = data.frame(ID = hazen2013$Accession.Number,
                             Name = hazen2013$Strain,
                             Pathotype = hazen2013$Phylogenetic.lineage,
                             Country = hazen2013$Location,
                             Continent = hazen2013$Location,
                             Year = rep("ND", dim(hazen2013)[1]),
                             Patient_ID = rep("ND", dim(hazen2013)[1]),
                             Source = rep("human", dim(hazen2013)[1]),
                             Isolation = rep("stool", dim(hazen2013)[1]),
                             Phylgroup = hazen2013$Phylogroup,
                             Patient_age = rep("ND", dim(hazen2013)[1]),
                             Patient_sex = rep("ND", dim(hazen2013)[1]),
                             Publication = rep("Hazen2013", dim(hazen2013)[1]),
                             Notes = rep("-", dim(hazen2013)[1]),
                             Removed = rep("no", dim(hazen2013)[1]),
                             stringsAsFactors = F)



#### murrays collection
murray = read.table("Murray/metadata_complete.csv", sep = ",", header = T,
                   stringsAsFactors = F, comment.char = "")
murray_metadata = data.frame(ID = murray$WTSI.sequence.tag,
                                Name = murray$Murray.Collection.Identifier,
                                Pathotype = rep("ND", dim(murray)[1]),
                                Country = rep("ND", dim(murray)[1]),
                                Continent = rep("ND", dim(murray)[1]),
                                Year = murray$Year,
                                Patient_ID =rep("ND", dim(murray)[1]),
                                Source = rep("human", dim(murray)[1]),
                                Isolation = murray$Source,
                                Phylgroup = rep("ND", dim(murray)[1]),
                                Patient_age = rep("ND", dim(murray)[1]),
                                Patient_sex = rep("ND", dim(murray)[1]),
                                Publication = rep("Baker2015_Murray", dim(murray)[1]),
                                Notes = rep("-", dim(murray)[1]),
                                Removed = rep("no", dim(murray)[1]))
### NCTC
nctc = read.table("NCTC/metadata.csv", sep = ",", header = T,
                    stringsAsFactors = F, comment.char = "")
### get names for ERS_to_lane file
lane_ids = read.table("NCTC/ERS_to_Lane.txt", sep = ",",
                      header = T, stringsAsFactors = F, comment.char = "")
lane_ids = lane_ids[match(lane_ids$Lane.accession , lane_ids$Lane.accession),]
nctc = lane_ids
nctc_metadata = data.frame(ID = nctc$Lane.name,
                             Name = nctc$Sample.accession,
                             Pathotype = rep("ND", dim(nctc)[1]),
                             Country = rep("UK", dim(nctc)[1]),
                             Continent = rep("Europe", dim(nctc)[1]),
                             Year = rep("ND", dim(nctc)[1] ),
                             Patient_ID =rep("ND", dim(nctc)[1]),
                             Source = rep("human", dim(nctc)[1]),
                             Isolation = rep("ND", dim(nctc)[1]),
                             Phylgroup = rep("ND", dim(nctc)[1]),
                             Patient_age = rep("ND", dim(nctc)[1]),
                             Patient_sex = rep("ND", dim(nctc)[1]),
                             Publication = rep("NCTC", dim(nctc)[1]),
                             Notes = rep("pacBio", dim(nctc)[1]),
                             Removed = rep("no", dim(nctc)[1]))


### Other - environment
env = read.table("Other_Ref_environment/metadata.csv", sep = ",", header = T,
                 stringsAsFactors = F, comment.char = "")
env_metadata = data.frame(ID = env$ID,
                          Name = env$Name,
                          Pathotype = env$Pathovar,
                          Country = env$Country,
                          Continent = env$Continent,
                          Year = env$Year,
                          Patient_ID =rep("ND", dim(env)[1]),
                          Source = env$Year,
                          Isolation = env$Isolation,
                          Phylgroup = rep("ND", dim(env)[1]),
                          Patient_age = rep("ND", dim(env)[1]),
                          Patient_sex = rep("ND", dim(env)[1]),
                          Publication = env$Publication,
                          Notes = rep("-", dim(env)[1]),
                          Removed = rep("no", dim(env)[1]),
                          stringsAsFactors = F)


### Hayley's addenbrookes collection
hayley = read.table("brodrick2017/13073_2017_457_MOESM1_ESM.csv", sep = ",", header = T,
                    stringsAsFactors = F, comment.char = "")
hayley_metadata = data.frame(ID = hayley$Sequence.ID,
                             Name = hayley$ENA.Accessions.2,
                             Pathotype = rep("ND", dim(hayley)[1]),
                             Country = rep("united kingdom", dim(hayley)[1]),
                             Continent =rep("europe", dim(hayley)[1]),
                             Year = hayley$Year.of.Isolation,
                             Patient_ID = hayley$Participant,
                             Source = rep("human", dim(hayley)[1]),
                             Isolation = hayley$Sample.Type,
                             Phylgroup = rep("ND", dim(hayley)[1]),
                             Patient_age = rep("ND", dim(hayley)[1]),
                             Patient_sex = rep("ND", dim(hayley)[1]),
                             Publication = rep("brodrick2017", dim(hayley)[1]),
                             Notes = rep("-", dim(hayley)[1]),
                             Removed = rep("no", dim(hayley)[1]),
                             stringsAsFactors = F)


### Astrid's ETEC study
etec = read.table("ETEC/metadata.csv", sep = ",", header = T,
                    stringsAsFactors = F, comment.char = "")
etec$Diarrhea.Asymptomatic.Environmental[which(etec$Diarrhea.Asymptomatic.Environmental == "D" |
                                           etec$Diarrhea.Asymptomatic.Environmental == "AS" |
                                             etec$Diarrhea.Asymptomatic.Environmental == "D/AS" ) ] = "stool"
etec$Diarrhea.Asymptomatic.Environmental[which(etec$Diarrhea.Asymptomatic.Environmental == "ENV") ] = "environment"

etec$Subjesct[which(etec$Subjet != "ENV" & etec$Subjet != "unknown") ] = "human"
etec$Subjesct[which(etec$Subjet == "ENV") ] = "environment"

### get names for ERS_to_lane file
names = read.table("ETEC/accessions_names.csv", sep = ",",
                      header = T, stringsAsFactors = F, comment.char = "")
names_match = names[match(etec$Isolate, names$Isolate..UG.designation. ),]
etec = data.frame(names_match, etec)
### get names for ERS_to_lane file
lane_ids = read.table("ETEC/ERS_to_Lane.txt", sep = ",",
                      header = T, stringsAsFactors = F, comment.char = "")
lane_ids = lane_ids[-which(lane_ids$Lane.accession == "not found"),]
lane_ids = lane_ids[match( etec$Isolate, lane_ids$Sample.name),]
etec = data.frame(lane_ids, etec)

etec_metadata = data.frame(ID = etec$Lane.name,
                             Name = etec$Isolate,
                             Pathotype = rep("ETEC", dim(etec)[1]),
                             Country = etec$Country,
                             Continent = etec$Continent,
                             Year = etec$Continent,
                             Patient_ID = rep("ND", dim(etec)[1]),
                             Source = etec$Subjesct,
                             Isolation = etec$Diarrhea.Asymptomatic.Environmental,
                             Phylgroup = etec$Phylogroup,
                             Patient_age = etec$Age,
                             Patient_sex = rep("ND", dim(etec)[1]),
                             Publication = rep("vonMentzer2014", dim(etec)[1]),
                             Notes = rep("-", dim(etec)[1]),
                             Removed = rep("no", dim(etec)[1]),
                             stringsAsFactors = F)

## PHE
phe = read.table("NCBI/PHE/phe_minimal.csv", sep = ",", stringsAsFactors = F,
                 header = T, comment.char = "", quote = "")
phe_metadata = data.frame(ID = phe$Sample.Accession,
                           Name = phe$Sample.Accession,
                           Pathotype = rep("ND", dim(phe)[1]),
                           Country = rep("UK", dim(phe)[1]),
                           Continent = rep("Europe", dim(phe)[1]),
                           Year = rep("ND", dim(phe)[1]),
                           Patient_ID = rep("ND", dim(phe)[1]),
                           Source = rep("Human", dim(phe)[1]),
                           Isolation = rep("ND", dim(phe)[1]),
                           Phylgroup = rep("ND", dim(phe)[1]),
                           Patient_age = rep("ND", dim(phe)[1]),
                           Patient_sex = rep("ND", dim(phe)[1]),
                           Publication = rep("PHE", dim(phe)[1]),
                           Notes = rep("-", dim(phe)[1]),
                           Removed = rep("no", dim(phe)[1]),
                           stringsAsFactors = F)

## Pathogen Detection
pd = read.table("NCBI/pathogen_detection/pathogen_detection_minimal.txt", sep = ",",
                stringsAsFactors = F, header = T, comment.char = "")
names = pd$strain
for (i in length(names)){
  if (names[i] == "NULL") {
    names[i] = pd$wgs_master_acc[i]
  }
}
pd_metadata = data.frame(ID = pd$Run,
                           Name = names,
                           Pathotype = rep("ND", dim(pd)[1]),
                           Country = pd$geo_loc_name,
                           Continent = rep("ND", dim(pd)[1]),
                           Year = pd$collection_date,
                           Patient_ID = rep("ND", dim(pd)[1]),
                           Source = rep("human", dim(pd)[1]),
                           Isolation = pd$isolation_source,
                           Phylgroup =  rep("ND", dim(pd)[1]),
                           Patient_age =  rep("ND", dim(pd)[1]),
                           Patient_sex = rep("ND", dim(pd)[1]),
                           Publication = rep("pathogen_detection", dim(pd)[1]),
                           Notes = rep("-", dim(pd)[1]),
                           Removed = rep("no", dim(pd)[1]),
                           stringsAsFactors = F)

complete_metadata = rbind(ecor_metadata,
                          human_eb,
                          hazen2016_metadata,
                          ingle2016_metadata,
                          bsac_metadata,
                          chen2013_metadata,
                          salipante2014_metadata,
                          others_metadata,
                          hazen2013_metadata,
                          murray_metadata,
                          nctc_metadata,
                          env_metadata,
                          hayley_metadata, 
                          etec_metadata,
                          phe_metadata)

write.table(file = "complete_metadata_to_fix.csv", x = complete_metadata,
            quote = F, row.names = F, col.names = T, sep = "\t")


## after running this, run the python script to fix all the data so it's more uniform and can be summarised




