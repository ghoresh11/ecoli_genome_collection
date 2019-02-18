

setwd("/Users/gh11/e_colis/genomes/READS/")

all_reads = read.table("all_run_ids.txt", header = T,
                       stringsAsFactors = F, sep = "\t")


### create conversion table from reads to identfiers that I've been using
conversion = data.frame(identifier = character(0), reads = character(0), source = character(0))

## chen 2013
chen2013 = read.table("../ExPEC/Chen2013.csv", stringsAsFactors = F, header = T, sep = ",",
                      comment.char = "")
runs_ids = all_reads$run[all_reads$identifier %in% chen2013$SRA.accession]
accessions = tolower(chen2013$WGS.accession[match(all_reads$identifier[all_reads$identifier %in% chen2013$SRA.accession],
                             chen2013$SRA.accession)])
conversion = rbind(conversion,data.frame(identifier = accessions, reads = runs_ids, source = rep("chen2013", length(runs_ids))))


## enterobase
enterobase = read.table("../NCBI/Enterobase/enterobase_minimal.csv", sep = "\t",
                        header = T) ## the enterobase file with the minimal number of genomes
runs_ids = all_reads$run[all_reads$identifier %in% toupper(enterobase$secondary.sample.id)]
accessions = tolower(enterobase$assembly.barcode[match(all_reads$identifier[all_reads$identifier %in%  toupper(enterobase$secondary.sample.id)],
                                                 toupper(enterobase$secondary.sample.id))])
conversion = rbind(conversion,data.frame(identifier = accessions, reads = runs_ids, source = rep("enterobase", length(runs_ids))))


## salipante
salipante = read.table("../ExPEC/salipante_sra_result.csv", sep = "\t", header = T)
conversion = rbind(conversion, 
                   data.frame(identifier = salipante$Library_Name,
                              reads = salipante$Run,
                              source = rep("salipante", dim(salipante)[1])))

## get ALL of pathogen detection up to find the relevant names
pathogen_detection =read.csv("../NCBI/pathogen_detection/pathogens.metadata.tsv", comment.char = "",
                             header = T, sep = "\t", stringsAsFactors = F)

merge_names <- function(conversion, reads_file, conv_file, name) {
  reads = read.table(reads_file, header = F, stringsAsFactors = F)
  
  ## conversion list with anything I could retrieve from NCBI using the "summary" from batch
  conv = read.table(conv_file, sep = ",",
                          header = T, comment.char = "", stringsAsFactors = F)
  ## other accessions  from pathogen detection
  master = pathogen_detection$wgs_master_acc[match(reads$V1,pathogen_detection$bioproject_acc)] ## either master accession
  strain = pathogen_detection$strain[match(reads$V1,pathogen_detection$bioproject_acc)] ## or strain name
  for (i in 1:length(master)){
    if (!is.na(master[i]) && master[i] == "NULL") {
      master[i] = strain[i]
    } else {
      master[i] = (strsplit(x = master[i], split = ".", fixed = T)[[1]][1])
    }
  }
  for (i in 1:length(reads$V1)) {
    if (master[i] %in% conv$reads) {
      index = which(conv$reads == master[i])
      conv[index,1] = tolower(master[i])
      conv[index,2] = reads$V1[i]
      conv[index,3] = name
    }
  } ## anything that remains with a "?" I couldn't connect back to the projects and will need artificial reads or to be
  ## looked up manually
  conversion = rbind(conversion, conv)
  return(conversion)
}


conversion = merge_names(conversion, "hazen2013_reads.txt","../Hazen2013_EHEC_EPEC/conversion.csv", "hazen2013")
conversion = merge_names(conversion, "ingle2016_reads.txt", "../EPEC/Ingle2016/conversion.csv", "ingle2016")
conversion = merge_names(conversion, "others_reads.txt","../Other_Ref_environment/conversion.csv", "others")

## no metadata file for others000, therefore just take what's in pathoge detection
reads = read.table("others000_reads.txt", header = F, stringsAsFactors = F)
## other accessions  from pathogen detection
master = pathogen_detection$wgs_master_acc[match(reads$V1,pathogen_detection$bioproject_acc)] ## either master accession
strain = pathogen_detection$strain[match(reads$V1,pathogen_detection$bioproject_acc)] ## or strain name
for (i in 1:length(master)){
  if (!is.na(master[i]) && master[i] == "NULL") {
    master[i] = strain[i]
  } else {
    master[i] = tolower((strsplit(x = master[i], split = ".", fixed = T)[[1]][1]))
  }
}
conv = data.frame(identifier = master,
                  reads = reads$V1,
                  source = rep("others000", length(master)))
conversion = rbind(conversion, conv)

## Pathogen detection -> no metadata, information matches the same IDs
pd_reads = read.table("pathogen_detection_reads.txt", stringsAsFactors = F, header = F, sep = ",")
conversion = rbind(conversion, 
                   data.frame(identifier = pd_reads$V1,
                              reads = pd_reads$V1,
                              source = rep("pathogen_detection", length(pd_reads$V1))))

## PHE -> no metadata, match SRS to SRR
phe_reads = read.table("PHE_reads.txt", stringsAsFactors = F, header = F, sep = ",")
runs_ids = all_reads$run[all_reads$identifier %in% phe_reads$V1]
accessions = phe_reads$V1[match(all_reads$identifier[all_reads$identifier %in% phe_reads$V1],
                                                 phe_reads$V1)]
conversion = rbind(conversion, 
                   data.frame(identifier = accessions,
                              reads = runs_ids,
                              source = rep("PHE", dim(phe_reads)[1])))

## add all the names of the Sanger identifiers
## ingle2016
ingle2016 = read.table("../EPEC/Ingle2016/metadata.csv", sep = ",",
                       header = T, stringsAsFactors = F, comment.char = "")
lane_ids = read.table("../EPEC/Ingle2016/ERS_to_Lane.txt", sep = ",",
                      header = T, stringsAsFactors = F, comment.char = "")
lane_ids = lane_ids[-which(lane_ids$Lane.accession == "not found"),]
lane_ids = lane_ids[match(ingle2016$Accession, lane_ids$Sample.accession ),]

lane_ids$Sample.name[which(is.na(lane_ids$Lane.accession))] = ingle2016$Accession[which(is.na(lane_ids$Lane.name))]
lane_ids$Sample.accession[which(is.na(lane_ids$Lane.accession))] = ingle2016$Accession[which(is.na(lane_ids$Lane.name))]
lane_ids$Lane.name[which(is.na(lane_ids$Lane.accession))] = ingle2016$Accession[which(is.na(lane_ids$Lane.name))]
lane_ids$Lane.accession[which(is.na(lane_ids$Sample.name))] = ingle2016$Accession[which(is.na(lane_ids$Lane.name))]
conversion = rbind(conversion,
                   data.frame(
                     identifier = lane_ids$Sample.accession,
                     reads = lane_ids$Lane.name,
                     source = rep("ingle2016", dim(lane_ids)[1])
                   ))


## bsac
bsac = read.table("../ExPEC/BSAC/BSAC_paper.csv", sep = ",", header = T,
                  stringsAsFactors = F, comment.char = "")
### get names for ERS_to_lane file
lane_ids = read.table("../ExPEC/BSAC/ERS_to_Lane.txt", sep = ",",
                      header = T, stringsAsFactors = F, comment.char = "")
lane_ids = lane_ids$Lane.name[match(bsac$ERS.accession, lane_ids$Sample.accession )]
conversion = rbind(conversion,
                   data.frame(
                     identifier = lane_ids,
                     reads = lane_ids,
                     source = rep("bsac", length(lane_ids))
                   ))

## murray
murray = read.table("../Murray/metadata_complete.csv", sep = ",", header = T,
                    stringsAsFactors = F, comment.char = "")
conversion = rbind(conversion,
                   data.frame(
                     identifier = murray$Murray.Collection.Identifier,
                     reads = murray$WTSI.sequence.tag,
                     source = rep("murray", dim(murray)[1])
                   ))
## Hayley
hayley = read.table("../brodrick2017/13073_2017_457_MOESM1_ESM.csv", sep = ",", header = T,
                    stringsAsFactors = F, comment.char = "")
conversion = rbind(conversion,
                   data.frame(
                     identifier = hayley$Sequence.ID,
                     reads = hayley$Sequence.ID,
                     source = rep("brodrick2017", dim(hayley)[1])
                   ))

### Astrid's ETEC study
etec = read.table("../ETEC/metadata.csv", sep = ",", header = T,
                  stringsAsFactors = F, comment.char = "")
### get names for ERS_to_lane file
names = read.table("../ETEC/accessions_names.csv", sep = ",",
                   header = T, stringsAsFactors = F, comment.char = "")
names_match = names[match(etec$Isolate, names$Isolate..UG.designation. ),]
etec = data.frame(names_match, etec)
### get names for ERS_to_lane file
lane_ids = read.table("../ETEC/ERS_to_Lane.txt", sep = ",",
                      header = T, stringsAsFactors = F, comment.char = "")
lane_ids = lane_ids[-which(lane_ids$Lane.accession == "not found"),]
lane_ids = lane_ids[match( etec$Isolate, lane_ids$Sample.name),]
etec = data.frame(lane_ids, etec)
conversion = rbind(conversion,
                   data.frame(
                     identifier = etec$Lane.name,
                     reads = etec$Lane.name,
                     source = rep("etec", dim(etec)[1])
                   ))



write.table(x = conversion, file = "conversion_all.csv", sep = ",", quote = F,
            col.names = T, row.names = F)
