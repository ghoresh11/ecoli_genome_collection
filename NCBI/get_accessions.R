library(ggplot2)

setwd("/Users/gh11/e_colis/genomes/NCBI//")


## sequences from "pathogen detection"
all_pathogen_detection = read.csv("pathogen_detection/pathogens.metadata.tsv", comment.char = "",
                                    header = T, sep = "\t", stringsAsFactors = F)

## sequences from Enterobase
all_enterobase = read.table("Enterobase/enterobase_fixed.txt", sep = "\t", stringsAsFactors = F,
                            header = T, comment.char = "", quote = "")

## sequences from PHE
all_PHE = read.table("PHE/sra_result.csv", sep = ",", header = T, stringsAsFactors = F, comment.char = "")
all_PHE = all_PHE[which(all_PHE$Organism.Name == 'Escherichia coli'),]
write.table(x = all_PHE$Sample.Accession, file = "PHE/phe_reads.txt",
            quote = F, col.names = F, row.names = F)
write.table(x = all_PHE, file = "PHE/phe_minimal.csv", col.names = T, row.names = F,
            quote = F, sep = ",")

## reads from pathogen detection
all_pathogen_detection = all_pathogen_detection[which(all_pathogen_detection$epi_type == "clinical"),]
all_pathogen_detection = all_pathogen_detection[which(!grepl(x = tolower(all_pathogen_detection$scientific_name), 
                                                            pattern = "shigella")), ]

all_pathogen_detection = all_pathogen_detection[which(all_pathogen_detection$Run!="NULL"), ]

write.table(x = all_pathogen_detection$Run, file = "pathogen_detection/pathogen_detection_reads.txt",
            row.names = F, col.names = F,
            quote = F)
write.table(x = all_pathogen_detection, file = "pathogen_detection/pathogen_detection_minimal.txt",
            row.names = F, col.names = T, sep = ",",
            quote = T)

## enterobase, only genomes from sources I can use
### have a look at the different sources
sources = data.frame(table(all_enterobase$lab.contact))
sources = sources[order(sources$Freq, decreasing = T),]
sources = sources[which(sources$Freq>50),]
sources$Var1 =  factor(sources$Var1, sources$Var1)
ggplot(sources, aes(x = Var1, y = Freq)) + geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


oked_sources = c("centers for disease control and prevention enteric diseases laboratory branch",
                 "public health england",
                 "fda center for food safety and applied nutrition",
                 "the sanger center",
                 "edlb-cdc",
                 "broad institute", ## oked, except for 41 strains which I was told not to use
                 "the wellcome trust sanger institute",
                 "united states department of agriculture, food safety and inspection service",
                 "gastrointestinal bacteria reference unit, public health england",
                 "food and drug administration, center for food safety and applied nutrition",
                 "enteric diseases laboratory branch, centers for disease control and prevention",
                 "fda/cfsan", "fda", "centers for disease control and prevention", 
                 "centers for disease control and prevention - clinical and environmental branch")
all_enterobase = all_enterobase[which(all_enterobase$lab.contact %in% oked_sources),]
all_enterobase = all_enterobase[which(all_enterobase$source.niche == "human"),] ## take only human isolates

write.table(x = all_enterobase, file = "Enterobase/enterobase_minimal.csv", sep = "\t", row.names = F, col.names = T,
            quote = F)
write.table(x = toupper(all_enterobase$secondary.sample.id), file = "Enterobase/enterobase_reads.txt", 
            col.names = F, row.names = F, quote = F)

