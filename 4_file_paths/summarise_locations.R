library(ggplot2)
library(RColorBrewer)

setwd("/Users/gh11/e_colis/genomes/4_file_paths/")

md_with_locs_prev = read.table("metadata_fixed_with_loc_prev.csv", 
                          sep = "\t", stringsAsFactors = F, comment.char = "", header =T)

md_with_locs = read.table("metadata_fixed_with_loc.csv", 
                               sep = "\t", stringsAsFactors = F, comment.char = "", header =T)


## print out what the options are for the location columns
assembly_fails = sort(table(md_with_locs$Assembly_Location), decreasing = T)[1:2]

annot_fails = sort(table(md_with_locs$Annotation_Location), decreasing = T)[1]
read_fails = sort(table(md_with_locs$Reads_Location), decreasing = T)[1]


### summary of the assemblies
summary = data.frame( ID = md_with_locs$ID, 
                      name = md_with_locs$Name,
                      make_artificial = md_with_locs$Make_artificial,
                      read_ids = md_with_locs$Run_ID,
                      source = md_with_locs$Publication,
                      loc = md_with_locs$Assembly_Location,
                      cat = rep("NA", dim(md_with_locs)[1]), stringsAsFactors = F)

summary$cat[which(summary$make_artificial == "No" &
                    summary$loc != "No assembly")] = "Reads and assembly available"

summary$cat[which(summary$make_artificial == "No" &
                    summary$loc == "No assembly")] = "Reads available, assembly required" ##these should be gone

summary$cat[which(summary$make_artificial == "Yes" &
                    summary$loc != "Not found")] = "Assembly available, artifical reads required"

summary$cat[which(summary$make_artificial == "Yes" &
                    summary$loc == "Not found")] = "Problem: neither reads or assembly available"  ## these were removed entirely

for_plot = data.frame(table(summary$cat))
for_plot$Var1 = factor(for_plot$Var1,
                       c( "Reads and assembly available", "Reads available, assembly required",
                          "Assembly available, artifical reads required","Problem: neither reads or assembly available"))
ggplot(for_plot, aes(x = Var1, y = Freq)) + geom_bar(stat = "identity") +
  xlab("") + ylab("Number of genomes")+ theme_bw(base_size = 16) + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) + theme(axis.text.x = element_text(angle = 40, hjust = 1)) +
  scale_y_continuous(expand = c(0,0)) + ggtitle(paste("Total number of genomes:", sum(for_plot$Freq)))

### here: make sure I have the location of everything!! There should only be a few that failed and cannot be located!!
## how many are still missing assemblies?
missing_loc = md_with_locs[which(md_with_locs$Assembly_Location == "No assembly" |
                                   md_with_locs$Annotation_Location == "No annotation" |
                                   md_with_locs$Reads_Location == "No reads"),c(1, 2, 16, 18, 19, 20)]
write.table(x = no_assembly, file = "no_assembly.csv", sep = ",",
            col.names = T, row.names = F, quote = F)


missing_loc = which(md_with_locs$Assembly_Location == "No assembly" |
                      md_with_locs$Annotation_Location == "No annotation" |
                      md_with_locs$Reads_Location == "No reads")
md_with_locs = md_with_locs[-missing_loc,]

### final list of genomes with all the locations of the relevant files
write.table(x = md_with_locs, file = "../final_metadata_with_loc.csv", sep = "\t",
            col.names = T, row.names = F, quote = F)


#### ALl of the below add up EXACTLY to the dimentions of the missing files DF ####
### 1. Have assembly + annotation but no reads
no_reads = which(missing_loc$Assembly_Location != "No assembly" & 
                   missing_loc$Annotation_Location != "No annotation")
length(no_reads)
write.table(unlist(strsplit(missing_loc$Assembly_Location[no_reads], split = ",")), 
            file = "missing/artificial_reads.txt", col.names = F,
            row.names = F, quote = F)

## have annotation but no assembly
no_assembly_w_annot = which(missing_loc$Assembly_Location == "No assembly" &
                              missing_loc$Annotation_Location != "No annotation")
length(no_assembly_w_annot)
write.table(unlist(strsplit(missing_loc$Annotation_Location[no_assembly_w_annot], split = ",")), 
            file = "missing/retrieve_assembly.txt", col.names = F,
            row.names = F, quote = F)

### have annotation but no assembly -> require to run prokka
no_annot_w_assembly = which(missing_loc$Assembly_Location != "No assembly" &
                              missing_loc$Assembly_Location != "Not found" &
                              missing_loc$Annotation_Location == "No annotation")
length(no_annot_w_assembly)
write.table(unlist(strsplit(missing_loc$Assembly_Location[no_annot_w_assembly], split = ",")), 
            file = "missing/run_prokka.txt", col.names = F,
            row.names = F, quote = F)

## have reads but no assembly or annotation -> need to run assembly and annotation
no_assembly_or_annot_but_reads = which(missing_loc$Assembly_Location == "No assembly" & 
                                         missing_loc$Annotation_Location == "No annotation" &
                                         missing_loc$Reads_Location != "No reads")
length(no_assembly_or_annot_but_reads)
write.table(unlist(strsplit(missing_loc$Reads_Location[no_assembly_or_annot_but_reads], split = ",")), 
            file = "missing/run_assembly_prokka.txt", col.names = F,
            row.names = F, quote = F)

## can't help with these
nothing = which(missing_loc$Assembly_Location == "No assembly" & 
                             missing_loc$Annotation_Location == "No annotation" &
                             missing_loc$Reads_Location == "No reads")
print(length(nothing))


