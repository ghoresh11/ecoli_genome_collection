library(ggplot2)
library(RColorBrewer)



setwd("/Users/gh11/e_colis/genomes/5_QC/")

metadata = read.table("beep_mash/FINAL_METADATA_MASH.csv", sep = "\t",
                      stringsAsFactors = F, comment.char = "", header = T)


qc = read.table("report.csv", sep = ",", header = T, 
                stringsAsFactors = F, comment.char = "")

## Cutoffs for fail/pass according to pathpipelines
kraken="?"
mapped_bases=80 # minimum 80%
error_rate=0.02 # smaller than 0.02
indel_ratio=0.02 # also smaller than 0.02
insert_size=80 #should be less than 80%
insert_size_rev=20 # not sure


for (i in 2:dim(qc)[2]){
  qc[,i] = as.numeric(qc[,i], na.rm = T)
  plot(density(qc[,i], na.rm = T), main = colnames(qc)[i], xlab =  colnames(qc)[i])
}

# for correlations between the different measures
for (i in 2:(dim(qc)[2]-1)){
  for (j in (i+1):dim(qc)[2]) {
    plot(qc[,i],qc[,j],
         xlab = colnames(qc)[i], ylab = colnames(qc)[j])
  }
}


plot(qc$mapped_bases, qc$kraken)
abline(v = 60, col = "red")
abline(h = 30, col = "red")

plot(qc$mapped_bases, qc$error_rate)
abline(v = 60, col = "red")
abline(h = 0.03, col = "red")

qc$indel_ratio[which(qc$indel_ratio>1)] = 1
plot(qc$error_rate, qc$indel_ratio)
abline(v = 60, col = "red")
abline(h = 0.03, col = "red")


length(which(qc$heterozygous_snps > 2))
## altogether: there are approximately 1,172 sequences that should be removed for low quality
remove = qc$ID[which(qc$kraken < 30 | qc$mapped_bases < 60 |
               qc$error_rate > 0.03 |
               qc$indel_ratio > 0.03 |
               qc$insert_size > 80 | 
               qc$heterozygous_snps > 3)]


### add st and MASH data to the QC dataframe
# qc = data.frame(qc, mash_cluster = rep(NA, dim(qc)[1]), st = rep(NA, dim(qc)[1]))
# for (i in 1:dim(metadata)[1]){
#   index = which(metadata$Run_ID == qc$ID[i])
#   if (length(index) == 0){
#     next
#   }
#   qc$mash_cluster[i] = metadata$MASH[index[1]]
#   qc$st[i] = metadata$ST[index[1]]
#   
# }
## The NA values are ones that have already been removed due to
## number of contigs and genome length!




#### FUNCTIONS TO REMOVE QCed GENOMES FROM DATAFRAME ####
split <- function(x){
  return(unlist(strsplit(x = x, split = ",", fixed = T)))
}

remove_from_vec <- function(run_id, vec){
  new_vec = c()
  for (item in vec){
    if (is.na(item) || item == "NA") {
      next
    }
    if (! grepl(pattern = tolower(run_id), x = tolower(item))){
      new_vec = c(new_vec, item)
    }
  }
  return(new_vec)
}

create_new_output <- function(vec) {
  if (length(vec) == 0) {
    return(NA)
  } 
  return(paste(vec, sep = "", collapse = ","))
}

remove_run <- function(metadata, colname, row) {
  col = which(colnames(metadata) == colname)
  vars = split(metadata[row, col])
  vars = remove_from_vec(run, vars)
  metadata[row,col] = create_new_output(vars) 
  return(metadata)
} 


### go over the the metadata file (all the READ locations)
## remove any READ/Assemblies/annotations associated with a bad quality run ID
for (run in remove){
  row = which(grepl(pattern = run, x= metadata$Run_ID))
  metadata = remove_run(metadata, "Assembly_Location", row)
  metadata = remove_run(metadata, "Annotation_Location", row)
  metadata = remove_run(metadata, "Reads_Location", row)

} 


## find removed - there are 795 which have NA now in all their fields
index_to_remove = which(is.na(metadata$Reads_Location)) ## remove anything where the reads were bad
index_to_remove = c(index_to_remove, which(metadata$Source != "human"))
## removed 1,460 genomes due to the QC reports

metadata_final = metadata[-index_to_remove,]
write.table(file = "../../FINAL_METADATA_MASH_CLEANED.csv",
            x = metadata_final, sep = "\t", quote = F, col.names = T, row.names = F)
