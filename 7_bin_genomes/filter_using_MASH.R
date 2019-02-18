

setwd("/Users/gh11/e_colis/genomes/5_QC/beep_mash/")


## re-read the file
metadata = read.table("/Users/gh11/e_colis/mlst/metadata_w_mlst.csv", sep = "\t", header = T,
                      comment.char = "", stringsAsFactors = F)


# ## ADD MASH CLUSTER TO METADATA
mash_clusters = read.table("final_binned.txt", sep = ",",
                           header = F, comment.char = "", stringsAsFactors = F)

mash_to_metadata = rep(NA, dim(metadata)[1])
for (i in 1:dim(metadata)[1]) {
  assemblies = strsplit(metadata$Annotation_Location[i], split = ",", fixed = T)[[1]]
  for (a in assemblies){
    index = which(mash_clusters$V1 == a)
    if (length(index) > 0) {
      mash_to_metadata[i] = index
    }
  }
}

## This works -> the only ones that are missing are the ones that have been removed
## because they are contaminated!

## Remove contaminents/wrong length/num_contigs:
metadata = metadata[-which(is.na(mash_to_metadata)),]
mash_to_metadata = mash_to_metadata[-which(is.na(mash_to_metadata))]
mash_clusters = mash_clusters[mash_to_metadata, ]

metadata = data.frame(metadata, MASH = mash_clusters$V2)
write.table(x = metadata, file = "FINAL_METADATA_MASH.csv", sep = "\t", col.names = T,
            row.names = F, quote = F)


# ### Compare MASH and MLST
# metadata = read.table("FINAL_METADATA_MASH.csv", sep = "\t", header = T,
#                       comment.char = "", stringsAsFactors = F)
st_v_mash = data.frame(table(metadata$MASH, metadata$ST))
st_v_mash = st_v_mash[-which(st_v_mash$Freq == 0),]
colnames(st_v_mash) = c("MASH", "ST", "Freq")



## how many unique STs are there?
print("Number of unique MLSTs")
num_STs = length(unique(metadata$ST))
print(num_STs)

counts = sort(table(all_mlst$ST), decreasing = T)
counts = melt(counts)
counts$Var1 = factor(counts$Var1, counts$Var1)
counts$value = as.numeric(counts$value)
