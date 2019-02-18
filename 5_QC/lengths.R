library(ggplot2)
library(RColorBrewer)

setwd("/Users/gh11/e_colis/genomes/5_QC/")

qc = read.table("lengths.csv", sep = ",", header = T, 
                stringsAsFactors = F, comment.char = "")

qc$length = qc$length/1000000

ggplot(qc, aes(x = length)) + geom_histogram(bins = 70) +
  theme_classic(base_size = 16) + xlab("Length (MBP)") + ylab("Count") +
  geom_vline(xintercept=c(4, 6), color = "#de4fe7", lwd = 1)

length(which(qc$num_contigs > 600))

ggplot(qc, aes(x = num_contigs)) + geom_histogram(bins = 70) +
  theme_classic(base_size = 16) + xlab("Length (MBP)") + ylab("Count")+ 
  geom_vline(xintercept=600, color = "#de4fe7", lwd = 1)

metadata = read.table("/Users/gh11/e_colis/genomes/final_metadata_with_loc.csv", sep = "\t", header = T, 
                      stringsAsFactors = F, comment.char = "")

qc = qc[match(metadata$ID,qc$ID),]
metadata = data.frame(metadata, length = qc$length, num_contigs = qc$num_contigs)
## genomes that need to be removed because of length or number of contigs
rm = which(qc$num_contigs > 600 | qc$length<4 | qc$length>6)
metadata = metadata[-rm,]

write.table(file = "final_metadata_post_QC.csv", x = metadata,
            col.names = T, row.names = F, quote = F, sep = "\t")
### This is how many were removed for length and number of contigs