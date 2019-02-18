library(ggplot2)

## summarise platforms to know what to trim and how

setwd("/Users/gh11/e_colis/genomes/READS/")


platforms = read.table("sequencing_platforms.csv", stringsAsFactors = F,
                       header = F, comment.char = "", sep = ",")

## what platforms are in my list
unique_platforms = unique(platforms$V2)


## mostly illumina
## a few minION
## some from publications, for Hazen 2013 and Luo 2011 I can't say for sure what the sequencing platform was
# I might remove them from trimming (Luo 2011 are environmental which I decided not to use anyway)

## let's see how many of each
counts = data.frame(table(platforms$V2)) ## in total I need to remove 6 sequences from my collection and 3 mioion from trimming

counts$Var1 = factor(counts$Var1, counts$Var1[order(counts$Freq, decreasing = T)])
ggplot(counts, aes(x = Var1, y = Freq)) + geom_bar(stat = "identity")


## save a file on which ones not to trim
## only 13 marked for no trimming, the rest can be trimmed as usual
dont_trim = platforms$V1[which(platforms$V2 %in%
                                 c("hazen2013","MinION","luo2011"))]
write.table(dont_trim, "dont_trim.txt", quote = F, col.names = F, row.names = F)
