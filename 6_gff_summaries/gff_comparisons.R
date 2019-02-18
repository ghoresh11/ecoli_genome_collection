library(ggplot2)
library(RColorBrewer)


setwd("/Users/gh11/e_colis/gff_summaries/")


metadata = read.table("../FINAL_METADATA_CLEANED.csv", sep = "\t", quote = "",
                      comment.char = "", header = T, stringsAsFactors = F)
gff_summary = read.table("gff_summary.csv", sep = ",", comment.char = "",
                         quote = "", header = T, stringsAsFactors = F)

gff_summary = gff_summary[match(metadata$ID, gff_summary$ID),]

gff_summary = data.frame(ID = metadata$ID, 
                         Pathotype = metadata$Pathotype,
                         Continent =  metadata$Continent, 
                         Year = metadata$Year, 
                         Isolation = metadata$Isolation,
                         MASH  = metadata$MASH, 
                         ST = metadata$ST, 
                         Length = metadata$length,
                         Contigs = metadata$num_contigs, 
                         gff_summary[,seq(3,15)],
                         stringsAsFactors = F)

common_STs = data.frame(table(gff_summary$ST), stringsAsFactors = F)
common_STs = common_STs[order(common_STs$Freq, decreasing = T),]
common_STs = common_STs$Var1[1:20]
common_STs = gff_summary[which(gff_summary$ST %in% common_STs),]

common_MASH = data.frame(table(gff_summary$MASH), stringsAsFactors = F)
common_MASH = common_MASH[order(common_MASH$Freq, decreasing = T),]
common_MASH = common_MASH$Var1[1:20]
common_MASH = gff_summary[which(gff_summary$MASH %in% common_MASH),]
common_MASH$MASH = as.character(common_MASH$MASH)


### Connection between the size of the genome and number of genes
## number of contigs and number of genes


### correlation between length and number of genes
a = lm(formula =cds~Length, data = gff_summary)
res = abs(residuals(a))
res[which(res > 500)] = "no_fit"
res[which(res != "no_fit")] = "fit"
test = data.frame(gff_summary, col = res, stringsAsFactors = F)
ggplot(test, aes(x = Length, y = cds, color = col)) + geom_point() + theme_bw(base_size = 20) +
  scale_color_manual(values = c("black", "#ff6666"), guide = F) + ylab("Number of CDSs") +
  xlab("Length (MBP)")

ggplot(test, aes(x = Contigs, y = cds, color = col)) + geom_point() + theme_bw(base_size = 20)+
  scale_color_manual(values = c("black", "#ff6666"), guide = F) + ylab("Number of CDSs") +
  xlab("Number of Contigs")


ggplot(test, aes(x = Length, y = Contigs, color = col)) + geom_point() + theme_bw(base_size = 20) + 
  scale_color_manual(values = c("black", "#ff6666"), guide = F) + xlab("Length (MBP)") +
  ylab("Number of Contigs")

### save the final metadata file without the weird genomes
final = metadata[which(test$col == "fit"),]
write.table("../VERY_FINAL_METADATA_CLEANED.csv", x = final, sep = "\t",
            quote = F, row.names = F, col.names = T)


## Connection between genome length and tyoes
ggplot(gff_summary, aes(x = Pathotype, y = Length)) + geom_violin() + theme_bw(base_size = 20)+
  stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.1 ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(gff_summary, aes(x = Continent, y = Length)) + geom_violin() + theme_bw(base_size = 20)+
  stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.1 ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))




ggplot(common_STs, aes(x = ST, y = Length)) + geom_violin() + theme_bw(base_size = 20)+
  stat_summary(fun.data="mean_sdl",  geom="crossbar", width=0.1 ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

## to go to apairwise wilcox test
pairwise.wilcox.test(as.numeric(common_STs$Length), common_STs$ST, paired = F,p.adjust.method = "fdr")


ggplot(common_MASH, aes(x = MASH, y = Length)) + geom_violin() + theme_bw(base_size = 20)+
  stat_summary(fun.data="mean_sdl",  geom="crossbar", width=0.1 ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
pairwise.wilcox.test(as.numeric(common_MASH$Length), common_MASH$MASH, paired = F,p.adjust.method = "fdr")

# 
# barplot_length_by_variable<- function(df){
#   length_bins = seq(min(df$Length),max(df$Length), by = 0.2)
#   length_bins = c(0, length_bins, 100)
#   vec = rep("",dim(df)[1])
#   for (i in 2:length(length_bins)){
#     vec[which(df$Length > length_bins[i-1] & df$Length <= length_bins[i])] = length_bins[i]
#   }
#   
#   ### barplot colored by a specific variable
#   df = cbind(df, length_bins = vec)
#   df2 = data.frame(table(df$length_bins, df$ST))
#   df2$Var1 = factor(df2$Var1, as.character(seq(min(df$Length),max(df$Length), by = 0.2)))
#   ggplot(df2, aes(x = Var1, y = Freq, fill = Var2)) + geom_bar(stat = "identity")
# }
# barplot_length_by_variable(common_STs)
# barplot_length_by_variable(common_MASH)


## How many genes does E. coli have in general?
ggplot(gff_summary, aes(x = cds)) + geom_histogram(bins = 40) + theme_bw(base_size = 20)

ggplot(gff_summary, aes(x = Pathotype, y = cds)) + geom_violin() + theme_bw(base_size = 20)+
  stat_summary(fun.data="mean_sdl", mult=1, geom="crossbar", width=0.1 ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(gff_summary, aes(x = Continent, y = cds)) + geom_violin() + theme_bw(base_size = 20)+
  stat_summary(fun.data="mean_sdl",  geom="crossbar", width=0.1 ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(common_STs, aes(x = ST, y = cds)) + geom_violin() + theme_bw(base_size = 20)+
  stat_summary(fun.data="mean_sdl",  geom="crossbar", width=0.1 ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(common_MASH, aes(x = MASH, y = cds)) + geom_violin() + theme_bw(base_size = 20)+
  stat_summary(fun.data="mean_sdl",  geom="crossbar", width=0.1 ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



### rather than looking at CDS -> look at the specific elements I dug out
ggplot(gff_summary, aes(x = phage)) + geom_histogram(bins = 40) + theme_bw(base_size = 20)

ggplot(gff_summary, aes(x = Pathotype, y = phage)) + geom_violin() + theme_bw(base_size = 20)+
  stat_summary(fun.data="mean_sdl", mult=1, geom="crossbar", width=0.1 ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(gff_summary, aes(x = Continent, y = phage)) + geom_violin() + theme_bw(base_size = 20)+
  stat_summary(fun.data="mean_sdl", mult=1, geom="crossbar", width=0.1 ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(common_STs, aes(x = ST, y = phage)) + geom_violin() + theme_bw(base_size = 20)+
  stat_summary(fun.data="mean_sdl",  geom="crossbar", width=0.1 ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(common_MASH, aes(x = MASH, y = phage)) + geom_violin() + theme_bw(base_size = 20)+
  stat_summary(fun.data="mean_sdl",  geom="crossbar", width=0.1 ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


### see if there is a correlation between the presence of two genes
for (i in 8:(dim(gff_summary)[2]-1)){
  for (j in (i+1):dim(gff_summary)[2]){
    correlation = cor(gff_summary[,i], gff_summary[,j])
    
    plot(gff_summary[,i],gff_summary[,j], pch = 16,
         xlab = colnames(gff_summary)[i], ylab = colnames(gff_summary)[j], 
         main = paste("Pearson:", correlation))
    
  }
}

## see if one group is enriched to having more of a specific kind of gene?
### PATHOTYPE
for (i in 8:(dim(gff_summary)[2])){
  boxplot(gff_summary[,i] ~ gff_summary$Pathotype, pch = 16,
          xlab = colnames(gff_summary)[i], ylab = "Pathotype",
          main = colnames(gff_summary)[i])
}

## see if one group is enriched to having more of a specific kind of gene?
### GEOGRAPHY
for (i in 8:(dim(gff_summary)[2])){
  boxplot(gff_summary[,i] ~ gff_summary$Continent, pch = 16,
          xlab = colnames(gff_summary)[i], ylab = "Continent",
          main = colnames(gff_summary)[i])
}

### MASH cluster
for (i in 8:(dim(common_MASH)[2])){
  boxplot(common_MASH[,i] ~ common_MASH$MASH, pch = 16,
          xlab = colnames(common_MASH)[i], ylab = "MASH cluster",
          main = colnames(common_MASH)[i])
}



### ST
for (i in 8:(dim(common_STs)[2])){
  boxplot(common_STs[,i] ~ common_STs$ST, pch = 16,
          xlab = colnames(common_STs)[i], ylab = "ST",
          main = colnames(common_STs)[i])
}
