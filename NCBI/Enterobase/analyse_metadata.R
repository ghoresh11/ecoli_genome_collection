library(ggplot2)
library(reshape)


setwd("/Users/gh11/e_colis/genomes/Enterobase/")
metadata = read.table(file = "enterobase_fixed.txt", 
                      sep = "\t", header = T, comment.char = "", 
                      quote = "", stringsAsFactors = F)


plot_dat
a <- function(vec,name){
  t = melt(table(vec))

  
  if (name == "Year"){
    t = t[-which(t$vec == "nd"),]
    t$vec = as.numeric(as.character(unlist(t$vec)))
   # t = t[-which(t$vec < 1960),]
  } else {
    o = t$vec[order(t$value, decreasing = T)]
    t$vec = factor(t$vec, o)
  }

  
  p = ggplot(t, aes(x = vec, y = value)) + geom_bar(stat = "identity", fill = "#0b426b") +
    xlab("") + ylab("Number of genomes")+ theme_bw(base_size = 16) + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
    scale_y_continuous(expand = c(0,0)) + theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
    ggtitle(name)
  ggsave(filename = paste("figures/",name,".png",sep=""), dpi = 300,
         height = 5, width = max(0.5*dim(t)[1],4))
  return(p)
}

### strains which have no metadata associated with them at all
no_data = length(which(metadata$source.niche == "ND" &
                         metadata$source.type == "ND" &
                         metadata$collection.year == "nd" &
                         metadata$country == "nd" &
                         metadata$continent == "nd" &
                         metadata$disease == "nd" &
                         metadata$pathogen.non.pathogen == "nd"))
no_data / dim(metadata)[1]
## 31% of the sequences in enterobase have no metadata associated with them



### focus on human and environmental isolates?
plot_data(metadata$source.niche, "Source_Niche")
plot_data(metadata$source.type, "Source_Type")

metadata = metadata[which(metadata$source.niche == "human" |
                            metadata$source.niche == "environment"  |
                            grepl("ecor", metadata$name, ignore.case = T)),]


plot_data(metadata$collection.year, "Year")
plot_data(metadata$country, "Country")

metadata$continent[which(metadata$continent == "unresolved")] = "nd"
plot_data(metadata$continent, "Continent")
plot_data(metadata$disease, "Disease")
plot_data(metadata$pathogen.non.pathogen, "Pathogenicity")

labs = data.frame(table(metadata$lab.contact))
labs = labs[order(labs$Freq, decreasing = T),]
write.table(x = labs, file = "lab_contants.csv", sep = ",", quote = F,
          col.names = T, row.names = F)

labs = labs[which(labs$Freq > 50),]
labs$Var1 = factor(labs$Var1, labs$Var1)
ggplot(labs, aes(x = Var1, y = Freq)) + geom_bar(stat = "identity", fill = "#0b426b") +
         xlab("") + ylab("Number of genomes")+ theme_bw(base_size = 16) + 
         theme(panel.border = element_blank(), panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
         scale_y_continuous(expand = c(0,0)) + theme(axis.text.x = element_text(angle = 70, hjust = 1)) +
         ggtitle("Labs")


### check release dates
metadata$release.date[which(metadata$release.date == "nd")] = "2100-01-01"
metadata$release.date = as.Date(metadata$release.date , format = "%Y-%m-%d")
metadata = metadata[order(metadata$release.date, decreasing = T),]

length(metadata$release.date[which(metadata$release.date > as.Date("2019-01-01", format = "%Y-%m-%d"))])


qplot(metadata$contig.number...200.bp., geom = "histogram", bins = 100)+ theme_bw(base_size = 16) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  scale_y_continuous(expand = c(0,0)) + theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
  ggtitle("Number of contigs") + xlab("Number of contigs") + ylab("Number of genomes")
ggsave(filename = paste("figures/","num_contigs",".png",sep=""), dpi = 300,
       height = 5, width = 5)

qplot(metadata$length, geom = "histogram", bins = 100)+ theme_bw(base_size = 16) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  scale_y_continuous(expand = c(0,0)) + theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
  ggtitle("Genome Length") + xlab("Genome Length (bp)") + ylab("Number of genomes")
ggsave(filename = paste("figures/","genome_lengths",".png",sep=""), dpi = 300,
       height = 5, width = 5)
