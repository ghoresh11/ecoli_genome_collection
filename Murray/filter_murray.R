setwd("/Users/gh11/e_colis/Murray/")

murray = read.table("metadata.csv", sep = ",", header = T, comment.char = "", quote = "", stringsAsFactors = F)
murray = murray[which(murray$WGS.genus == "Escherichia/Shigella"),]


murray2 = read.table("metadata2.txt", sep = "\t", header = T,
                     stringsAsFactors = F, comment.char = "", quote = "")

murray2 = murray2[match(murray$Murray.Collection.Identifier,murray2$Murray.Collection.Identifier),]

murray = cbind(murray, Date = murray2$Date.on..tube, Year = murray2$Year, Source =  murray2$Origin)

write.table(file = "metadata_complete.csv", x = murray,
            col.names = T, row.names = F, quote = F, sep = ",")
