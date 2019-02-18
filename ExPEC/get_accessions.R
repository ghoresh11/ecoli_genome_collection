

setwd("/Users/gh11/e_colis/genomes/ExPEC/")

chen2013 = read.table("Chen2013.csv", stringsAsFactors = F, sep = ",", header = T,
                      comment.char = "")
chen2013 = chen2013$SRA.accession
write.table(chen2013, file = "../READS/chen2013_reads.txt", row.names = F, col.names = F,
            quote = F)

others000 = read.table("other.csv", stringsAsFactors = F, sep = ",", header = T,
                       comment.char = "")
write.table(others000$Accession, file = "others0000_accessions.txt", row.names = F,
            col.names = F, quote = F)
