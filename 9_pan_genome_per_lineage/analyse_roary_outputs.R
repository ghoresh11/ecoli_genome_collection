library(ggplot2)
library(RColorBrewer)
library(vegan)
library(ggpubr)
library(gridExtra)
library(micropan)

setwd("/Users/gh11/poppunk_pangenome/2_dists_roary_analysis//")


cluster_sizes = read.table("cluster_sizes_updated.csv", sep = ",",
                           header = T, stringsAsFactors = F)

## draw a point plot with the cluster sizes
cluster_sizes = cluster_sizes[order(cluster_sizes$Size, decreasing = T),]
cluster_sizes = cluster_sizes[cluster_sizes$Size > 20, ]

### gene size and number of genes
filtered_md = read.table("/Users/gh11/e_colis/FILTERED_MD_FINAL_ALL.tab", sep = "\t", 
                         comment.char = "", stringsAsFactors = F, quote = "", header = T)

## genome size
filtered_md$Poppunk_cluster = factor(filtered_md$Poppunk_cluster, cluster_sizes$Cluster)

graphics = read.table("/Users/gh11/Submissions/my_thesis/Chapter4//figures/cluster_graphics.csv", sep = ",",
                      comment.char = "", stringsAsFactors = F, header = T)
graphics = graphics[match(cluster_sizes$Cluster,graphics$Cluster),]

A = ggplot(filtered_md, aes( x = length, y = num_new_genes, color = Poppunk_cluster, shape = Poppunk_cluster))  + 
 geom_point(alpha = 0.6, size = 1.5, stroke = 1) +
  theme_bw(base_size = 12) + ylab("Genes") + xlab("Genome length (Mbp)") +
  scale_color_manual(values = graphics$Color, name = "") + 
  scale_shape_manual(values = graphics$Shape, name = "") + ggtitle("B") +
  guides(shape=guide_legend(nrow=5, byrow = T), color=guide_legend(nrow=5, byrow =T))
legendA = as_ggplot(get_legend(A))
legendA
A = A +  theme(legend.position="none")
A
legendA
## get the weighted mean of genome length and number of genes
means = aggregate(filtered_md[,36], list(filtered_md$Poppunk_cluster), median)
mean_num_genes = mean(means$x)
min(filtered_md$num_new_genes)
max(filtered_md$num_new_genes)

filtered_md$Poppunk_cluster = factor(filtered_md$Poppunk_cluster, graphics$Cluster)
A0 = ggplot(filtered_md, aes( x = Poppunk_cluster, y = num_new_genes))  + 
#  geom_jitter(height = 0, width = 0.01, alpha = 0.6, size = 1.5, stroke = 1) +
  theme_bw(base_size = 12) + ylab("Genes") + xlab("PopPUNK cluster") +
  scale_color_manual(values = graphics$Color, name = "", guide = F) + 
  scale_shape_manual(values = graphics$Shape, name = "", guide = F) + ggtitle("A") +
  geom_boxplot(color = "black", outlier.shape = NA, fill = NA) + geom_hline(yintercept = mean_num_genes, lty = 2, color = "#ff0080") + ggtitle("A")+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))





# ## gene accumilation curves using vegan
# df = data.frame(cluster = character(0),
#                 richness = numeric(0),
#                 genomes = numeric(0),
#                 sd = numeric(0))
# 
# files = list.files(path = "gene_presence_absences/", full.names = T, pattern = "*.Rtab")
# for (f in files){
#   cluster = strsplit(basename(f), split = "_", fixed = T)[[1]][1]
#   print(cluster)
#   mydata <- data.frame(t(read.table(f, header = T, row.names = 1, comment.char = "", quote = )))
#   sp <- specaccum(mydata, "random", permutations=100)
#   
#   df = rbind(df, data.frame(cluster = rep(cluster, length(sp$sites)),
#                                  richness = sp$richness,
#                                  genomes = sp$sites,
#                                  sd = sp$sd))
# }
## save the output

# ## estimate the slope?
# # ## gene accumilation curves using vegan
# df = data.frame(cluster = character(0),
#                 alpha = numeric(0),
#                 intercept = numeric(0), stringsAsFactors = F)
# # 
# files = list.files(path = "gene_presence_absences/", full.names = T, pattern = "*.Rtab")
# for (f in files){
#   cluster = strsplit(basename(f), split = "_", fixed = T)[[1]][1]
#   print(cluster)
#   mydata <- data.frame(t(read.table(f, header = T, row.names = 1, comment.char = "", quote = )))
#   h.est <- heaps(mydata)
# 
#   df = rbind(df, data.frame(cluster = cluster,
#                             alpha = h.est[[2]],
#                             intercept = h.est[[1]],
#                             stringsAsFactors = F))
# }
# ## save the heaps result
# write.table(x = df, file = "heaps.csv", col.names = T, row.names = F, quote = F, sep = ",")

#write.table(df, "accumilation_curves.csv", col.names = T, row.names = F, quote = F, sep = ',')
df = read.table("accumilation_curves.csv", header = T, comment.char = "", sep = ",", stringsAsFactors = F)
df = cbind(df, min= df$richness-df$sd, max = df$richness+df$sd )
df$cluster = factor(df$cluster,cluster_sizes$Cluster)
df = df[which(df$genomes %% 16 == 0),] ## to visualise more clearly

  B = ggplot(df, aes(x = genomes, y = richness, color = cluster, shape = cluster))+ geom_line(alpha = 0.8)+ geom_point(size = 2, alpha = 0.8) +
  theme_bw(base_size = 12) + xlab("Genomes") +
  ylab("Genes") +
  geom_errorbar(aes(ymin=min, ymax=max), width=.2,alpha = 0.2,
                position=position_dodge(0.05)) +
  scale_color_manual(values = graphics$Color, guide = F) + scale_shape_manual(values = graphics$Shape, guide = F) + ggtitle("C") +
    scale_x_continuous(limits = c(0,900)) 

### 
#write.table(tree, "../7_tree/tree.txt", col.names = F, row.names = F, quote = F, sep = ",")

phylogroup_shapes = read.table("/Users/gh11/Submissions/my_thesis/Chapter3/figures/phylogroup_shapes.csv",
                               sep = ",", header = T, stringsAsFactors = F)
plot.new()
legend("topleft", legend =  phylogroup_shapes$Phylogroup, pch = as.numeric(phylogroup_shapes$Shape))


## combine this with the plots in "analyse_dists.R"
lay = rbind(c(1,1,1,2,2,3,3),
            c(1,1,1,2,2,3,3),
            c(4,4,4,4,4,4,4),
            c(4,4,4,4,4,4,4),
            c(4,4,4,4,4,4,4))
D
grid.arrange(  B + ggtitle("B"), C, layout_matrix = lay)
## Add the legends separately...
## To save: height = 630 width = 880
# 
# 
# cluster_12 = read.table("/Users/gh11/e_colis/UPDATED_FINAL_METADATA_CLEANED.csv", sep = "\t",
#                         header = T, stringsAsFactors = F, comment.char = "", quote = "")
# cluster_12 = cluster_12[which(cluster_12$Poppunk_cluster == 35),]


### calculate whether some clusters have fewer or more genes relative to general population!
calculate_stats <- function(vec, md) {
  rs = sample(x  = unique(md$Poppunk_cluster),size = 10000, replace = T)
  res = rep(0, 10000)
  for (index in 1:10000) {
    cluster = rs[index]
    curr = vec[which(md$Poppunk_cluster == cluster)]
    res[index] = sample(size = 1, x = curr)
  } 
  print(boxplot(res))
  signif = c()
  for (cluster in unique(md$Poppunk_cluster)) {
    curr_test = wilcox.test(res, vec[which(md$Poppunk_cluster == cluster)])
    if (curr_test$p.value<(0.05/50)) {
      signif = c(signif, cluster)
    }
  }
  return(signif)
}
calculate_stats(filtered_md$num_new_genes, filtered_md)
