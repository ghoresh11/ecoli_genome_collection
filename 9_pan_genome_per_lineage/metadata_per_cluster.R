library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(gridExtra)
library(ggpubr)

setwd("/Users/gh11/poppunk_pangenome/2_dists_roary_analysis/")

metadata = read.table("metadata_per_cluster.csv", sep = ",",
                      header = T, stringsAsFactors = F, comment.char = "")

cluster_order = read.table("cluster_sizes_updated.csv", sep = ",", header = T, stringsAsFactors = F)
sizes = rep(0, dim(metadata)[1])
for (i in 1:length(sizes)){
  sizes[i] = cluster_order$Size[which(cluster_order$Cluster == metadata$cluster[i])]
}
metadata = cbind(metadata, size = sizes, orig = sizes * metadata$count)

cluster_sizes = read.table("cluster_sizes_updated.csv", sep = ",", header = T, stringsAsFactors = F)
total_count = sum(cluster_order$Size)

cluster_order = cluster_order$Cluster
clusters = unique(metadata$cluster)


md = read.table("/Users/gh11/e_colis/UPDATED_FINAL_METADATA_CLEANED.csv",sep = "\t", header = T, stringsAsFactors = F, comment.char = "", quote = "")
md_filtered = read.table("/Users/gh11/e_colis/FILTERED_MD_FINAL_ALL.tab",sep = "\t", header = T, stringsAsFactors = F, comment.char = "", quote = "")


### ST
ST = metadata[which(metadata$variable == "ST"),]
## count = data.frame(table(ST$value)) ## other than ST131,10,59 -> each ST is confined to one poppunk cluster

plot_piechart_for_ST <- function(c,min, title){
  curr_df = ST[which(ST$cluster == c),]
  #return(curr_df$value[which(curr_df$count == max(curr_df$count))])
  #  curr_df$cluster = as.character(curr_df$cluster)
  small = which(curr_df$count < min)
  if (length(small) > 0){
    curr_df = curr_df[-small,]
    curr_df = rbind(curr_df,
                    data.frame(cluster = c, variable = "ST", count = 1 - sum(curr_df$count), value =  "Other",
                               stringsAsFactors = F))
  }
  o =  curr_df$value[order(curr_df$count, decreasing = T)]
  if ("Other" %in% o){
    o = o[-which(o == "Other")]
    o = c(o, "Other")
  }
  curr_df$value = factor(curr_df$value,o)
  num_cols = length(unique(curr_df$value))
  colors = brewer.pal(n = 8, "Set2")[1:num_cols]
  colors = c(colors, rep("white", dim(curr_df)[1]))
  p =   ggplot(curr_df, aes(x = cluster, y = count, fill = value)) +
    geom_bar(stat = "identity", width = 1, color = "black")+ coord_polar(theta="y") +
    theme_minimal(base_size = 12) +
    xlab("") + ylab("")+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid  = element_blank()) +
    ggtitle(paste("Cluster:",c)) +
    scale_fill_manual(values = colors, name = "")+ theme(legend.position="bottom")  +
    guides(fill=guide_legend(nrow = 4, byrow=T)) + ggtitle(title) +
    labs(subtitle = paste("Cluster:", c))
  
  return(p)
}

pie1 = plot_piechart_for_ST(1, 0.001, "E") ## FOR THIS ONE ADD ST AS NAME
pie2 = plot_piechart_for_ST(2, 0.001, "F")
pie3 = plot_piechart_for_ST(3, 0.0025, "G")
pie4 = plot_piechart_for_ST(17, 0.0025, "H")
pie5 = plot_piechart_for_ST(40, 0.0025, "I")

plot_piechart_for_ST(12, 0.001, "hi")
### TO make a file of the most common ST in each cluster (can also look at the porportion)
# df = data.frame(i = 1:51, ans = rep("", 51), stringsAsFactors = F)
# for (i in 1:51){
#   df$ans[i] =  plot_piechart_for_ST(i, 0.0025, "I")
# }
# write.table(x = df, file = "/Users/gh11/Desktop/common_sts.csv", row.names = F, col.names = F, quote = F, sep=",")
# plot_piechart_for_ST(51)

# 
## MASH
# MASH = metadata[which(metadata$variable == "MASH"),]
# MASH$cluster = factor(MASH$cluster, cluster_order)
# ggplot(MASH, aes( x = cluster, y = count, fill = value)) + geom_bar(stat = "identity", color = "black") +
#   theme_bw(base_size = 16) +
#   xlab("Cluster") + ylab("% of isolates")

outpath = "/Users/gh11/Submissions/my_thesis/Chapter3/figures/4_clusters_metadata/"

plot_metadata_per_cluster<- function(name, order, cols, title){
  Continent = metadata[which(metadata$variable == name),]
  Continent$cluster = factor(Continent$cluster, cluster_order)
  Continent$value = factor(Continent$value, order)
  p = ggplot(Continent, aes( x = cluster, y = count, fill = value)) + geom_bar(stat = "identity", color = "black") +
    theme_classic(base_size = 12) + scale_fill_manual(values = cols, name = name)+
    xlab("Lineage") + ylab("Fraction of isolates") + scale_y_continuous(expand = c(0,0)) + ggtitle(title)+ 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme(legend.position = "bottom")
  
  return(p)
}


## Try to understand the proportion of isolates remaining in my dataset compared to how much there was originally
compare_filtered_to_complete <- function(name){
  j = which(colnames(md) == name)
  isolation_orig = data.frame(table(md[,j]))
  isolation_filtered = data.frame(table(md_filtered[,j]))
  isolation_filtered = isolation_filtered[match( isolation_orig$Var1, isolation_filtered$Var1),]
  isolation_orig = cbind(isolation_orig, filtered = isolation_filtered$Freq)
  isolation_orig$filtered[is.na(isolation_orig$filtered)] = 0
  isolation_orig = cbind(isolation_orig, diff = as.numeric(isolation_orig$Freq) - as.numeric(isolation_orig$filtered))
  isolation_orig = cbind(isolation_orig, percent = isolation_orig$diff / as.numeric(isolation_orig$Freq))
  return(isolation_orig)
}

statistically_test <- function(name, significance = 0.05){
  column = which(colnames(md_filtered) == name)
  uniques = unique(metadata$value[metadata$variable == name])
  poppunk_clusters = unique(md_filtered$Poppunk_cluster)
  num_tests = length(uniques)*length(poppunk_clusters)
  significance = significance / num_tests
  res = data.frame(cluster = character(0), value = character(0), pval = numeric(0), stringsAsFactors = F)
  for (u in uniques) {
    for (p in poppunk_clusters){
      all_vals = metadata[(metadata$variable == name & metadata$value == u),]
      all_true = sum(all_vals$orig)
      all_false = total_count - all_true
      curr_true = all_vals$orig[all_vals$cluster == p]
      curr_false = cluster_sizes$Size[cluster_sizes$Cluster == p] - curr_true
      curr = matrix(data = c(all_true, all_false, curr_true, curr_false), nrow = 2, ncol = 2)
      row.names(curr) = c("Africa", "Not Africa")
      colnames(curr) = c("All", "34")
      
     # if (dim(curr)[1] != 2 || dim(curr)[2] != 2) { next }
      test = fisher.test(curr, alternative = "less")
      
      if (test$p.value <= significance) {
       res = rbind(res, 
                   data.frame(cluster = p, value = u, pval = test$p.value, stringsAsFactors = F))
      }
    }
  }
  return(res)
}

path_order = c("nd", "epec", "etec","ehec","stec","expec", "eaec","epec/eaec","commensal")
path_cols =  c("#dddddd", brewer.pal(n = 7, "Set2"), "#4c4cff")
plot_metadata_per_cluster("Pathotype", path_order, path_cols, "")
path_comp = compare_filtered_to_complete("Pathotype")
path_stats = statistically_test("Pathotype")

iso_order = rev(c("feces","blood","urine","other/unknown"))
iso_cols = rev( c(brewer.pal(n = 3, "Set2"),"#dddddd"))
B= plot_metadata_per_cluster("Isolation", iso_order, iso_cols, "A")
iso_comp = compare_filtered_to_complete("Isolation")
iso_stat = statistically_test("Isolation")


cont_order = c("Europe", "North America","Africa","Asia","Oceania", "South America","nd")
cont_cols = c( brewer.pal(n = 6, "Set2"), "#dddddd")
C = plot_metadata_per_cluster("Continent", cont_order, cont_cols, "B")
cont_comp = compare_filtered_to_complete("Continent")
cont_stat = statistically_test("Continent")

pub_order = data.frame(table(metadata$value[metadata$variable == "Publication"]))
pub_order = as.character(pub_order$Var1[order(pub_order$Freq, decreasing = T)])
pub_cols = c(brewer.pal(n=8, "Dark2"), brewer.pal(n = 8, "Set3"), brewer.pal(n = 8, "Set1"))
plot_metadata_per_cluster("Publication", pub_order, pub_cols, "")
pub_comp = compare_filtered_to_complete("Publication", unique(md$Publication), rainbow(n = length(unique(md$Publication))))
pub_stat = statistically_test("Publication")



#### YEAR ###
Year = metadata[which(metadata$variable == "Year" & metadata$value != "nd"),]
length(which(Year$value>2010))
Year$cluster = factor(Year$cluster, cluster_order)
Year$value = as.numeric(Year$value)
graphics = read.table("/Users/gh11/Submissions/my_thesis/Chapter4//figures/cluster_graphics.csv", sep = ",",
                      header = T, comment.char = "")
graphics = graphics[match(cluster_order, graphics$Cluster),]

D = ggplot(Year, aes(x = value, y = cluster, size = count)) + geom_point(fill =  "#a9a9a9", pch = 21, color = "black") +
  scale_size_continuous(range = c(1,8), name = "Fraction of Isolates", breaks = seq(from=0.2, to =1, by =0.2)) +
  xlab("Year") + ylab("Lineage") + theme_bw(base_size = 12) + ggtitle("C")
# 
# D = ggplot(Year, aes(x = value, y = count, color = cluster, shape = cluster)) + geom_point(size = 3, alpha = 0.8, stroke = 1)+ 
#   scale_color_manual(values = as.character(graphics$Color), name = "") + scale_shape_manual(values = graphics$Shape, name = "") +
#   theme_bw(base_size = 12) + xlab("Year") + ylab("Fraction of isolates")+  
#   guides(color=guide_legend(ncol=7), shape =guide_legend(ncol = 7)) + ggtitle("D") +guides(color=guide_legend(nrow=6,byrow=TRUE))

#ggsave(plot =  p, file = paste(outpath, "Year.pdf", sep = ""), height = 5, width = 7)
# 
# legend = as_ggplot(get_legend(D))
# legend
# D = D + theme(legend.position = "None")

all_years = Year$value
year_stat = data.frame(cluster = unique(Year$cluster), pval = rep(0, length(unique(Year$cluster))), stringsAsFactors = F)
signif = 0.05 / dim(year_stat)[1]
for (i in 1:dim(year_stat)[1]){
  curr_cluster = year_stat$clusteruster[i]
  years_of_cluser = Year$value[Year$cluster == curr_cluster]
  pval = wilcox.test(all_years, years_of_cluser)$p.value
  year_stat$pval[i] = pval
}
year_stat = year_stat[which(year_stat$pval < 0.05),]

plot(hist(all_years))
plot(hist(years_of_cluser))
wilcox.test(all_years, years_of_cluser)

lay = rbind(c(NA,NA,NA,NA,3,3),
            c(NA,NA,NA,NA,3,3),
            c(1,1,1,1,1,1),
            c(2,2,2,2,2,2))

grid.arrange(B,C,D, layout_matrix = lay)
## to save as an image: width = 800, height = 800
##length(which(md$Poppunk_cluster == "2"))

grid.arrange(B + coord_flip() + theme(legend.position = "None"),
             C + coord_flip()+ theme(legend.position = "None"),  D +  theme(legend.position = "None"),
             layout_matrix = rbind(c(1,2,3,3)))

grid.arrange(as_ggplot(get_legend(B)),as_ggplot(get_legend(C)), as_ggplot(get_legend(D + theme(legend.position = "bottom"))))

### save 1000x1300
# ### Look at the metadata per cluster, stratified
# 
# loc = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/dists_analysis/metadata_within_cluster"
# plots = paste(loc, "plots/", sep = "/")
# #loc = "/Users/gh11/poppunk_pangenome/dists_analysis/"
# #plots = loc
# setwd(loc)
# 
# num_tests = 546 ## I had a look how many files were created
# required_pval = 0.05/ num_tests ## bonferroni corrected
# ## The hypothesis is that x<y (i.e. distance within a group is smaller than between groups)
# 
# for (cluster in 1:39){
#   curr_cluster = read.table(paste(cluster,"_metadata_within_cluster.csv", sep = ""), sep = ",",
#                             stringsAsFactors = F, header = T, comment.char = "", quote = "")
#   same = rep("Different", dim(curr_cluster)[1])
#   curr_cluster = cbind(curr_cluster, same)
#   columns = c(4:10)
#   for (c in columns){
#     curr_cluster$same = rep("Different", dim(curr_cluster)[1])
#     curr_cluster$same[which(curr_cluster[,c] == curr_cluster[,c+7])] = "Same"
#     title = paste("Cluster:", cluster )
#     name = strsplit(x = colnames(curr_cluster)[c], split = "_", fixed = T)[[1]][1]
#     
#     if (length(which(curr_cluster$same == "Different")) == 0  || length(which(curr_cluster$same == "Same")) == 0 ) {
#       next
#     }
#     
#     core_test = wilcox.test(x = curr_cluster$Core_dist[which(curr_cluster$same == "Same")],
#                             y = curr_cluster$Core_dist[which(curr_cluster$same == "Different")], alternative = "less")
#     
#     acc_test = wilcox.test(x = curr_cluster$Acc_dist[which(curr_cluster$same == "Same")],
#                            y = curr_cluster$Acc_dist[which(curr_cluster$same == "Different")], alternative = "less")
#     
#     if (core_test$p.value < required_pval || acc_test$p.value < required_pval) {
#       p = ggplot(curr_cluster, aes(x = Core_dist, y = Acc_dist, color = same)) + 
#         geom_point( size = 3, alpha = 0.7) +
#         theme_bw(base_size = 16) + 
#         scale_color_manual(values = c("#909090","#7c26cb")) +
#         xlab("Core distance") + ylab("Accessory distance") + ggtitle(paste(title, name, sep = "\n")) 
#       ggsave(p, height = 5, width = 6,
#                filename = paste(plots, cluster, "_", name, ".pdf", sep = ""))
#       
#       p = ggplot(curr_cluster, aes(x = same, y = Core_dist)) + geom_violin() +
#         geom_boxplot(width = 0.1)  +  ggtitle(paste(title, "\n",name, "\n",  core_test$p.value , sep = "")) +
#         ylab("Core distance") + xlab("") + theme_classic(base_size = 16)
#       ggsave(p, height = 4, width = 4, 
#              filename = paste(plots, cluster, "_", name, "_core_dists.pdf", sep = ""))
#       p = ggplot(curr_cluster, aes(x = same, y = Acc_dist)) + geom_violin() +
#         geom_boxplot(width = 0.1)  + ggtitle(paste(title, "\n",name, "\n",  acc_test$p.value , sep = "")) +
#         ylab("Accessory distance") + xlab("") + theme_classic(base_size = 16)
#       ggsave(p, height = 4, width = 4, 
#              filename = paste(plots, cluster, "_", name, "_acc_dists.pdf", sep = ""))
#     }
#   }
# }
# 
# 
# df = data.frame(gene_systems = c(1, 60000), genomes = c(300, 10000), type = c("ta", "ecoli"))
# ggplot(df, aes(x = gene_systems, y = genomes, color = type)) + geom_point(size = 6, pch = 13, stroke = 2)  + 
#   scale_y_continuous(limits = c(0, 10000)) + theme_bw(base_size = 20) + scale_alpha_continuous(limits = c(0, 60000)) + 
#   xlab("Genes") + ylab("Genomes") + scale_color_manual(values = c("red", "#d3d3d3"), guide = F)
