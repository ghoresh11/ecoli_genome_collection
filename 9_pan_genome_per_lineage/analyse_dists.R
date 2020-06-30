library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(ggpubr)

setwd("/Users/gh11/poppunk_pangenome/2_dists_roary_analysis//")


cluster_sizes = read.table("cluster_sizes_updated.csv", sep = ",",
                           header = T, stringsAsFactors = F)

## draw a point plot with the cluster sizes
cluster_sizes = cluster_sizes[order(cluster_sizes$Size, decreasing = T),]
cluster_sizes = cbind(ID = 1:dim(cluster_sizes)[1], cluster_sizes, Cumsum = cumsum(cluster_sizes$Size))

ggplot(cluster_sizes, aes(x = ID ,y = cluster_sizes$Cumsum)) + geom_point(alpha = 0.5, size = 3) +
  xlab("poppunk cluster") + ylab("Total number of genomes") +
  scale_y_continuous(limits = c(0, 13500)) + theme_bw(base_size = 16)

temp = cluster_sizes[which(cluster_sizes$Size>20),]
cluster_order = temp$Cluster[order(as.numeric(temp$Size), decreasing = T)]
temp$Cluster = factor(temp$Cluster, cluster_order)
ggplot(temp, aes(x = Cluster, y = Size)) + geom_bar(stat = "identity") +
  theme_classic(base_size = 18)  +ylab("Genomes") + xlab("Group")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



## this goes up to 519, because from 520 the clusters are of size 1 so there's no distance
within_distances = read.table("all/within_cluster_dist.csv", sep = ",",
                           header = T, stringsAsFactors = F, comment.char = "")


between_distances = read.table("all/between_cluster_dist.csv", sep = ",",
                               header = T, stringsAsFactors = F, comment.char = "")

## analysing cluster 12
cluster_12 = between_distances[which(between_distances$cluster1 == 12 | between_distances$cluster2 == 12),]
cluster_12 = between_distances[which(between_distances$cluster1 %in% c(12,28,35) & between_distances$cluster2 %in% c(12,28,35)),]



chosen = cluster_sizes$Cluster[which(cluster_sizes$Size >= 20)]

## add a color specifying which clusters I kept in the analysis
chosen_within = rep("no",dim(within_distances)[1])
chosen_within[which(within_distances$Cluster %in% chosen)] = "yes"
within_distances = cbind(within_distances, chosen = chosen_within)

chosen_between = rep("no",dim(between_distances)[1])
chosen_between[which(between_distances$cluster1 %in% chosen & between_distances$cluster2 %in% chosen)] = "yes"
between_distances = cbind(between_distances, chosen = chosen_between)

### see the range of distance as boxplots
cols = brewer.pal(n = 3, "Greys")[c(2,3)]
within_core = melt(within_distances[,c(1:4, 8)], id.vars = c("Cluster", "chosen"))
ggviolin(within_core, "variable", "value", fill = "chosen", 
         palette =cols, add = "boxplot", xlab = "", ylab = "Distance")
within_acc =  melt(within_distances[,c(1, 5:8)], id.vars = c("Cluster", "chosen"))
ggviolin(within_acc, "variable", "value", fill = "chosen", 
         palette =cols, add = "boxplot", xlab = "", ylab = "Distance")
### As density plots
ggplot(within_distances, aes(x = Core_median, fill = chosen)) + geom_density(alpha = 0.6) +
  theme_classic(base_size = 16) + scale_fill_manual(values = cols) + 
  scale_y_continuous(expand = c(0,0)) + ylab("Density") + xlab("Within cluster mean core distance")
ggplot(within_distances, aes(x = Acc_median, fill = chosen)) + geom_density(alpha = 0.6) +
  theme_classic(base_size = 16) + scale_fill_manual(values = cols) + 
  scale_y_continuous(expand = c(0,0)) + ylab("Density") + xlab("Within cluster mean accessory distance")


### As boxplots
between_core = melt(between_distances[,c(1:5, 9)], id.vars = c("cluster1", "cluster2", "chosen"))
ggviolin(between_core, "variable", "value", fill = "chosen", 
         palette =cols, add = "boxplot", xlab = "", ylab = "Distance")
between_acc =  melt(between_distances[,c(1,2, 6:9)], id.vars = c("cluster1", "cluster2", "chosen"))
ggviolin(between_acc, "variable", "value", fill = "chosen", 
         palette =cols, add = "boxplot", xlab = "", ylab = "Distance")

# As density/histograms plots (copy paste from here or abov, but the density plots are NON misleading)
ggplot(between_distances, aes(x = core, fill = chosen)) + geom_density(alpha = 0.6) +
  theme_classic(base_size = 16) + scale_fill_manual(values =cols) + 
  scale_y_continuous(expand = c(0,0)) + ylab("Density") + xlab("Between cluster median core distance")
ggplot(between_distances, aes(x = accessory, fill = chosen)) + geom_density(alpha = 0.6) +
  theme_classic(base_size = 16) + scale_fill_manual(values =cols) + 
  scale_y_continuous(expand = c(0,0)) + ylab("Density") + xlab("Between cluster median accessory distance")

### get general values
max(within_distances$Core_max)
max(within_distances$Acc_max)
max(between_distances$core_max)
max(between_distances$accessory_max)


## Does the size of the cluster affect the within core/acc distances?
## I have this again below in a better form
within_distances = cbind(within_distances, size = cluster_sizes$Size[match(within_distances$Cluster, cluster_sizes$Cluster)])
ggplot(within_distances, aes(y = Core_median, x = size)) + geom_point(alpha = 0.4, size = 2) +
  theme_bw(base_size = 16) + xlab("Cluster size") + ylab("Within cluster core distance")

ggplot(within_distances, aes(y = Acc_median, x = size)) + geom_point(alpha = 0.4, size = 2) +
  theme_bw(base_size = 16) + xlab("Cluster size") + ylab("Within cluster accessory distance")


## It's too big and not important enough to check for the between cluster distances
## (Does having a bigger cluster mean you're more different from all the rest? NO (if run this))
# between_distances = cbind(between_distances, size1 = cluster_sizes$Size[match(between_distances$cluster1, cluster_sizes$Cluster)],
#              size2 = cluster_sizes$Size[match(between_distances$cluster2, cluster_sizes$Cluster)])
# ggplot(between_distances, aes(y = core_median, x = size1)) + geom_point(alpha = 0.4, size = 2) +
#   theme_bw(base_size = 16) + xlab("Cluster size") + ylab("Between cluster core distance")
# ggplot(between_distances, aes(y = core_median, x = size2)) + geom_point(alpha = 0.4, size = 2) +
#   theme_bw(base_size = 16) + xlab("Cluster size") + ylab("Between cluster core distance")

### Combining poppunk output with the roary outputs
variable_order = c("rare", "inter", "soft_core", "core")
variable_2_order = c("small", "not_small")
cols = brewer.pal(n = 5, "Blues")[-1]

roary_outputs = read.table("roary_summary_file.csv", sep = ",", header = T, stringsAsFactors = F)

## add info to the roary outputs CSV
sizes = cluster_sizes$Size[match(roary_outputs$cluster, cluster_sizes$Cluster)]
roary_outputs = cbind(roary_outputs, size = sizes)

core_dist = within_distances$Core_median[match(roary_outputs$cluster, within_distances$Cluster)]
core_max_dist = within_distances$Core_max[match(roary_outputs$cluster, within_distances$Cluster)]
acc_dist = within_distances$Acc_median[match(roary_outputs$cluster, within_distances$Cluster)]
acc_max_dist = within_distances$Acc_max[match(roary_outputs$cluster, within_distances$Cluster)]
roary_outputs = cbind(roary_outputs, core_dist, core_max_dist, acc_dist, acc_max_dist)


roary_outputs = cbind(roary_outputs,
                        size = roary_outputs$sizes[match(roary_outputs$cluster, roary_outputs$cluster)])



labs = c("Rare (present in fewer than <15% of isolates)", "Intermediate (present in 15% to 95%)", "Soft core (present in 95% to 99%)", "Core (present in >99% of isolates)")

#roary_outputs$cluster = factor(roary_outputs$cluster, cluster_order)



## boxplot
roary_outputs$variable = factor(roary_outputs$variable, rev(variable_order))
ggplot(roary_outputs, aes(x = variable, y = count, color = length, fill = variable))+
  geom_boxplot(width=0.4) +
  scale_color_manual(values = c("black", "#909090"), labels = c(">=100aa", "<100aa"), name = "Length") +
  theme_bw(base_size = 16) + ylab("Genes") + xlab("") +
  scale_x_discrete(labels = rev(labs)) +
  scale_fill_manual(values = rev(cols), guide = F)

## Collapse the length column
roary_output_collapsed = data.frame(stringsAsFactors = F)
for (i in unique(roary_outputs$cluster)){
  for (j in variable_order) {
    lines = roary_outputs[which(roary_outputs$cluster == i & roary_outputs$variable == j),]
    new_line = lines[1,]
    new_line$count = sum(lines$count)
    roary_output_collapsed = rbind(roary_output_collapsed, new_line)
  }
}
roary_outputs = roary_output_collapsed

## Barplot
cluster_order = unique(roary_outputs$cluster[order(roary_outputs$size, decreasing = T)])
roary_outputs$variable = factor(roary_outputs$variable, variable_order)
roary_outputs$cluster = factor(roary_outputs$cluster , graphics$Cluster)
C = ggplot(roary_outputs, aes(x = cluster, y = count, fill = variable)) + 
  geom_bar(stat = "identity", color = "black") + 
  scale_fill_manual(values = cols, labels = labs, name = "") + theme_classic(base_size = 12) +
  scale_y_continuous(expand = c(0,0)) + xlab("PopPUNK Cluster") + 
  ylab("Genes") + ggtitle("D") + theme(legend.position = "bottom") + guides(fill=guide_legend(nrow=3, byrow = T))
C


## General values to keep track of
median(roary_outputs$count[roary_outputs$variable == "core"])
sd(roary_outputs$count[roary_outputs$variable == "core"])
median(roary_outputs$count[roary_outputs$variable == "soft_core"])
sd(roary_outputs$count[roary_outputs$variable == "soft_core"])
median(roary_outputs$count[roary_outputs$variable == "inter"])
sd(roary_outputs$count[roary_outputs$variable == "inter"])
median(roary_outputs$count[roary_outputs$variable == "rare"])
sd(roary_outputs$count[roary_outputs$variable == "rare"])



## connection between cluster size and number of genes from each category
roary_outputs$variable = factor(roary_outputs$variable, variable_order)
roary_outputs$cluster = factor(roary_outputs$cluster , graphics$Cluster)
ggplot(roary_outputs, aes(y = count, x = size, fill = variable)) + 
  geom_point(size = 3, pch=21, color = "black", alpha = 0.7) +
  theme_bw(base_size = 16) + xlab("Cluster size") + ylab("Genes") +
  scale_fill_manual(values = cols, labels = labs)
D = ggplot(roary_outputs, aes(y = count, x = size, fill = variable)) + 
  geom_point(size = 3, pch=21, color = "black", alpha = 0.7) +
  theme_bw(base_size = 12) + xlab("Cluster size") + ylab("Genes") +
  scale_fill_manual(values = cols, guide = F) +  geom_smooth(method = "lm", aes(color = variable)) +
  scale_color_manual(values = cols, guide = F) + ggtitle("D") + scale_x_log10()
D

## add distances
ggplot(roary_outputs, aes(y = count, x = acc_dist, fill = variable)) + 
  geom_point(size = 3, pch=21, color = "black", alpha = 0.7) +
  theme_bw(base_size = 16) + xlab("Within cluster accessory distance") + ylab("Genes") +
  scale_fill_manual(values = cols, guide = F) 

ggplot(roary_outputs, aes(y = count, x = core_dist, fill = variable)) + geom_point(size = 3, pch=21, color = "black", alpha = 0.7) +
  theme_bw(base_size = 16) + xlab("Within cluster core distance") + ylab("Genes") +
  scale_fill_manual(values = cols, guide = F) 




### Piarwise distances between every two clusters
pairwise_core_matrix = matrix(0, nrow = 50, ncol = 50)
pairwise_acc_matrix = matrix(0, nrow = 50, ncol = 50)
for (k_i in 1:49) {
  for (k_j in (k_i+1):50) {
    i = cluster_sizes$Cluster[k_i]
    j = cluster_sizes$Cluster[k_j]
    index = which((between_distances$cluster1 == i & between_distances$cluster2 == j) |
             (between_distances$cluster1 == j & between_distances$cluster2 == i))
    pairwise_core_matrix[k_i,k_j] = between_distances$core_median[index]
    pairwise_core_matrix[k_j,k_i] = between_distances$core_median[index]
    pairwise_acc_matrix[k_i,k_j] = between_distances$accessory_median[index]
    pairwise_acc_matrix[k_j,k_i] = between_distances$accessory_median[index]
  }
}


get_lower_tri<-function(m){
  m[upper.tri(m)] <- NA
  return(m)
}

reorder_matrix <- function(m){
  hc <- hclust(dist(m))
  m <-m[hc$order, hc$order]
  return(m)
}

plot_matrix <- function(m, legend, names){
  colnames(m) = names
  rownames(m) = names
  m <- reorder_matrix(m)
  lower_tri = get_lower_tri(m)
#  lower_tri[which(lower_tri == 0)] = NA
  melted_m <- melt(lower_tri, na.rm = TRUE)
  melted_m$Var1 = factor(melted_m$Var1, melted_m$Var1[1:length(unique(melted_m$Var1))])
  melted_m$Var2 = factor(melted_m$Var2, melted_m$Var1[1:length(unique(melted_m$Var1))])
  p = ggplot(data = melted_m, aes(Var2, Var1, fill = value)) +
    geom_tile(color = "black")+
    scale_fill_gradient( high = "white", low = "blue", space = "Lab", name = legend) +
    theme_bw(base_size = 16)+
    coord_fixed() + xlab("Cluster") + ylab("Cluster") + theme(axis.text.x = element_text(angle = 90, vjust = 1))
  return (p)
}                                            

plot_matrix(pairwise_core_matrix, "Core", cluster_sizes$Cluster)
plot_matrix(pairwise_acc_matrix, "Accessory", cluster_sizes$Cluster)

plot(between_distances$core_median[chosen], between_distances$accessory_median[chosen])
