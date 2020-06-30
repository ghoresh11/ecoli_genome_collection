library(ggplot2)
library(RColorBrewer)
library(vegan)
library(gridExtra)
library(ggpubr)
library(ggtree)
library(reshape2)
library(ape)

setwd("/Users/gh11/poppunk_pangenome/4_pairwise_roary/")

variable_order = c("rare","inter", "core")


## graphics
graphics = read.table("/Users/gh11/Submissions/my_thesis/Chapter4//figures/cluster_graphics.csv", sep = ",",
                      header = T, comment.char = "", stringsAsFactors = F)
freqs = read.table("231019_corrected//freqs.csv", header = T,row.names = 1,
                   stringsAsFactors = F, comment.char = "", quote = "", sep =",")
freqs = freqs[,-which(colnames(freqs) == "X50")]
clusters = sapply(X = colnames(freqs), FUN = gsub, pattern = "X", replacement = "")
graphics = graphics[match(clusters, graphics$Cluster),]
o = clusters
phylo_order = c("B1","C", "Shigella","A","E","D","B2","F")

create_pca_plot <- function(curr_freqs){
  for_pca = t(curr_freqs)
  remove = c()
  for (i in 1:dim(for_pca)[2]) {
    if (length(unique(for_pca[,i])) == 1) { ## no variation in gene
      remove = c(remove, i)
    }
  }
  for_pca = for_pca[,-remove]
  ### PCA plot of the clusters -> what are the relationships between the clusters based on the frequencies of all genes
  freqs.pca = prcomp(for_pca , center = T)
  summary(freqs.pca)
  freqs.pca = data.frame(freqs.pca$x)
  freqs.pca = cbind(freqs.pca, Cluster = as.character(o))
  freqs.pca$Cluster = factor(freqs.pca$Cluster, o)
  return (freqs.pca)
}

freqs.pca = create_pca_plot(freqs)
## to focus only on particular gene
# counts = apply(FUN = function(v) length(which(v < 0.15)), 2, X = for_pca)
# for_pca = for_pca[,which(counts==47)]
freqs.pca$phylogroup = graphics$Phylogroup[match(freqs.pca$Cluster, graphics$Cluster)]


## PCA plot
C = ggplot(freqs.pca, aes(x = PC1, y = PC2, label = Cluster,colour = phylogroup)) + geom_text(size =6,nudge_x=0.5,nudge_y = 0.5 )+ geom_point()+
  scale_color_manual(values = graphics$Color, guide = F) + theme_bw(base_size = 12) +
  scale_shape_manual(values =  graphics$Shape, guide = F) + ggtitle("C")  +
  annotate("text", x = -12, y = 1, label = "B1", size = 5, parse = T) +
  annotate("text", x = 0.16, y = 0.02 , label = "A", size = 5, parse = T) + 
  annotate("text", x = -8, y = -8 , label = "E", size = 5,parse = T)+ 
  annotate("text", x = 8, y = -11 , label = "F", size = 5,parse = T)+ 
  annotate("text", x = -3, y = -13 , label = "D", size = 5,parse = T)+ 
  annotate("text", x = 10, y = 4 , label = "B2", size = 5,parse = T)
C

## for my interview
freqs.pca$phylogroup = graphics$Phylogroup[match(freqs.pca$Cluster, graphics$Cluster)]
C = ggplot(freqs.pca, aes(x = PC1, y = PC2, fill = phylogroup)) + geom_point(size = 3, pch = 21, color = 'black') +
  theme_bw(base_size = 14) + scale_fill_brewer(palette = "Set2", name = "Phylogroup") +
  theme(legend.position = "bottom")

# ## create PCA plots for genes which are core (i.e. core in 47 of 47) -> plotting according to core is meaningless because it's the same for all...
# keep = c()
# for (i in 1:dim(freqs)[1]) {
#   non_inter = length(which(freqs[i,] < 0.15 | freqs[i,] >= 0.9))
#   if (non_inter == 0) { keep = c(keep, i)}
# }
# curr_freqs = freqs[keep,]
# freqs.core.pca = create_pca_plot(curr_freqs)
# ggplot(freqs.core.pca, aes(x = PC1, y = PC2, color = Cluster, shape = Cluster)) + geom_point(size = 3.5, stroke = 1, alpha = 0.7) +
#   scale_color_manual(values = graphics$Color, guide = F) + theme_bw(base_size = 12) +
#   geom_text(aes(label=Cluster),hjust=-0.3, vjust=-0.3) +
#   scale_shape_manual(values =  graphics$Shape, guide = F) + ggtitle("C")
# 
# ## create PCA plots for genes which are rare (always under 15%)
# keep = c()
# for (i in 1:dim(freqs)[1]) {
#   non_rare = length(which(freqs[i,] > 0.15))
#   if (non_rare == 0) { keep = c(keep, i)}
# }
# curr_freqs = freqs[keep,]
# freqs.rare.pca = create_pca_plot(curr_freqs)
# ggplot(freqs.rare.pca, aes(x = PC4, y = PC3, color = Cluster, shape = Cluster)) + geom_point(size = 3.5, stroke = 1, alpha = 0.7) +
#   scale_color_manual(values = graphics$Color, guide = F) + theme_bw(base_size = 12) +
#   geom_text(aes(label=Cluster),hjust=-0.3, vjust=-0.3) +
#   scale_shape_manual(values =  graphics$Shape, guide = F) + ggtitle("C")
# 


# 
# # Trying to answer these questions:
# # 1. How many rare/core genes are shared?
# # 2. How many core genes are shared?
# gene_classes = data.frame(gene = rownames(freqs),
#                           core = rep(0, dim(freqs)[1]),
#                           inter = rep(0, dim(freqs)[1]),
#                           rare = rep(0, dim(freqs)[1]), stringsAsFactors = F)
# for (i in 1:dim(freqs)[1]){
#   curr_vec = freqs[i,]
#   gene_classes$core[i] = length(which(curr_vec >= 0.95))
#   gene_classes$inter[i] = length(which(curr_vec < 0.95 & curr_vec >= 0.15))
#   gene_classes$rare[i] =  length(which(curr_vec < 0.15 & curr_vec > 0))
# }
# write.table(x = gene_classes, file  = "231019_corrected///gene_classes.csv", sep = ",", col.names = T, row.names = F, quote = F)

gene_classes = read.table("231019_corrected//gene_classes.csv", sep = ",", header = T, stringsAsFactors = F, comment.char = "", quote = "", row.names = 1)
total_presence = data.frame(table(gene_classes$total_presence))

## examples of genes
gene_name = "intA_1" ## change here to check the frequency of a gene across the popPUNK clusters
curr = data.frame(name = o,
                  freq = unlist(freqs[which(rownames(freqs) == gene_name),]), stringsAsFactors = F)
curr = cbind(curr, phylo = graphics$Phylogroup[match(curr$name, graphics$Cluster)])
curr$name = factor(curr$name, o)
curr$phylo = factor(curr$phylo, phylo_order)
A = ggplot(curr, aes(x = name, y = freq)) + geom_bar(stat = "identity") +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) + theme_bw(base_size = 12) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Frequency\nin lineage") + xlab("Lineage") + 
  facet_grid(~phylo, scales = "free",  space = "free",switch = "x")  +
  ggtitle("A - intA")


gene_name = "wzyE" ## change here to check the frequency of a gene across the popPUNK clusters
curr = data.frame(name = o,
                  freq = unlist(freqs[which(rownames(freqs) == gene_name),]), stringsAsFactors = F)
curr = cbind(curr, phylo = graphics$Phylogroup[match(curr$name, graphics$Cluster)])
curr$name = factor(curr$name, o)
curr$phylo = factor(curr$phylo, phylo_order)
B = ggplot(curr, aes(x = name, y = freq)) + geom_bar(stat = "identity") +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) + theme_bw(base_size = 12) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Frequency\nin lineage") + xlab("Lineage") + 
  facet_grid(~phylo, scales = "free",  space = "free",switch = "x")  +
  ggtitle("B - wzyE")


grid.arrange(A + ggtitle("A"), B + ggtitle("B"), 
             C + ggtitle("C") + theme(legend.position = "right"), layout_matrix = rbind(c(1),
                                                     c(2),
                                                     c(3),
                                                     c(3)))

#grid.arrange(A,B,C + ggtitle("C"),ncol = 1)

## plot presence and absence of two genes -> i.e. which clusters have one gene and not another - on a tree
#gene_names = c("group_5441", "group_2329*************", "group_1059**")
### get the tree object
tree = read.tree("../7_AMR_vir_plasmid/smaller_tree/raxml_tree_mod.nwk")
tree = root(tree,outgroup = "NC_011740")
tree = drop.tip(tree, tip =  "NC_011740")
for (i in c(21,43,49)){
  tree = drop.tip(tree, tip =  as.character(i))
}
p = ggtree(tree)  +
  theme(legend.position="right") + geom_tiplab(align = T)
gene_names = c("sacA","group_1059**")
curr = data.frame(name = rep(o, length(gene_names)),
                  gene = rep("", length.out = length(o)*length(gene_names)),
                  freq = rep(0, length.out = length(o)*length(gene_names)), stringsAsFactors = F)
for (i in 1:length(gene_names)) {
  curr_gene = gene_names[i]
  curr$gene[(length(o)*(i-1) + 1): (length(o)*i)] = rep(curr_gene, length(o))
  curr$freq[(length(o)*(i-1) + 1): (length(o)*i)] = unlist(freqs[which(rownames(freqs) == curr_gene),])
}
curr$name = factor(curr$name, o)
curr = cbind(curr, phylo = graphics$Phylogroup[match(curr$name, graphics$Cluster)])

res = dcast(curr, formula = name~gene, value.var = "freq")
rownames(res) = res[,1]
res = res[,-1]

gheatmap(p = p, data = res, offset = 0.1, width=5, color = "black", 
         high = "#023858", low = "#fff7fb", colnames_angle = 90, colnames = T, colnames_position = "top",
         font.size = 3, colnames_offset_y = 5)
ggplot(curr,aes(x = name, y = gene, fill = freq)) + geom_tile(stat = "identity", color = "black") + 
  theme_classic(base_size = 11) + scale_fill_gradient(low = "white", high = "black") +
  scale_y_discrete(labels = c("effector","hypothetical\nprotein","immunity")) +
  xlab("PopPUNK cluster") + 
  facet_grid(~phylo, scales = "free",  space = "free",switch = "x")


## this is what I would call the core -> these genes are found across all the cluster
## but there are a few clusters that are missing some key genes
df = data.frame(variable = character(0), type = character(0),value = numeric(0),stringsAsFactors = F)
for (size in 1:max(gene_classes$core)) {
  omni = gene_classes[which(gene_classes$total_presence == size),]
  omni = data.frame(table(omni$core, omni$inter, omni$rare), stringsAsFactors = F)
  colnames(omni) = c("core","inter","rare","freq")
  #omni = omni[-which(omni$freq == 0),]
  for (i in 1:3){
    omni[,i] = as.numeric(as.character(omni[,i]))
    omni[,i] = (omni[,i]/size)*omni$freq
  }
  omni_percent = colSums(omni)[-4] / length(which(gene_classes$total_presence == size))
  omni_percent = data.frame(variable = rep(paste(size,"/47",sep=""),3),type = names(omni_percent), value = omni_percent)
  df = rbind(df, omni_percent)
}
p2 =  ggplot(df, aes(x = variable, y = value, fill = type)) + geom_bar(stat = "identity", color = "black") +
  theme_classic(base_size = 12) + scale_fill_manual(values = rev(brewer.pal(n=4,"Blues"))) + theme(legend.position = "bottom")+
  xlab("Number of clusters in which gene is present") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Fraction of genes")
p2
total_presence = cbind(total_presence, lab = paste(total_presence$Var1,"/47",sep=""))
#total_presence = total_presence[-which(total_presence$Var1 == 0),]
p1 =  ggplot(total_presence, aes(x = total_presence$Var1, y = total_presence$Freq)) + geom_bar(stat = "identity") +
  xlab("Number of clusters in which gene is present") + ylab("Number of genes") +
  theme_bw(base_size = 12)+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_x_discrete(labels = total_presence$lab) + ggtitle("D") +
  annotate("text", x = 4, y = 28000, label = format(total_presence$Freq[which(total_presence$Var1 == 1)], big.mark = ",", scientific = F)) + 
  annotate("text", x = 45, y = 3500, label = format(total_presence$Freq[which(total_presence$Var1 == 47)], big.mark = ",", scientific = F))

p1

### Add information of assignment of truncated genes here so I'm aware
## how they distribute
# truncation_assignments = read.table("../4.1_correct_pan_genome/truncated_genes.csv", sep = ",",
#                                     header = T, comment.char = "", stringsAsFactors = F, quote = "")
# truncation_assignments = truncation_assignments$Assignment[match(rownames(gene_classes), truncation_assignments$Gene)]
# truncation_assignments[is.na(truncation_assignments)] = "normal"
# gene_classes = cbind(gene_classes, truncation_assignments)
# gene_classes = gene_classes[-which(is.na(gene_classes$total_presence)),]
# order_truncs = c("normal","common variant (long)", "common variant (short)", "long variant","short variant",
#                  "secondary variant (short)", "secondary variant (long)")
# #gene_classes$truncation_assignments = factor(gene_classes$truncation_assignments, order_truncs)
# gene_classes$total_presence = factor(gene_classes$total_presence, 1:47)
# p1 = ggplot(gene_classes, aes(x = total_presence, fill = truncation_assignments)) + geom_bar(color = "black") +
#   xlab("Number of clusters in which gene is present") + ylab("Number of genes") +
#   theme_bw(base_size = 12)+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("D") +
#   annotate("text", x = 4, y = 28000, label = format(total_presence$Freq[which(total_presence$Var1 == 1)], big.mark = ",", scientific = F)) + 
#   annotate("text", x = 45, y = 3500, label = format(total_presence$Freq[which(total_presence$Var1 == 47)], big.mark = ",", scientific = F)) +
#   scale_fill_manual(values = c(rev(brewer.pal(4, "Greens")[-1]),
#                                rev(brewer.pal(4, "Oranges")[-(1:2)]),
#                                "#FFD92F","#E5C494")) + scale_x_discrete(labels = total_presence$lab)
# legend = as_ggplot(get_legend(p1))
# p1 = p1 + theme(legend.position = "None")
# p1
grid.arrange(intA,wzyE, p1, B, layout_matrix = rbind(c(1,1,1,2,2,2),
                                                     c(1,1,1,2,2,2),
                                                     c(4,4,3,3,3,3),
                                                     c(4,4,3,3,3,3),
                                                     c(4,4,3,3,3,3)))
#write.table(x = gene_classes, file  = "231019_corrected///gene_classes.csv", sep = ",", col.names = T, row.names = T,quote = F)

### ubiquitous genes
ubiq = gene_classes[which(gene_classes$total_presence == 47),]
length(which(ubiq$core == 47))
length(which(ubiq$core >= 40)) / dim(ubiq)[1]

## rare genes
rare = gene_classes[which(gene_classes$total_presence == 1),]
length(which(rare$core == 1))
length(which(rare$inter == 1))
length(which(rare$rare == 1))

length(which(rare$core == 1 & rare$truncation_assignments == 	"secondary variant (short)"))
length(which(rare$core == 1 & rare$truncation_assignments == 	"secondary variant (long)"))

length(which(rare$truncation_assignments != 	"normal"))

