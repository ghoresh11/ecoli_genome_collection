library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(ggpubr)
library(ggtree)
library(reshape2)
library(ape)
library(dplyr)

######## General things ###########
setwd("/Users/gh11/poppunk_pangenome/7_AMR_vir_plasmid//")

orig_md = read.table("/Users/gh11/e_colis/FILTERED_MD_FINAL_ALL.tab", sep = "\t",
                     header = T, comment.char = "", quote = "", stringsAsFactors = F)
md_names = as.character(sapply(sapply(sapply(sapply(orig_md$New_annot_loc, strsplit, split = "/", fixed = T),tail, n = 1), 
                                      strsplit, split = ".",fixed = T), head, n=1))
rownames(orig_md) = md_names

cluster_graphics = read.table("/Users/gh11/Submissions/my_thesis/Chapter4//figures/cluster_graphics.csv",
                              sep = ",", comment.char = "", header = T, stringsAsFactors = F)

### plotting heatmaps of genes along the tree
tree = read.tree("smaller_tree/raxml_tree_mod.nwk")
tree = root(tree,outgroup = "NC_011740")
#plot(tree)
tree = drop.tip(tree, tip =  "NC_011740")
#plot(tree)

write.tree(tree, file = "/Users/gh11/Submissions/bioresource/data/tree_50.nwk")



cluster_order = # (taken from tree)
  rev(c(8,36,43,24,42,32,26,49,9,34,48,5,21,15,6,14,39,22,40,45,30,28,12,35,23,1,16,10,51,19,7,37,46,29,
        4,18,25,13,3,41,17,11,27,20,2,47,38,33,31,44))

############## FUNCTIONS #############

get_gene_for_one_cluster <- function(cluster, md, df){
  df = df[which(md$Poppunk_cluster == cluster),]
  df[df=="1*"] = 0 ## here: decide what to do with proteins which appear to be truncated (maybe want to look at them)
  df = data.frame(lapply(df, as.numeric))
  freq = colSums(df, na.rm = T)/dim(df)[1]
  # df[df=="1"] = 0
  # df[is.na(df)] = 1
  # freq_trunc = colSums(df, na.rm = T)/dim(df)[1]
  res = data.frame(gene = rep(names(freq), 1),
                   cluster = rep(cluster,length(freq)),
                   #  type = c(rep("full", length(freq)), rep("trunc",length(freq))),
                   freq = c(freq), stringsAsFactors = F)
  return(res)
}

plot_on_tree <- function(filename, desc_file, tree, signif, display_names = T){
  df = read.table(filename, sep = ",", header = T, comment.char = "",
                  quote = "", stringsAsFactors = F, row.names = 1)
  md = orig_md[match(rownames(df), rownames(orig_md)),]
  dfs = lapply(X = unique(md$Poppunk_cluster), FUN = get_gene_for_one_cluster, md = md, df = df)
  res = do.call(rbind, dfs)
  res = dcast(res, cluster ~ gene, value.var = "freq")
  rownames(res) = res[,1]
  res = res[,-1]
  res = res[match(tree$tip.label, rownames(res)),]
  
  ## change the order of the columns
  res = res[,order(colSums(res), decreasing = T)]
  res = res[,-which(sapply(FUN = max, X = res)<0.1)] ## to remove very rare genes
  
  ## create an output file for the cluster summaries
  df_out = data.frame(cluster = rownames(res),
                      label = rep("-", dim(res)[1]), stringsAsFactors = F)
  for (i in 1:dim(res)[1]) {
    props = sort(round(res[i,][which(res[i,] >0)], digits = 2), decreasing = T)
    df_out$label[i] = paste(names(props), props, sep = ":", collapse = "/")
  }
  write.table(df_out, file = paste(filename, "_per_cluster.csv", sep = ""), col.names = T, row.names = F, quote = F, sep = ",")
  
  descs = read.table(desc_file, sep = ",", header = T, stringsAsFactors = F, quote = "")
  new_labs = descs$gene[match(colnames(res), descs$identifier)]
  labs = c()
  for (l in new_labs) {
    while (l %in% labs) {
      l = paste(l, "^", sep = "")
    }
    labs = c(labs,l)
  }
  colnames(res) = labs
  col.order <- rev(hclust(dist(t(res)))$order)
  res = res[,col.order]
  
  tip_labels = tree$tip.label
  for (i in 1:length(tip_labels)) {
    if (tip_labels[i] %in% signif){
      tip_labels[i] = paste(tip_labels[i], "*", sep = "")
    }
  }
  tree$tip.label = tip_labels
  rownames(res) = tip_labels
  p = ggtree(tree)  +
    theme(legend.position="right") + geom_tiplab(align = T)
  # tree$tip.label = factor(tree$tip.label, tree$tip.label)
  p2 = gheatmap(p = p, data = res, offset = 0.1, width=5, color = "black", 
                high = "#023858", low = "#fff7fb", colnames = display_names, 
                font.size = 4, colnames_angle = 90)
  
  
  
  return(p2)
}


### AMR genes
amr_results = read.table("/Users/gh11/poppunk_pangenome/7_AMR_vir_plasmid/results/resfinder.csv", sep = ",",
                         comment.char = "", stringsAsFactors = F, header = T, row.names = 1)
amr_results[amr_results == "1*"] = 0
amr_results_2 = mutate_all(amr_results, function(x) as.numeric(as.character(x)))
rownames(amr_results_2) = rownames(amr_results)
amr_results = amr_results_2
## all amr
all_amr = read.table("/Users/gh11/poppunk_pangenome/7_AMR_vir_plasmid/DBs/break_names/resfinder_genes_fixed.csv", sep = ",",
                     header = T, stringsAsFactors = F, comment.char = "", quote = "")
all_amr = all_amr[match(colnames(amr_results),all_amr$identifier),]
antibiotic_classes = read.table("antibiotic_classes.csv", sep = ",", quote = "", comment.char = "", header = T, stringsAsFactors = F)

ggplot(antibiotic_classes, aes(x = ))

# # 1.Count how many antibiotic classes each strain is resistant to
# antibiotic_classes = data.frame(strain = rownames(amr_results),
#                                 count = rep(0, dim(amr_results)[1]),
#                                 names =  rep("-", dim(amr_results)[1]), stringsAsFactors = F)
# 
# for (i in 1:dim(antibiotic_classes)[1]) {
#   print(i)
#   results =  unique(all_amr$antimicrobial.category[which(amr_results[i,]>0)])
#   if (length(results) == 0) {next}
#   antibiotic_classes$names[i] = paste(sort(results), collapse = ";")
#   antibiotic_classes$count[i] = length(results)
# }
# antibiotic_classes$PopPunk_cluster = orig_md$Poppunk_cluster[match(antibiotic_classes$strain, rownames(orig_md))]
# antibiotic_classes$Phylogroup = cluster_graphics$Phylogroup[match(antibiotic_classes$PopPunk_cluster, cluster_graphics$Cluster)]


median_classes = aggregate(x = antibiotic_classes$count, by = list(antibiotic_classes$PopPunk_cluster), median)
cluster_graphics$MDR = rep("No",dim(cluster_graphics)[1])
cluster_graphics$MDR[which(cluster_graphics$Cluster %in% median_classes$Group.1[median_classes$x >=3])] = "Yes"
# antibiotic_classes$MDR = cluster_graphics$MDR[match(antibiotic_classes$PopPunk_cluster, cluster_graphics$Cluster)]
# write.table(antibiotic_classes, file = "antibiotic_classes.csv", sep = ",", col.names = T, row.names = F, quote = F)

antibiotic_classes$PopPunk_cluster = factor(antibiotic_classes$PopPunk_cluster, cluster_order)
antibiotic_classes$fill = rep("No", dim(antibiotic_classes)[1])
antibiotic_classes$fill[antibiotic_classes$PopPunk_cluster %in% cluster_graphics$Cluster[cluster_graphics$MDR == "Yes"]] = "Yes"

A = ggplot(antibiotic_classes, aes(x = PopPunk_cluster, y = count, fill = fill)) + geom_boxplot(outlier.size = 0.5) +
  geom_hline(yintercept = 3, color = "red", lwd = 0.8) +
  theme_classic(base_size = 12) + xlab("PopPUNK Cluster") + ylab(ylab) +
  scale_fill_manual(values = c("#e0e0e0","#fddbc7"), name = "MDR")+ 
  scale_y_continuous(breaks = seq(from = 0, to = 8, by = 2)) + ylab("Antimicrobial classes") + xlab("Lineage") +
  coord_flip() 


# signif = cluster_graphics$Cluster[cluster_graphics$MDR == "Yes"]
# B = plot_on_tree("results/resfinder.csv","DBs/break_names/resfinder_genes_fixed.csv", tree, signif )
# grid.arrange(A + theme(legend.position = "bottom") + ylab("Antimicrobial Classes") + ggtitle("A"),
#              B + ggtitle("B"),  layout_matrix = rbind(c(1,2,2)))


A_amr = A


### virulence genes
vir_results = read.table("/Users/gh11/poppunk_pangenome/7_AMR_vir_plasmid/results/virulence.csv", sep = ",",
                         comment.char = "", stringsAsFactors = F, header = T, row.names = 1)
vir_results[vir_results == "1*"] = 0
vir_results_2 = mutate_all(vir_results, function(x) as.numeric(as.character(x)))
rownames(vir_results_2) = rownames(vir_results)
vir_results = vir_results_2


pathotype_genes = read.table("DBs/break_names/pathotypes.csv", sep = ",", header = T, stringsAsFactors = F)
all_virulence_genes = read.table("DBs/break_names/virulence.csv", sep = ",", header = T, stringsAsFactors = F)
identifiers = all_virulence_genes$identifier[which(all_virulence_genes$gene %in% pathotype_genes$Gene)]
filtered_vir_results =  vir_results[,colnames(vir_results) %in% identifiers]


## predict the pathotype of isolate
pathotypes = data.frame(strain = rownames(filtered_vir_results),
                        pathotype = orig_md$Pathotype[match(rownames(filtered_vir_results), rownames(orig_md))],
                        PopPUNK_cluster = orig_md$Poppunk_cluster[match(rownames(filtered_vir_results), rownames(orig_md))], 
                        stringsAsFactors = F)
pathotypes = cbind(pathotypes, filtered_vir_results)

## change the pathotype column to include only clear ones
pathotypes$pathotype[which(!pathotypes$pathotype %in% c("nd","expec","epec","stec","etec","eaec"))] = "nd"
pathotypes$pathotype = toupper(pathotypes$pathotype)

## this works for ExPEC
for (i in 1:length(pathotypes$strain)) {
  strain = pathotypes$strain[i]
  index = which(rownames(orig_md) == strain)
  if (orig_md$Isolation[index] %in% c("blood","urine")) {
    if (!pathotypes$pathotype[i] %in% c("EXPEC","ND")) { ## sanity check
      print(pathotypes$pathotype[i])
    }
    pathotypes$pathotype[i] = "ExPEC"
  }
}

for (gene in pathotype_genes$Gene) {
  print(gene)
  p = pathotype_genes$Pathotype[pathotype_genes$Gene == gene]
  curr_idents = all_virulence_genes$identifier[which(all_virulence_genes$gene == gene)]
  for (ident in curr_idents) {
    for (i in 1:dim(pathotypes)[1]) {
      if (pathotypes$pathotype[i] == "ExPEC") {next}
      if (pathotypes[i,which(colnames(pathotypes) == ident)] == 1) {
        if (pathotypes$pathotype[i] == "ND") {
          pathotypes$pathotype[i] = p
        } else if (!grepl(x = pathotypes$pathotype[i], pattern = p, fixed = T)) {
          pathotypes$pathotype[i] = paste(pathotypes$pathotype[i], "/",p)
        }
      }
    }
  }
}
## correct to include pathotypes that are combinations of other pathotypes
for (i in 1:dim(pathotypes)[1]) {
  if (grepl(x = pathotypes$pathotype[i], pattern = "EPEC", fixed = T) && 
      grepl(x = pathotypes$pathotype[i], pattern = "STEC", fixed = T) ) {
    pathotypes$pathotype[i] = "EHEC"
  }
  else if (grepl(x = pathotypes$pathotype[i], pattern = "STEC", fixed = T) && 
           grepl(x = pathotypes$pathotype[i], pattern = "EAEC", fixed = T) ) {
    pathotypes$pathotype[i] = "EAEC+STEC"
  }
  elements = unlist(strsplit(x = pathotypes$pathotype[i], split = " / ", fixed = T))
  if (("EPEC" %in% elements) && ("aEPEC" %in% elements) && length(elements) == 2) {
    pathotypes$pathotype[i] = "EPEC"
  }
  else if ( ("ETEC" %in% elements) &&
           (("EPEC" %in% elements) || ("aEPEC" %in% elements))) {
    pathotypes$pathotype[i] = "EPEC/ETEC"
  }
}

pathotypes$pathotype[pathotypes$pathotype %in% c("EPEC","aEPEC")] = "aEPEC/EPEC"
write.table(x = pathotypes, file = "pathotype_per_isolate.csv", sep = ",",
            col.names = T, row.names = F, quote = F)



counts_per_cluster = table(pathotypes$PopPUNK_cluster,pathotypes$pathotype)
pathotype_per_cluster = data.frame(PopPUNK_cluster = rownames(counts_per_cluster),
                                   Pathotype_details = rep("", dim(counts_per_cluster)[1]),
                                   Pathotype_label =  rep("", dim(counts_per_cluster)[1]),
                                   stringsAsFactors = F)

for (i in 1:dim(counts_per_cluster)[1]) {
  counts_per_cluster[i,] = round(counts_per_cluster[i,]/sum(counts_per_cluster[i,]), digits = 2)
  chosen = which(counts_per_cluster[i,] > 0)
  pathotype_per_cluster$Pathotype_details[i] = paste(colnames(counts_per_cluster)[chosen], "(",
  counts_per_cluster[i,chosen],")", sep = "", collapse = ";")
  
  chosen = which(counts_per_cluster[i,] > 0.5)
  if (length(chosen) == 0) {
    pathotype_per_cluster$Pathotype_label[i] = "No prevalent pathotype"
    next
  }
  pathotype_per_cluster$Pathotype_label[i] = colnames(counts_per_cluster)[chosen]
}
write.table(x = pathotype_per_cluster, file = "pathotypes_per_cluster.csv", sep = ",",
            col.names = T, row.names = F, quote = F)


### create the boxplot
vir_per_isolate = data.frame(strain = rownames(vir_results),num_genes = rowSums(vir_results), stringsAsFactors = F)
vir_per_isolate$PopPUNK_cluster = orig_md$Poppunk_cluster[match(vir_per_isolate$strain, rownames(orig_md))]
vir_per_isolate$pathotype = pathotype_per_cluster$Pathotype_label[match(vir_per_isolate$PopPUNK_cluster, pathotype_per_cluster$PopPUNK_cluster)]
vir_per_isolate$pathotype[vir_per_isolate$PopPUNK_cluster %in% c(30, 45, 23)] = "EIEC/Shigella"


vir_per_isolate$PopPUNK_cluster = factor(vir_per_isolate$PopPUNK_cluster, cluster_order)

vir_per_isolate$pathotype = factor(vir_per_isolate$pathotype, c("ExPEC", "EHEC", "STEC","aEPEC/EPEC","EAEC","EIEC/Shigella","No prevalent pathotype","ND"))



A = ggplot(vir_per_isolate, aes(x = PopPUNK_cluster, y = num_genes, fill = pathotype)) + geom_boxplot(outlier.size = 0.5) +
  coord_flip() + theme_classic(base_size = 12) + ylab("Virulence factors") +xlab("Lineage") +
  scale_y_continuous(breaks = seq(from= 0, to = 25, by = 4))+
  scale_fill_manual(values = c( brewer.pal(n=7, "Set2"), NA), name = "Pathotype") 

A_vir = A


vir_per_isolate$phylogroup = cluster_graphics$Phylogroup[match(vir_per_isolate$PopPUNK_cluster, cluster_graphics$Cluster)]
ggplot(vir_per_isolate, aes( y = num_genes, x = pathotype)) + geom_boxplot(outlier.size = 0.5) +
  coord_flip() + theme_classic(base_size = 14) + ylab("Virulence genes") +xlab("Pathotype") +
  scale_y_continuous(breaks = seq(from= 0, to = 25, by = 4)) 

# B = plot_on_tree("results/virulence.csv","DBs/break_names/virulence.csv", tree, c(),display_names = F)
# B
# 
# grid.arrange(A, B, layout_matrix = rbind(c(1,1,2,2,2)))


### relationship between resistance and virulence

amr =  aggregate(x = antibiotic_classes$count, by = list(antibiotic_classes$PopPunk_cluster), FUN = median)
vir = aggregate(x = vir_per_isolate$num_genes, by = list(vir_per_isolate$PopPUNK_cluster), FUN = median)
amr = amr[match(vir$Group.1, amr$Group.1),]
vir_and_amr = data.frame(PopPUNK_cluster = amr$Group.1,
                         amr = amr$x,
                         vir = vir$x,
                         pathotype = pathotype_per_cluster$Pathotype_label[match(as.character(amr$Group.1), pathotype_per_cluster$PopPUNK_cluster)],
                         stringsAsFactors = F)

vir_and_amr$pathotype[vir_and_amr$PopPUNK_cluster %in% c(30, 45, 23)] = "EIEC/Shigella"
vir_and_amr$pathotype = factor(vir_and_amr$pathotype, c("ExPEC", "EHEC", "STEC","aEPEC/EPEC","EAEC","EIEC/Shigella","No prevalent pathotype","ND"))
pd <- position_dodge(0.3)
C = ggplot(vir_and_amr, aes(x = amr, y = vir, fill = pathotype, label = PopPUNK_cluster)) + geom_point(pch = 21, size = 8) +
  geom_text()+
  theme_bw(base_size = 12) +
  scale_fill_manual(values = c( brewer.pal(n=7, "Set2"), NA), name = "Pathotype") +
  scale_x_continuous(breaks = seq(from=0, to = 9, by = 1)) +
  scale_y_continuous(breaks = seq(from = 0, to = max(vir$x), by = 2)) + xlab("Median antimicrobial resistance classes") +
  ylab("Median virulence factors per isolate") +
  geom_vline(xintercept = 3, lty = 2, colour = "red")

grid.arrange(A_amr + theme(legend.position = "None") + ggtitle("A"), 
             A_vir  + theme(legend.position = "None") + ggtitle("B"), 
             C  + theme(legend.position = "None") + ggtitle("C"), layout_matrix = rbind(c(1,2),
                                                    c(1,2),
                                                    c(1,2),
                                                    c(3,3),
                                                    c(3,3)))
amr_legend = as_ggplot(get_legend(A_amr))
vir1_legend = as_ggplot(get_legend(A_vir))
vir2_legend = as_ggplot(get_legend(C))

calc_props_one_cluster <- function(curr_cluster) {
  curr = antibiotic_classes[which(antibiotic_classes$PopPunk_cluster == curr_cluster),]
  curr$path = vir_per_isolate$num_genes[vir_per_isolate$PopPUNK_cluster == curr_cluster]
  res = data.frame(table(curr$count, curr$path),stringsAsFactors = F)
  res$Freq = res$Freq/dim(curr)[1]
  res$name = paste(res$Var1, res$Var2, sep ="-")
  return(res)
}

all_res = calc_props_one_cluster(pathotype_per_cluster$PopPUNK_cluster[1])
for (i in 2:length(pathotype_per_cluster$PopPUNK_cluster)){
  res = calc_props_one_cluster(pathotype_per_cluster$PopPUNK_cluster[i])
  in_both = res$Freq[match(all_res$name,res$name)]
  in_both[which(is.na(in_both))] = 0
  all_res$Freq = all_res$Freq + in_both
  all_res = rbind(all_res, res[which(!res$name %in% all_res$name),])
}

all_res$Var1 = as.numeric(all_res$Var1)
all_res$Var2 = as.numeric(all_res$Var2)
ggplot(all_res, aes(x = Var1, y = Var2, size = Freq)) + geom_point() +
  scale_size_continuous(breaks = seq(from = 0.1, to = 2, by = 0.1), range = c(0,9))


