
library(foreach)
library(SNPRelate)
library(ape)
library(stringr)
library(phangorn)
library(ggtree)
library(treeio)
library(ggplot2)
Div4 = c("BDN","CER","CET")
tree_384 = read.raxml("data/trees//RAxML_bipartitionsBranchLabels.merged")

#svg("figures/test.svg")
#ggtree(tree_384,layout = "circular") + geom_tiplab()
#dev.off()
strain_annot_gal = read.csv("data/balancing_selection//suptable_annotations2.csv", header=T,stringsAsFactors = F)
with_under = grep("_",tree_384@phylo$tip.label)
#tree_384@phylo$tip.label[tree_384@phylo$tip.label == "CHN"]
tree_384@phylo$tip.label[with_under] = unlist(lapply(str_split(tree_384@phylo$tip.label[with_under],"_"), function(x){x[1]}))
with_perm = grep(".perm",tree_384@phylo$tip.label)
tree_384@phylo$tip.label[with_perm] = unlist(lapply(str_split(tree_384@phylo$tip.label[with_perm],".perm"), function(x){x[1]}))
g1 = strain_annot_gal$Sandardized.name[strain_annot_gal$GAL1 == "REF"]
g2 = strain_annot_gal$Standardized.name[strain_annot_gal$GAL1 == "DIV"]
g2 = c(g2,"DIV","CHIN")
g1 = c(g1, "YCL")
g1 = c(g1,"YBR018C")
g3 = tree_384@phylo$tip.label[!(tree_384@phylo$tip.label %in% c(g1 , g2))]
aa = groupOTU(tree_384,list( g1,g2,g3))
ggtree(aa,layout = "circular") + geom_tiplab()
g2[!(g2 %in% tree_384@phylo$tip.label)]
svg("figures/S25.svg",width=8,height=8)
ggtree(aa, aes(color=group),size=.4,layout="fan") + scale_color_manual(values = c("grey50","#e41a1c","black")) + geom_treescale()
dev.off()
