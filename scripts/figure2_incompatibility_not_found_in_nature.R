### Makes figures 2a. First reads in the VCF,  and then makes the figure.
library(foreach)
library(SNPRelate)
library(ape)
library(stringr)
library(phangorn)
library(ggtree)
library(treeio)
library(ggplot2)
genofile <- snpgdsOpen("data/vcf/1012.gds")
dissim = snpgdsDiss(genofile,autosome.only = F)
colnames(dissim$diss) = dissim$sample.id
rownames(dissim$diss) = dissim$sample.id
chinese_annotations1 = read.csv("data/annotations/chinese_annotations.csv",header=T)
chinese_annotations2 = read.delim("data/annotations/chinese_annotations_2.txt", header=T, sep="\t")
chinese_annotations_merged = merge(chinese_annotations1,chinese_annotations2, by.x="BioSample.ID",by.y="BioSample")
chinese_annotations_merged$name = substr(chinese_annotations_merged$Genome.Accession,0,4)
annotation_gal_genes = read.csv("data/annotations/gal_annotations_cbs_complete.csv",stringsAsFactors = F)
colnames(annotation_gal_genes) = c("Standardized.name","GAL1","GAL10","GAL7","GAL2","PGM1")
### Get strain with missing GAL genes.
## NPPO has no copies of gal.
### Renome IDS.
chinese_idx_dissim = na.omit(match(chinese_annotations_merged$Run ,dissim$sample.id))
idx_in_chinese_annotations = sort(na.omit((match(dissim$sample.id,chinese_annotations_merged$Run))))
aa = dissim$sample.id[chinese_idx_dissim] 
aa = chinese_annotations_merged$name[idx_in_chinese_annotations]  
dissim$sample.id[chinese_idx_dissim] = chinese_annotations_merged$name[idx_in_chinese_annotations]  
dissim$sample.id = str_replace(dissim$sample.id,pattern = "SACE_",replacement = "")
rownames(dissim$diss) = dissim$sample.id
colnames(dissim$diss) = dissim$sample.id
# Calculate a neighbour joining tree.
bn = bionj(dissim$diss)
## Midpoint th three
bn_n = midpoint(bn)
aa_2 = as.treedata(bn_n)
actual_groups = rep(NA,length(aa_2@phylo$tip.label))
Div3 = c("BTA","CLG")
Div4 = c("BDN","CER","CET")
Div5 = c("NPPO")
annotations_gal_genes$pgm1[annotations_gal_genes$sample_name == "BAM" | annotations_gal_genes$sample_name == "BAQ"] = "chin"
## Make all the groups
REF = aa_2@phylo$tip.label[aa_2@phylo$tip.label %in% annotations_gal_genes[annotations_gal_genes$gal1 == "wt",]$sample_name]
REF = REF[!(REF %in% c(Div3,Div4))]
actual_groups[aa_2@phylo$tip.label %in% REF] = "Reference"
actual_groups[aa_2@phylo$tip.label %in% annotations_gal_genes$sample_name[annotations_gal_genes$pgm1 == "chin" | annotations_gal_genes$pgm1 == "cbs/wt"]] = "Div2"
actual_groups[aa_2@phylo$tip.label %in% annotations_gal_genes$sample_name[annotations_gal_genes$pgm1 == "cbs"]] = "Div1"
actual_groups[aa_2@phylo$tip.label %in% Div3] = "Div3"
actual_groups[aa_2@phylo$tip.label %in% Div4] = "Div4"
actual_groups[aa_2@phylo$tip.label %in% Div5] = "Div5"
actual_groups[is.na(actual_groups)] = "Reference"
group1 = aa_2@phylo$tip.label[actual_groups == "Reference"] 
group2 = aa_2@phylo$tip.label[actual_groups == "Div1"] 
group3 = aa_2@phylo$tip.label[actual_groups == "Div2"] 
group4 =aa_2@phylo$tip.label[actual_groups == "Div3"]
group5 = aa_2@phylo$tip.label[actual_groups == "Div4"]
group6 = aa_2@phylo$tip.label[actual_groups == "Div5"]
aa_22 = groupOTU(aa_2@phylo,list(group2,group3,group1,group4,group5,group6))
## Figure 2a ##
p = ggtree(aa_22,aes(color=group),layout="fan",open.angle = 180)+scale_color_manual(values = c("#e41a1c","#377eb8","grey85","#f781bf","#984ea3","#ff7f00")) + theme(legend.position="right")
#### ### 
## Figure 2B ##
gal1_10_7_raxml = read.raxml("data/trees/RAxML_bipartitionsBranchLabels.merged_gal1_10_7.raxml.tree")
gal1_10_7_raxml = groupOTU(gal1_10_7_raxml,c("KAFR0F0242","KNAG0K0138"))
ggtree(gal1_10_7_raxml,aes(linetype=group)) +  geom_tiplab(aes(label=label)) + geom_nodelab(aes(label=bootstrap,x=branch),vjust=-0.5) + theme_tree2() -> p2
p2$data[p2$data$node %in% c(7, 8), "x"] = mean(p2$data$x)
svg("figures/2.svg")
plot_grid(p,p2)
dev.off()