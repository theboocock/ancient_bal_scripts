### Background distribution identity.
dn_ds = read.table("data/popgen_in/dn_ds/dn_ds.txt", header=F)
#dn_ds = read.table("popgen_phylo_notebooks/tmp/genes_cbs_new/dn_ds.txt", header=F)
#library(tidyverse)
#distances_genes = read.table("popgen_phylo_notebooks/data/identity_genes_cbs.txt", header=F)
#distances = read.table("popgen_phylo_notebooks/outputs/intergenic_alignments/distances.txt", header=F)
gal_genes = c("YBR018C","YLR081W","YBR019C","YBR020W")
#without_pgm1 = distances$V1[distances$V1 <= 0.46]
#pgm1_df = data.frame(distances=without_pgm1)
#p1 = pgm1_df %>% ggplot(aes(x=(1-distances)*100)) + theme_bw() + geom_histogram() + coord_cartesian(xlim=c(.40*100,1*100)) + 
 # geom_vline(xintercept = (1-0.5180023228803716)*100,color="red") + xlab("Percent identity") + ylab("Count") + 
 # theme(axis.title = element_text(size=22), axis.text = element_text(size=22),title =  element_text(size=22))+ ggtitle("Promoters") 
#ab_line_genes =100 *mean(distances_genes$V2[distances_genes$V1 %in% gal_genes])
#genes_df = data.frame(distances=distances_genes$V2[!(distances_genes$V1 %in% gal_genes)])     
#genes_df$distances = genes_df$distances * 100
##p2 = genes_df %>% ggplot(aes(x=distances)) +  theme_bw() + geom_histogram() + coord_cartesian(xlim=c(.7*100,1*100)) +
#  geom_vline(xintercept = ab_line_genes,color="red") + xlab("Percent identity")   + ylab("Count") + 
#  theme(axis.title = element_text(size=22), axis.text = element_text(size=22),title =  element_text(size=22))+ ggtitle("Genes") 
#svg("figures/S2.svg",width=8,height=6)
#cowplot::plot_grid(p2,p1,labels=c("a)","b)"))
#dev.off()
umean_gal_genes = 2.428821
#(dn_ds %>% filter(V1 %in% gal_genes) %>% summarise(mean(V4)))
#svg("figures/S8.svg",width=10, height=6)
svg("figures/S27.svg")
dn_ds %>% filter(!(V1 %in% gal_genes)) %>% ggplot(aes(x=V4)) +  theme_bw() + geom_histogram() + xlab("Synonymous substitions per site")   + ylab("Count") + 
  theme(axis.title = element_text(size=22), axis.text = element_text(size=22),title =  element_text(size=22)) + geom_vline(xintercept = 2.40,color="red",size=1) 
dev.off()
genes_locations = merge(as.data.frame(gff),dn_ds,by.x="Name",by.y=1)
ds_table = genes_locations[,c("Name","gene","V4")]
write.csv(ds_table,file="figures/ds_table.csv")
