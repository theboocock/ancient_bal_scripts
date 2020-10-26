# Plots figure S3.
library(tidyverse)
igff= rtracklayer::import.gff3("data/saccharomyces_cerevisiae.gff")
gff = gff[which(gff$orf_classification == "Verified" & gff$type == "gene"),]
distances_genes = read.table("data/popgen_in/identity_genes_cbs.txt", header=F)
distances = read.table("data/popgen_in/distances.txt", header=F)
gal_genes = c("YBR018C","YLR081W","YBR019C","YBR020W")
without_pgm1 = distances$V1[distances$V1 <= 0.46]
pgm1_df = data.frame(distances=without_pgm1)
p1 = pgm1_df %>% ggplot(aes(x=(1-distances)*100)) + theme_bw() + geom_histogram() + coord_cartesian(xlim=c(.40*100,1*100)) + 
  geom_vline(xintercept = (1-0.5180023228803716)*100,color="red") + xlab("Percent identity") + ylab("Count") + 
  theme(axis.title = element_text(size=22), axis.text = element_text(size=22),title =  element_text(size=22))+ ggtitle("Promoters") 
ab_line_genes =100 *mean(distances_genes$V2[distances_genes$V1 %in% gal_genes])
genes_df = data.frame(distances=distances_genes$V2[!(distances_genes$V1 %in% gal_genes)])     
genes_df$distances = genes_df$distances * 100
p2 = genes_df %>% ggplot(aes(x=distances)) +  theme_bw() + geom_histogram() + coord_cartesian(xlim=c(.7*100,1*100)) +
  geom_vline(xintercept = ab_line_genes,color="red") + xlab("Percent identity")   + ylab("Count") + 
  theme(axis.title = element_text(size=22), axis.text = element_text(size=22),title =  element_text(size=22))+ ggtitle("Genes") 
svg("figures/S3.svg")
cowplot::plot_grid(p2,p1,labels=c("a)","b)"))
dev.off()

