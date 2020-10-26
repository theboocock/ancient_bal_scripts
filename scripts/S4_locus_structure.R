library(gggenes)
gal4binding_sites = c()
gene_locations = c()
by_genes = read.csv("data/ggenes//by_genes_gal2.csv", header=T)
by_genes$Gene = factor(by_genes$Gene)
by_genes$strand_t = by_genes$strand  == "forward"
myColors = RColorBrewer::brewer.pal(11,"Set3")
names(myColors) = levels(by_genes$Gene)
by_bs = by_genes[by_genes$type == "BS" & by_genes$locus == "GAL1",]
#by_genes %>% filter(type == "gene" & locus == "GAL1" ) %>% ggplot(aes(xmin=start,xmax=end,y=Strain,fill=Gene,label=Gene,forward=strand_t)) +   geom_gene_arrow(arrow_body_height  = unit(10, "mm"), arrowhead_height = unit(10,"mm"), arrowhead_width = unit(5, "mm"))  + geom_gene_arrow(data=by_bs, aes(xmin=start,xmax=end,y=Strain,fill=Gene,label=Gene,forward=strand_t),arrow_body_height  = unit(10, "mm"),arrowhead_height = unit(3, "mm"), arrowhead_width = unit(0, "mm")) +   geom_gene_label(height = unit(10,"mm"),grow=T)  + theme_genes() +  scale_fill_manual(name="Gene",values=myColors) + ggtitle("PGM1") +  theme(plot.title = element_text(hjust = 0.5))
p1= by_genes %>% filter(type == "gene" & locus == "GAL1" ) %>% ggplot(aes(xmin=start,xmax=end,y=Strain,fill=Gene,label=Gene,forward=strand_t)) +   geom_gene_arrow(arrow_body_height  = unit(10, "mm"), arrowhead_height = unit(10,"mm"), arrowhead_width = unit(5, "mm")) + geom_gene_arrow(data=by_bs, aes(xmin=start,xmax=end,y=Strain,fill=Gene,label=Gene,forward=strand_t),arrow_body_height  = unit(10, "mm"),arrowhead_height = unit(3, "mm"), arrowhead_width = unit(0, "mm")) +   geom_gene_label(height = unit(10,"mm"),grow=T) + theme_genes() +scale_fill_manual(name="Gene",values=myColors) + ggtitle("GAL1/10/7") +  theme(plot.title = element_text(hjust = 0.5))
#    by_g = by_genes[by_genes$type == "gene" & by_genes$locus == "PGM1",]
by_bs = by_genes[by_genes$type == "BS" & by_genes$locus == "PGM1",]
p2 = by_genes %>% filter(type == "gene" & locus == "PGM1" ) %>% ggplot(aes(xmin=start,xmax=end,y=Strain,fill=Gene,label=Gene,forward=strand_t)) +   geom_gene_arrow(arrow_body_height  = unit(10, "mm"), arrowhead_height = unit(10,"mm"), arrowhead_width = unit(5, "mm")) + geom_gene_arrow(data=by_bs, aes(xmin=start,xmax=end,y=Strain,fill=Gene,label=Gene,forward=strand_t),arrow_body_height  = unit(10, "mm"),arrowhead_height = unit(3, "mm"), arrowhead_width = unit(0, "mm")) +   geom_gene_label(height = unit(10,"mm"),grow=T) + theme_genes() +  scale_fill_manual(name="Gene",values=myColors) + ggtitle("PGM1") +  theme(plot.title = element_text(hjust = 0.5))
by_g = by_genes[by_genes$type == "gene" & by_genes$locus == "GAL2",]
by_bs = by_genes[by_genes$type == "BS" & by_genes$locus == "GAL2",]
p3 = ggplot(by_g,aes(xmin=start,xmax=end,y=Strain,fill=Gene,label=Gene,forward=strand_t)) +   geom_gene_arrow(arrow_body_height  = unit(10, "mm"), arrowhead_height = unit(10,"mm"), arrowhead_width = unit(5, "mm")) + geom_gene_arrow(data=by_bs, aes(xmin=start,xmax=end,y=Strain,fill=Gene,label=Gene,forward=strand_t),arrow_body_height  = unit(10, "mm"),arrowhead_height = unit(3, "mm"), arrowhead_width = unit(0, "mm")) +   geom_gene_label(height = unit(10,"mm"),grow=T) + theme_genes() +scale_fill_manual(name="Gene",values=myColors)  + ggtitle("GAL2") +  theme(plot.title = element_text(hjust = 0.5))
svg("figures/S4.svg")
gridExtra::grid.arrange(p1,p2,p3,nrow=3)
dev.off()
## Figure not included in paper shows where primers were designed. 
ar_strains = read.csv("data/ggenes/ar_strains.csv",header=T)
ar_strains$strand_t = ar_strains$strand == "forward"   
by_g = ar_strains[(ar_strains$type == "Gene"  | ar_strains$type == "Repair Template Primers")& ar_strains$locus == "GAL1",]
by_bs = ar_strains[ar_strains$type == "GAL4BS" & ar_strains$locus == "GAL1",]
myColors = RColorBrewer::brewer.pal(14,"Set3")
names(myColors) = levels(by_genes$Gene)
p4 = ggplot(by_g,aes(xmin=start,xmax=end,y=Strain,fill=type,label=Gene,forward=strand_t)) +   geom_gene_arrow(arrow_body_height  = unit(10, "mm"), arrowhead_height = unit(10,"mm"), arrowhead_width = unit(5, "mm")) + geom_gene_arrow(data=by_bs, aes(xmin=start,xmax=end,y=Strain,fill=type,label=Gene,forward=strand_t),arrow_body_height  = unit(10, "mm"),arrowhead_height = unit(3, "mm"), arrowhead_width = unit(0, "mm")) +   geom_gene_label(height = unit(10,"mm"),grow=T) + theme_genes() +scale_fill_brewer(palette = "Set3")+ ggtitle("GAL1/10/7") +  theme(plot.title = element_text(hjust = 0.5))
by_g = ar_strains[(ar_strains$type == "Gene"  | ar_strains$type == "Repair Template Primers")& ar_strains$locus == "PGM1",]
by_bs = ar_strains[ar_strains$type == "GAL4BS" & ar_strains$locus == "PGM1",]
p5 = ggplot(by_g,aes(xmin=start,xmax=end,y=Strain,fill=type,label=Gene,forward=strand_t)) +   geom_gene_arrow(arrow_body_height  = unit(10, "mm"), arrowhead_height = unit(10,"mm"), arrowhead_width = unit(5, "mm")) + geom_gene_arrow(data=by_bs, aes(xmin=start,xmax=end,y=Strain,fill=type,label=Gene,forward=strand_t),arrow_body_height  = unit(10, "mm"),arrowhead_height = unit(3, "mm"), arrowhead_width = unit(0, "mm")) +   geom_gene_label(height = unit(10,"mm"),grow=T) + theme_genes() +scale_fill_brewer(palette = "Set3")+ ggtitle("PGM1") +  theme(plot.title = element_text(hjust = 0.5))
by_g = ar_strains[(ar_strains$type == "Gene"  | ar_strains$type == "Repair Template Primers")& ar_strains$locus == "GAL2",]
by_bs = ar_strains[ar_strains$type == "GAL4BS" & ar_strains$locus == "GAL2",]
p6 = ggplot(by_g,aes(xmin=start,xmax=end,y=Strain,fill=type,label=Gene,forward=strand_t)) +   geom_gene_arrow(arrow_body_height  = unit(10, "mm"), arrowhead_height = unit(10,"mm"), arrowhead_width = unit(5, "mm")) + geom_gene_arrow(data=by_bs, aes(xmin=start,xmax=end,y=Strain,fill=type,label=Gene,forward=strand_t),arrow_body_height  = unit(10, "mm"),arrowhead_height = unit(3, "mm"), arrowhead_width = unit(0, "mm")) +   geom_gene_label(height = unit(10,"mm"),grow=T) + theme_genes() +scale_fill_brewer(palette = "Set3")+ ggtitle("GAL2") +  theme(plot.title = element_text(hjust = 0.5))
gridExtra::grid.arrange(p4,p5,p6)
