in_all = list.files("data/balancing_selection/gall_dn_ds_flanks_all/", pattern="*dn_ds",full.names = T)
in_cbs = list.files("data/balancing_selection/gall_dn_ds_flanks_all_cbs_ref/", pattern="*dn_ds", full.names = T)
strain_annot_gal = read.csv("data/balancing_selection//suptable_annotations2.csv", header=T)
gal_genes = c("YBR018C","YBR019C","YBR020W","YLR081W")
gff= import.gff3("data/saccharomyces_cerevisiae.gff")
gff = gff[which(gff$orf_classification == "Verified" & gff$type =="gene"),]
### TODO: DELETE CRL__
loci = 1
regions_overall = list()
regions_cbs = list()
for (file in in_all){
  print(file)
  a = read.table(file, header=F,stringsAsFactors = F)
  
  #print(head(a))
  sample_names = str_split(str_replace_all(basename(file),"_[0-9]",""),"\\.")[[1]][1]
  # sample_names_c = c(sample_names_c, sample_names)
  overall = a %>% filter(V2 == "WINDOW") %>% filter(!(V1 %in% gal_genes))
  annot_gal = strain_annot_gal %>% filter(strain_annot_gal$Standardized.name ==sample_names)
  if(nrow(overall) == 0){
    print(sample_names)
    next
  }
  overall$sample_name = sample_names
  overall$order = 1:nrow(overall)
  
  overall_df = merge(overall, annot_gal, by.x="sample_name",by.y=1)
  overall_df = overall_df[order(overall_df$order),]
  overall_df = merge(overall_df,as.data.frame(gff),by.x="V1",by.y="ID")
  create_indexes = lapply(split(as.data.frame(overall_df), f=factor(overall_df$V1)),function(x){
    #print(x)
    #print(seq(from=0,to=length(x$start)*30 - 1,by=30));
    #print(x$start + seq(from=0,to=length(x$start)*30-1,by=30));
    x$pos = x$start + seq(from=0,to=length(x$start)*30-1,by=30);
    ##print(x$pos)
    #print(x$pos);  print("HERE");print(x);return(x)}
    return(x)
  })
  
  gff_all_random = do.call("rbind",create_indexes)
  #gff_all_random$pos2 = gff_all_random_idx$pos - position
  # gff_all_random_idx$index = loci
  
  gff_all = (gff_all_random[,c(1:14)])
  gff_all$pos = gff_all_random$pos
  ### other file
  in_cbs_f = read.table(in_cbs[loci], stringsAsFactors = F, header=F)
  
  #print(head(a))
  sample_names = str_split(str_replace_all(basename(in_cbs[loci]),"_[0-9]",""),"\\.")[[1]][1]
  # sample_names_c = c(sample_names_c, sample_names)
  overall = in_cbs_f %>% filter(V2 == "WINDOW") %>% filter(!(V1 %in% gal_genes))
  annot_gal = strain_annot_gal %>% filter(strain_annot_gal$Standardized.name ==sample_names)
  if(nrow(overall) == 0){
    print(sample_names)
    next
  }
  overall$sample_name = sample_names
  overall$order = 1:nrow(overall)
  
  overall_df = merge(overall, annot_gal, by.x="sample_name",by.y=1)
  overall_df = overall_df[order(overall_df$order),]
  overall_df = merge(overall_df,as.data.frame(gff),by.x="V1",by.y="ID")
  create_indexes = lapply(split(as.data.frame(overall_df), f=factor(overall_df$V1)),function(x){
    #print(x)
    #print(seq(from=0,to=length(x$start)*30 - 1,by=30));
    #print(x$start + seq(from=0,to=length(x$start)*30-1,by=30));
    x$pos = x$start + seq(from=0,to=length(x$start)*30-1,by=30);
    ##print(x$pos)
    #print(x$pos);  print("HERE");print(x);return(x)}
    return(x)
  })
  
  gff_all_random = do.call("rbind",create_indexes)
  #gff_all_random$pos2 = gff_all_random_idx$pos - position
  # gff_all_random_idx$index = loci
  
  gff_all2 = (gff_all_random[,c(1:14)])
  gff_all2$pos = gff_all_random$pos
  regions_overall[[loci]] = gff_all
  regions_cbs[[loci]] = gff_all2
  loci = loci + 1
}

gff_df = gff_df %>% mutate(direction=ifelse(strand == "+",T,F))

#dn_ds_w = merge(dn_ds_w, as.data.frame(gff),by.y="Name",by.x="V1")
regions_overall_df = bind_rows(regions_overall)
regions_overall_df$pos2 = regions_overall_df$pos + 150

regions_cbs_df = bind_rows(regions_cbs)
chrom_xi_cbs = regions_cbs_df %>% filter(grepl("YK",V1)) %>% filter(PGM1 != "DIV/REF")
chrom_ii = regions_overall_df %>% filter(grepl("YB",V1))
p_chrii = chrom_ii %>% ggplot(aes(y=V4,x=pos2, color=GAL1,group=interaction(sample_name,V1))) + geom_line() + 
  theme_bw() + scale_color_manual(values=c("#e41a1c","grey50")) +xlab("") + ylab("Synonymous substitutions per site") + xlim(c(266725,287925))
p3 = gff_df %>% filter(seqnames == "chrII") %>% filter(start > 260725 & end < 290000) %>% ggplot(aes(xmin=start,xmax=end,y=seqnames,forward=direction))  +  geom_gene_arrow(arrow_body_height = unit(20,"mm"),arrowhead_height = unit(24,"mm"))+   geom_gene_label(align = "left",aes(label=gene),grow=T) +
  xlim(c(266725,287925)) + theme_genes() 
pg1 = cowplot::plot_grid(p_chrii + theme(legend.position = "none"),p3,nrow=2,rel_heights = c(1.0,0.3))



chrom_xi = regions_overall_df %>% filter(grepl("YK",V1)) %>% filter(PGM1 != "DIV/REF")

chrom_xi_cbs = regions_cbs_df %>% filter(grepl("YK",V1)) %>% filter(PGM1 != "DIV/REF")
p_chrxi = chrom_xi %>% filter(V4 < 1) %>% ggplot(aes(y=V4,x=pos, color=PGM1,group=interaction(sample_name,V1))) + geom_line() +
  theme_bw() + scale_color_manual(values=c("#377eb8","#e41a1c","grey50")) +xlab("") + ylab("Synonymous substitutions per site") + 
  coord_cartesian(xlim=c(196349,210130)) 

chrom_xi_cbs %>% filter(V4 < 1) %>% ggplot(aes(y=V4,x=pos, color=PGM1,group=interaction(sample_name,V1))) + geom_line() +
  theme_bw() + scale_color_manual(values=c("#377eb8","#e41a1c","grey50")) +xlab("") + ylab("Synonymous substitutions per site") + 
  coord_cartesian(xlim=c(196349,210130)) 
library(gggenes)
gff_df = as.data.frame(gff)
p2 = gff_df %>% filter(seqnames == "chrXI") %>% filter(start > 190000 & end < 220000) %>% ggplot(aes(xmin=start,xmax=end,y=seqnames,forward=direction))  +  geom_gene_arrow(arrow_body_height = unit(20,"mm"),arrowhead_height = unit(24,"mm"))+   geom_gene_label(align = "left",aes(label=gene),grow=T) +
  xlim(c(196349,210130)) + theme_genes() 

pg = cowplot::plot_grid(p_chrxi + theme(legend.position = "none"),p2,nrow=2,rel_heights = c(1.0,0.3))
pg

#legend <- get_legend(
# create some space to the left of the legend
# p_chrxi+ theme(legend.box.margin = margin(0, 0, 0, 12))
#)

cowplot::plot_grid(pg, legend,nrow=2)

chrom_xii = regions_overall_df %>% filter(grepl("YL",V1)) %>% filter(PGM1 != "DIV/REF")
chrom_xii_cbs = regions_cbs_df %>% filter(grepl("YL",V1)) %>% filter(PGM1 != "DIV/REF")


chrom_ii %>% filter(GAL1 == "REF" & V4 > .5)


chrom_xii = regions_overall_df %>%  filter(V1 != "YLR080W") %>% filter(grepl("YL",V1)) %>% filter(PGM1 != "DIV/REF")
p_chrxii = chrom_xii  %>% ggplot(aes(y=V4,x=pos2, color=PGM1,group=interaction(sample_name,V1))) + geom_line() +
  theme_bw() + scale_color_manual(values=c("#377eb8","#e41a1c","grey50")) +xlab("Position") + ylab("Synonymous substitutions per site") +
  coord_cartesian(xlim=c(283872,300251))

p3 = gff_df %>% filter(seqnames == "chrXII") %>% filter(start > 280000 & end < 310000) %>% ggplot(aes(xmin=start,xmax=end,y=seqnames,strand=strand,forward=direction))  +  geom_gene_arrow(arrow_body_height = unit(20,"mm"),arrowhead_height = unit(24,"mm"))+   geom_gene_label(align = "left",aes(label=gene),grow=T) +
  xlim(c(283872,300251)) + theme_genes() 
pg2 = cowplot::plot_grid(p_chrxii + theme(legend.position = "none"),p3,nrow=2,rel_heights = c(1.0,0.3))
#chrom_xii %>% filter(V4 >.4 & GAL2 == "REF" )

svg("figures/s24a.svg")
pg1
dev.off()

svg("figures/s24b.svg")
pg
dev.off()

svg("figures/s24c.svg")
pg2

dev.off()


#plot(regions_overall_df$V4,regions_overall_df$pos)


#gtrack = GenomeAxisTrack()
#gff_yeast_gr$gene = gff_yeast_gr$Name
#gff_yeast_gr$gene = ifelse(is.na(gff_yeast_gr$gene),gff_yeast_gr$Name, gff_yeast_gr$gene)
#grtrack = GeneRegionTrack(gff_yeast_gr,chromosome="chrII",name="Yeast genes",transcriptAnnotation="gene",shape="arrow",fontsize.group=48)
#gff_all_sorted_idx$V1
#chrom_ii_filt = chrom_ii %>% filter(!(V1 %in% gal_genes)) 
#gr = GRanges(seqnames="chrII",ranges=IRanges(start=chrom_ii$pos-1 + 150, end=chrom_ii$pos + 150,names=chrom_ii$V1))
#gr$V2 = chrom_ii$V4

#gr
#datatrack = DataTrack(gr,chromosome="chrII",name="Synonymous divergence (dS)", groups=chrom_ii$GAL1)
#dTrack <- DataTrack(twoGroups, name = "Syno")
#svg("Figures/balancing_gal1.svg",width=10,height=6)
#p1= plotTracks(list(gtrack,grtrack,datatrack),from=266725,to=287925,cex=1,lwd=1, cex.id=20,fontsize=24, col.id=10,fontcolor="black",col.axis="black",col.border.title="black",col.title="black")
#dev.off()
