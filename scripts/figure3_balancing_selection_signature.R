### 
### Figure 3.
library(tidyverse)
library(lemon)
library(rtracklayer)
gff= import.gff3("data/saccharomyces_cerevisiae.gff")
gff = gff[which(gff$orf_classification == "Verified"),]
gal1_10_7_dn = read.table("data/popgen_in/dn_ds/dn_ds//gal1107_region.txt", header=F)
#random_dn = read.table("data/dn_ds/random_region.txt", header=F)
#gdt = read.table("data/dn_ds/gdt.txt", header=F)
#gdt2 = read.table("data/dn_ds/gdt2.txt", header=F)


gal2 = read.table("data/popgen_in/dn_ds/dn_ds//gal2.txt",header=F)
pgm1 = read.table("data/popgen_in/dn_ds/dn_ds//pgm1.txt", header=F)

gal1_10_7_dn = rbind(gal1_10_7_dn,gal2,pgm1)
gal1_10_7_dn$order= 1:nrow(gal1_10_7_dn)

gff_all = merge(gal1_10_7_dn,gff,by.x="V1",by.y="ID")
gff_all_sort = gff_all[order(gff_all$order),]
gff_all_sort$order = 1:nrow(gff_all_sort)
gff_all_sort$V1 = factor(gff_all_sort$V1)

create_indexes = lapply(split(gff_all_sort, f=gff_all_sort$V1),function(x){x$pos = (x$start + seq(from=0,to=length(x$start)*30-1,by=30)) + 300;return(x)})
gff_all_sorted_idx = do.call("rbind",create_indexes)
#gff_all_sorted_idx %>% ggplot(aes(y=V2,x=pos)) + geom_point() + facet_rep_wrap( ~ seqnames,scales="free_x",repeat.tick.labels = T) + ylab("Synonomous per-site mutation (dS)") + xlab("Galactose loci")+  theme_bw()
library(Gviz)

gff_all_sorted_idx = gff_all_sorted_idx[which(gff_all_sorted_idx$gene != "EMP46"),]
# gal genes 
gal_genes = c("GAL2","GAL1","GAL10","GAL7")
gff_all_sorted_idx = gff_all_sorted_idx[which(!(gff_all_sorted_idx$gene %in% gal_genes)),]

gff_yeast_gr = gff[gff$type == "gene" & gff$orf_classification == "Verified",]

gtrack = GenomeAxisTrack()
#gff_yeast_gr$gene = gff_yeast_gr$Name
gff_yeast_gr$gene = ifelse(is.na(gff_yeast_gr$gene),gff_yeast_gr$Name, gff_yeast_gr$gene)
grtrack = GeneRegionTrack(gff_yeast_gr,chromosome="chrII",name="Yeast genes",transcriptAnnotation="gene",shape="arrow",fontsize.group=48)
#gff_all_sorted_idx$V1
gr = GRanges(seqnames=gff_all_sorted_idx$seqnames,ranges=IRanges(start=gff_all_sorted_idx$pos-1 + 150, end=gff_all_sorted_idx$pos + 150,names=gff_all_sorted_idx$V1))
gr$V2 = gff_all_sorted_idx$V2
#gr
datatrack = DataTrack(gr,chromosome="chrII",name="Synonymous divergence (dS)",baseline=0.018)
#dTrack <- DataTrack(twoGroups, name = "Syno")
svg("figures/3a.svg",width=10,height=6)
p1= plotTracks(list(gtrack,grtrack,datatrack),from=266725,to=287925,col="black",cex=1,lwd=1, cex.id=20,fontsize=24, col.id=10,fontcolor="black",col.axis="black",col.border.title="black",col.title="black")
dev.off()

grtrack = GeneRegionTrack(gff_yeast_gr,chromosome="chrXI",name="Yeast genes",transcriptAnnotation="gene",shape="arrow",fontsize.group=48)
gr$V2 = gff_all_sorted_idx$V2
datatrack = DataTrack(gr,chromosome="chrXI",name="Synonymous divergence (dS)",ylim=c(0,.5),baseline=0.018)
svg("figures/3b.svg",width=10,height=6)
p2 = plotTracks(list(gtrack,grtrack,datatrack),from=196349,to=210130,col="black",cex=1,lwd=1, cex.id=20,fontsize=24, col.id=20,fontcolor="black",col.axis="black",col.border.title="black",col.title="black")
dev.off()
gtrack = GenomeAxisTrack()
#gff_yeast_gr$gene = gff_yeast_gr$Name
grtrack = GeneRegionTrack(gff_yeast_gr,chromosome="chrXII",name="Yeast genes",transcriptAnnotation="gene", shape="arrow", fontsize.group=48)
gr$V2 = gff_all_sorted_idx$V2
datatrack = DataTrack(gr,chromosome="chrXII",name="Synonymous divergence (dS)",color="black",baseline=0.018)
svg("figures/3c.svg",width=10,height=6)
p3 = plotTracks(list(gtrack,grtrack,datatrack),from=283872,to=300251,col="black",cex=1,lwd=1, cex.id=20,fontsize=24, col.id=20,fontcolor="black",col.axis="black",col.border.title="black",col.title="black")
dev.off()


