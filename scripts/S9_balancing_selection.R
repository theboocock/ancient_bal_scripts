###
library(tidyverse)
library(lemon)
library(Gviz)
library(rtracklayer)
gff= import.gff3("data/saccharomyces_cerevisiae.gff")
gff = gff[which(gff$orf_classification == "Verified"),]


window_sizes_all = list.files("data/balancing_selection/windowless/",pattern="*txt2",full.names = T)
plots = list()
i = 1
gffs= list()
# Create Figure S9
for (window in window_sizes_all){
  print(window)
  
  window_in = read.table(window, header=F)
  print(dim(window_in))
  
  window_in = window_in %>% filter(window_in$V2 != "OVERALL")
  gff_all = merge(window_in,gff,by.x="V1",by.y="ID")
  gff_all$order= 1:nrow(gff_all)
  
  gff_all_sort = gff_all[order(gff_all$order),]
  gff_all_sort$order = 1:nrow(gff_all_sort)
  gff_all_sort$V1 = factor(gff_all_sort$V1)
  
  window_str = as.numeric(strsplit(unlist(strsplit(basename(window),".txt2")),"_")[[1]][3])
  step = as.numeric(strsplit(unlist(strsplit(basename(window),".txt2")),"_")[[1]][4])
  if(window_str != step){
    next
  }
  if(window_str != "75"){
    next
  }
  
  create_indexes = lapply(split(gff_all_sort, f=gff_all_sort$V1),function(x){x$pos = x$start + seq(from=0,to=length(x$start)*step*3-1,by=step*3);return(x)})
  gff_all_sorted_idx = do.call("rbind",create_indexes)
  gal_genes = c("GAL2","GAL1","GAL10","GAL7")
  gff_all_sorted_idx = gff_all_sorted_idx[which(!(gff_all_sorted_idx$gene %in% gal_genes)),]
  
  gffs[[i]] = gff_all_sorted_idx
  
  
  
  plots[[i]] = (gff_all_sorted_idx %>% ggplot(aes(y=V4,x=pos)) + geom_point() + facet_rep_wrap( ~ seqnames,scales="free_x",repeat.tick.labels = T) + 
                  ylab("Synonomous per-site mutation (dS)") + xlab("Galactose loci")+  theme_bw() + ggtitle(paste("window: ",window_str,"step:", step,sep=" ")))
  
  
  
  gff_yeast_gr = gff[gff$type == "gene" & gff$orf_classification == "Verified",]
  
  gtrack = GenomeAxisTrack()
  #gff_yeast_gr$gene = gff_yeast_gr$Name
  gff_yeast_gr$gene = ifelse(is.na(gff_yeast_gr$gene),gff_yeast_gr$Name, gff_yeast_gr$gene)
  grtrack = GeneRegionTrack(gff_yeast_gr,chromosome="chrII",name="Yeast genes",transcriptAnnotation="gene",shape="arrow",fontsize.group=48)
  #gff_all_sorted_idx$V1
  gr = GRanges(seqnames=gff_all_sorted_idx$seqnames,ranges=IRanges(start=gff_all_sorted_idx$pos-1 + 150, end=gff_all_sorted_idx$pos + 150,names=gff_all_sorted_idx$V4))
  gr$V4 = gff_all_sorted_idx$V4
  #gr
  datatrack = DataTrack(gr,chromosome="chrII",name="Synonymous divergence (dS)",baseline=0.018)
  #dTrack <- DataTrack(twoGroups, name = "Syno")
  #svg("Figures/balancing_gal1.svg",width=10,height=6)
  
  p1l = list(gtrack,grtrack,datatrack)
  
  #p1= plotTracks(list(gtrack,grtrack,datatrack),from=266725,to=287925,col="black",cex=1,lwd=1, cex.id=20,fontsize=24, col.id=10,fontcolor="black",col.axis="black",col.border.title="black",col.title="black")
  #dev.off()
  
  
  #svg("figures/gal7.svg",width=10,height=6)
  #p1= plotTracks(list(gtrack,grtrack,datatrack),from=266725,to=275000,col="black",cex=1,lwd=1, cex.id=21,fontsize=24, col.id=20,fontcolor="black",col.axis="black",col.border.title="black",col.title="black")
  #dev.off()
  gtrack = GenomeAxisTrack()
  gtrack = GenomeAxisTrack()
  grtrack = GeneRegionTrack(gff_yeast_gr,chromosome="chrXI",name="Yeast genes",transcriptAnnotation="gene",shape="arrow",fontsize.group=48)
  gr$V4 = gff_all_sorted_idx$V4
  datatrack = DataTrack(gr,chromosome="chrXI",name="Synonymous divergence (dS)",baseline=0.018)
  p2l = list(gtrack,grtrack,datatrack)
  
  
  #p2 = plotTracks(list(gtrack,grtrack,datatrack),from=196349,to=210130,col="black",cex=1,lwd=1, cex.id=20,fontsize=24, col.id=20,fontcolor="black",col.axis="black",col.border.title="black",col.title="black")
  #dev.off()
  gtrack = GenomeAxisTrack()
  #gff_yeast_gr$gene = gff_yeast_gr$Name
  grtrack = GeneRegionTrack(gff_yeast_gr,chromosome="chrXII",name="Yeast genes",transcriptAnnotation="gene", shape="arrow", fontsize.group=48)
  gr$V4 = gff_all_sorted_idx$V4
  datatrack = DataTrack(gr,chromosome="chrXII",name="Synonymous divergence (dS)",color="black",baseline=0.018)
  
  p3l = list(gtrack,grtrack,datatrack)
  
  #p3 = plotTracks(list(gtrack,grtrack,datatrack),from=283872,to=300251,col="black",cex=1,lwd=1, cex.id=20,fontsize=24, col.id=20,fontcolor="black",col.axis="black",col.border.title="black",col.title="black")
  #dev.off()
  
  
  
  library(grid)
  svg(paste("figures/","S9.svg",sep=""),width=10,height=12)
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(3, 1)))
  pushViewport(viewport(layout.pos.col=1,layout.pos.row=1))
  
  plotTracks(p1l,from=266725,to=287925,col="black",cex=1,lwd=1, cex.id=20,fontsize=24, col.id=20,fontcolor="black",col.axis="black",col.border.title="black",col.title="black",add=T)
  popViewport()
  pushViewport(viewport(layout.pos.col=1,layout.pos.row=2))
  plotTracks(p2l,from=196349,to=210130,col="black",cex=1,lwd=1, cex.id=20,fontsize=24, col.id=20,fontcolor="black",col.axis="black",col.border.title="black",col.title="black",add=T)
  
  popViewport()
  pushViewport(viewport(layout.pos.col=1,layout.pos.row=3))
  plotTracks(p3l,from=283872,to=300251,col="black",cex=1,lwd=1, cex.id=20,fontsize=24, col.id=20,fontcolor="black",col.axis="black",col.border.title="black",col.title="black",add=T)
  popViewport()
  dev.off()
  i = i + 1
  
}

# ### Not sure what this is for ###
# #library(tidyverse)
# #library(lemon)
# #library(rtracklayer)
# gff= import.gff3("data//saccharomyces_cerevisiae.gff")
# gff = gff[which(gff$orf_classification == "Verified"),]
# 
# 
# window_sizes_all = list.files("popgen_phylo_notebooks/outputs/windowed/",pattern="*pmu1",full.names = T)
# plots2 = list()
# i = 1
# gffs= list()
# x = 1
# for (window in window_sizes_all){
#   print(window)
#   
#   window_in = read.table(window, header=F)
#   print(dim(window_in))
#   
#   window_in = window_in %>% filter(window_in$V2 != "OVERALL")
#   print(window_in)
#   gff_all = merge(window_in,gff,by.x="V1",by.y="ID")
#   gff_all$order= 1:nrow(gff_all)
#   
#   gff_all_sort = gff_all[order(gff_all$order),]
#   gff_all_sort$order = 1:nrow(gff_all_sort)
#   gff_all_sort$V1 = factor(gff_all_sort$V1)
#   
#   window_str = as.numeric(strsplit(unlist(strsplit(basename(window),".pmu1")),"_")[[1]][3])
#   step = as.numeric(strsplit(unlist(strsplit(basename(window),".pmu1")),"_")[[1]][4])
#   if(window_str != step){
#     x = x + 1
#     next
#   }else{
#     j = x
#     x =x +1
#   }
#   # j = i
#   
#   create_indexes = lapply(split(gff_all_sort, f=gff_all_sort$V1),function(x){x$pos = x$start + seq(from=0,to=length(x$start)*step*3-1,by=step*3);return(x)})
#   gff_all_sorted_idx = do.call("rbind",create_indexes)
#   
#   gal_genes = c("GAL2","GAL1","GAL10","GAL7")
#   gff_all_sorted_idx = gff_all_sorted_idx[which(!(gff_all_sorted_idx$gene %in% gal_genes)),]
#   
#   gffs[[i]] = gff_all_sorted_idx
#   
#   
#   
#   plots2[[i]] = (gff_all_sorted_idx %>% ggplot(aes(y=V4,x=pos)) + geom_point() + facet_rep_wrap( ~ seqnames,scales="free_x",repeat.tick.labels = T) + 
#                    ylab("Synonomous per-site mutation (dS)") + xlab("Galactose loci")+  theme_bw() + ggtitle(paste("window: ",window_str,"step:", step,sep=" ")))
#   
#   i =i +1
# }
# 
# 
