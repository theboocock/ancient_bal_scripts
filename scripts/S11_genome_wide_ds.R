dn_ds_w = read.table("data/popgen_in/dn_ds/dn_ds_w.txt")
gff= import.gff3("data/saccharomyces_cerevisiae.gff")
gff = gff[which(gff$orf_classification == "Verified" & gff$type =="gene"),]
dn_ds_w = merge(dn_ds_w, as.data.frame(gff),by.y="Name",by.x="V1")
library(plyranges)
ranges = read.table("data/popgen_in/dn_ds/dn_ds//full_windows.txt", header=F)
colnames(ranges) = c("seqnames","start","end")
ranges_gr = ranges %>% as_granges()
dn_ds_gr= dn_ds_w %>% as_granges()
dn_ds_gr$order = 1:length(dn_ds_gr)
row_meta = as.data.frame(dn_ds_gr)
dn_ds_gr$start_gene  = row_meta$start
dn_ds_gr$end_gene = row_meta$end
dn_ds_gr$seqnames_gene = row_meta$seqnames
create_indexes = lapply(split(as.data.frame(dn_ds_gr), f=factor(dn_ds_gr$V1)),function(x){
  #print(x)
  #print(seq(from=0,to=length(x$start)*30 - 1,by=30));
  #print(x$start + seq(from=0,to=length(x$start)*30-1,by=30));
  x$pos = x$start_gene + seq(from=0,to=length(x$start_gene)*30-1,by=30);
  ##print(x$pos)
  #print(x$pos);  print("HERE");print(x);return(x)}
  return(x)
})


dn_ds_gr_two  = do.call("rbind",create_indexes)
gal_genes = c("YBR018C","YBR019C","YBR020W","YLR081W","YLR080W")
flankers = c("YKL130C",
             "YKL129C",
             "YKL128C",
             "YKL126W",
             "YKL125W",
             "YBR014C",
             "YKL127W",
             "YBR015C",
             "YBR016W",
             "YBR017C",
             "YBR021W",
             "YBR022W",
             "YBR023C",
             "YLR079W",
             "YLR080W",
             "YLR082C",
             "YLR083C",
             "YLR084C")
dn_ds_gr_two$gal_genes = rep(NA,nrow(dn_ds_gr_two))
dn_ds_gr_two$gal_genes[dn_ds_gr_two$V1 %in% gal_genes] = "Galalactose genes"
dn_ds_gr_two$gal_genes[dn_ds_gr_two$V1 %in% flankers] = "Flanking genes"
dn_ds_gr_two$gal_genes[is.na(dn_ds_gr_two$gal_genes)] = "Other genes"
dn_ds_gr_two = dn_ds_gr_two %>% filter(seqnames != "chrmt")
png("figures/S11.png", width=12*200,height=16 * 200)
dn_ds_gr_two %>% filter(V2 == "WINDOW") %>% ggplot(aes(y=(V4),x=pos,color=gal_genes)) + geom_point() + facet_wrap(~seqnames,nrow=16) + ylab("Synonymous substitutions per site") + theme_bw() +
  scale_color_manual(values=c("#e41a1c","#377eb8","grey50"))
dev.off()

