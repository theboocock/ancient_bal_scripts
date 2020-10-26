### Read data stats ###
## S17 ##
library(tidyverse)
library(lemon)
library(rtracklayer)
library(plyranges)
window_step = 30
genome_width = 600
library(GenomicRanges)
tile_genome = GRanges(seqnames="chr1",IRanges(c(seq(1,5000-genome_width + 1,by=window_step),seq(5001,10000-genome_width + 1,by=window_step),
                                                seq(10001,15000-genome_width +1,by=window_step),seq(15001,20000-genome_width +1,by=window_step),
                                                seq(20001,25000-genome_width +1,by=window_step), seq(25001,30000-genome_width+1,by=window_step)),
                                              end=c(seq(1 -1 + genome_width,5000,by=window_step),seq(5001 - 1 + genome_width,10000,by=window_step),
                                                    seq(10001 -1 + genome_width, 15000,by=window_step),seq(15001-1+genome_width,20000,by=window_step),
                                                    seq(20000 -1 + genome_width,25000,by=window_step), seq(25001-1 + genome_width,30000,by=window_step))))



#tile_g2 = GRanges(seqnames="chr1",IRanges(1,5001,10001,15001,20001,250001,end=c(5000,10000,15000,20000,25000,30000)))
gff= import.gff3("data/saccharomyces_cerevisiae.gff")
gff = gff[which(gff$orf_classification == "Verified"),]
gal1_10_7_dn = read.table("data/popgen_in/dn_ds/dn_ds///gal1107_region.txt", header=F)
#random_dn = read.table("~/Dropbox/PHDTHESIS/projects/gal_final_github_october/data/dn_ds/random_region.txt", header=F)
#gdt = read.table("~/Dropbox/PHDTHESIS/projects/gal_final_github_october/data/dn_ds/gdt.txt", header=F)
#gdt2 = read.table("~/Dropbox/PHDTHESIS/projects/gal_final_github_october/data/dn_ds/gdt2.txt", header=F)
gal2 = read.table("data/popgen_in/dn_ds/dn_ds//gal2.txt",header=F)
pgm1 = read.table("data/popgen_in/dn_ds/dn_ds/pgm1.txt", header=F)

gal1_10_7_dn = rbind(gal1_10_7_dn,gal2,pgm1)
gal1_10_7_dn$order= 1:nrow(gal1_10_7_dn)

gff_all = merge(gal1_10_7_dn,gff,by.x="V1",by.y="ID")
gff_all_sort = gff_all[order(gff_all$order),]
gff_all_sort$order = 1:nrow(gff_all_sort)
gff_all_sort$V1 = factor(gff_all_sort$V1)

create_indexes = lapply(split(gff_all_sort, f=gff_all_sort$V1),function(x){
  x$pos = x$start + seq(from=0 + 300,to=length(x$start)*30-1 + 300,by=30); 
  x$start2 = x$start + seq(from=0 ,to=length(x$start)*30-1,by=30) +1;
  x$end2 =  x$start + seq(from=0 ,to=length(x$start)*30-1,by=30) + 600;
  return(x)})
gff_all_sorted_idx = do.call("rbind",create_indexes)
gal_genes = c("GAL2","GAL1","GAL10","GAL7","EMP46")
gff_all_sorted_idx = gff_all_sorted_idx[which(!(gff_all_sorted_idx$gene %in% gal_genes)),]
g1 = gff_all_sorted_idx %>% filter(seqnames == "chrII") %>% mutate(new=ifelse(pos < 274427,T,F)) %>% mutate(pos_new=ifelse(new,pos-274427 + 5000,pos-280607 + 5000),
                                                                                                            start2=ifelse(new,start2-274427 + 5000,start2-280607 + 5000),
                                                                                                            end2=ifelse(new,end2-274427 + 5000,end2-280607 + 5000)) %>%  filter(pos_new >  0 & pos_new < 10000)

# %>% ggplot(aes(x=pos_new,y=V2))  + geom_point() +
#  geom_vline(xintercept = 5000) + xlim(c(0,10000)) + theme_bw()  + xlab("Position") + ylab("DS") + ggtitle("GAL1/10/7")
#gff_all_sorted_idx %>% filter(seqnames == "chrII") %>% mutate(new=ifelse(pos < 274427,T,F)) %>% mutate(pos_new=ifelse(new,pos-274427 + 5000,pos-280607 + 5000)) %>% 
#  filter(pos_new > 0  & pos_new<  10000) %>% ggplot(aes(x=pos_new,y=(V2)))  + geom_point() + geom_vline(xintercept = 5000) + theme_bw()
### Build two regressionn lines #### 


gff_all_sorted_idx = gff_all_sorted_idx[which(!(gff_all_sorted_idx$gene %in% gal_genes)),]
start_pgm1 = 201772
end_pgm1 = 203540
g2 = gff_all_sorted_idx %>% filter(seqnames == "chrXI") %>% mutate(new=ifelse(pos < start_pgm1,T,F)) %>% mutate(pos_new=ifelse(new,pos-start_pgm1 + 15000,pos-end_pgm1 + 15000),
                                                                                                                start2=ifelse(new,start2-start_pgm1 + 15000,start2-end_pgm1 + 15000),
                                                                                                                end2=ifelse(new,end2-start_pgm1 + 15000,end2-end_pgm1 + 15000)) %>% filter(pos_new >  10000 & pos_new < 20000)# %>% ggplot(aes(x=pos_new,y=V2))  + geom_point() + geom_vline(xintercept = 15000)

start_gal2 = 290212
end_gal2 = 291936


g3 = gff_all_sorted_idx %>% filter(seqnames == "chrXII") %>% mutate(new=ifelse(pos < start_gal2,T,F)) %>% mutate(pos_new=ifelse(new,pos-start_gal2 + 25000,pos-end_gal2 + 25000),
                                                                                                                 start2=ifelse(new,start2-start_gal2 + 25000,start2-end_gal2 + 25000),
                                                                                                                 end2=ifelse(new,end2-start_gal2 + 25000,end2-end_gal2 + 25000)) %>% filter(pos_new >  20000 & pos_new < 30000) #%>% ggplot(aes(x=pos_new,y=V2))  + geom_point() + geom_vline(xintercept = 25000)

#start_gal2 
#end_gal2 
#lm(gff_all_sorted_idx)
b1 = bind_rows(g1,g2,g3)
bind_rows(g1,g2,g3) %>% ggplot(aes(x=pos_new,y=V2))  + geom_point() +
  geom_vline(xintercept = c(5000,15000,25000)) + theme_bw()  + xlab("Position") + ylab("dS")

g = GRanges(seqnames="chr1",IRanges(b1$pos_new,end = b1$pos_new + 1))
#g$value = b1$V2
df2 = as.data.frame(findOverlaps(g,tile_genome,minoverlap = 0)) %>% group_by(queryHits) %>% summarise(max_hit=ceiling(mean((subjectHits))))
df2$value = b1$V2[df2$queryHits]
#df2 = as.data.frame(findOverlaps(g,tile_genome,minoverlap = 1))
#df2as.data.frame(b1$)

#df = as.data.frame(seq_id=1:nrow())
g = tile_genome %>% as.data.frame
g$V1 = NA
g$V1[df2$max_hit] = df2$value


#g$Value = b1$V2
tile_g2 = GRanges(seqnames="chr1",IRanges(c(1,5001,10001,15001,20001,25001),end=c(5000,10000,15000,20000,25000,30000)))
### Calculate eLD ###
#tile_genome %>% 
df_out = tile_g2 %>% group_by_overlaps(tile_genome) %>% as.data.frame
#pi_between_df$query = NULL
df_out = cbind(df_out,as.data.frame(g))
out_df = data.frame()
df_out$start2 = df_out[,8]


for (i in unique(df_out$query)){
  df2 = df_out[df_out$query == i,]
  df2$row_number = 1:nrow(df2)
  
  df2 = df2[50:147,]
  m1 = ((lm(log(V1 + 0.0001)~ row_number,data=df2)))
  intercept= (summary(m1)$coef[1,1])
  slope = (summary(m1)$coef[2,1])
  
  max_ds = max(df2$V1,na.rm=T)
  p = data.frame(row_number=max(df2$row_number))
  predict1 = (predict(m1,p))
  #print(exp(predict1))
  #print(max_ds)
  tmp_row = data.frame(max_ds=max_ds,i=i, slope=slope,predict=exp(predict1),intercept=intercept)
  out_df = rbind(out_df, tmp_row)
}

# Real data 
epsilon = 0.5916645

#data.frame(start=df_out$start2, value=df_out$V1)
#y = log(intro$pi_between_dfs$`785.out.sim.rds`$pi_between)

df_real_data = data.frame(window_position=df_out$start2 + 300, value=df_out$V1)
#smooth.spline()
svg("figures/S17.svg")
df_real_data %>% ggplot(aes(x=window_position,y=value)) + geom_point() + geom_vline(xintercept = c(5000),color="#e41a1c") + 
  geom_vline(xintercept = c(15000),color="#377eb8") + geom_vline(xintercept = c(25000),color="#984ea3")+ theme_bw() +
  xlab("Position") + ylab("Synonymous rate of divergence (dS)") + xlim(c(0,30000)) + theme(text=element_text(size=20))
dev.off()
