### S23, S6, and 4C.

library(tidyverse)
in_files = list.files("data/ase/",pattern = "*ase",full.names = T)
res_list = list()
j  = 1
### ASE analyses.

for(in_f in in_files){
  by_rm_ase = read.delim(in_f, header=F,sep="\t")
  colnames(by_rm_ase) = c("chrom_gene","gene_start","gene_end","gene_name","chrom_snp","start","end","chrom_snp2","snp_pos","snp_id","ref","alt","refc","altc","totalcount","lowMAPqdepth","lowBASEqdepth","rawdepth","otherbases","improperpairs")
  by_genes = by_rm_ase %>% group_by(gene_name)
  spli_df = split(by_rm_ase,by_rm_ase$gene_name)
  summarise_counts = by_genes %>% summarise(alt_c=sum(altc),ref_c=sum(refc)) %>% as.data.frame
  summarise_counts  = summarise_counts[!(summarise_counts$alt_c + summarise_counts$ref_c) == 0,]
  p_values = rep(NA,nrow(summarise_counts))
  for(i in 1:nrow(summarise_counts)){
    p_values[i] = binom.test(summarise_counts$alt_c[i],summarise_counts$alt_c[i] + summarise_counts$ref_c[i])$p.value
  }
  summarise_counts$log_two = log2(summarise_counts$alt_c)  - log2(summarise_counts$ref_c)
  summarise_counts$p_values = p_values
  res_list[[j]] = summarise_counts
  j = j + 1
}

names(res_list) = c("Mid-log","0.5 Hour", "1 Hour","2 Hours","4 Hours", "5.5 Hours")

diff_exp = res_list[[3]] %>% filter(alt_c > 5 & ref_c > 5) %>% mutate(p_adj=p.adjust(p_values,method="fdr")) %>% filter(p_adj < 0.05)

all_genes = bind_rows(res_list,.id="condition") %>% filter(alt_c > 5 & ref_c > 5) %>% mutate(p_adj=p.adjust(p_values,method="fdr")) %>% filter(p_adj < 0.05 & abs(log_two) > 2)

colnames(all_genes) = c("Growth condition","Gene name","Alternate count","Reference count","Log2 allelic fold-change","P-value (unadj)","P-value (adj)")
write.csv(all_genes,file="figures/table7a.csv")
library(rtracklayer)
yeast_gff3= import.gff3("data/saccharomyces_cerevisiae.gff")
yeast_gff3 = as.data.frame(yeast_gff3)
yeast_verified = yeast_gff3[yeast_gff3$type == "gene" & yeast_gff3$orf_classification == "Verified",]
gene_list = c("YBR017C","YKL127W","YKL128C","YBR021W","YKL035W","YIL099W","YDR074W","YKL128C")

ab = lapply(res_list,function(x){x[x$gene_name %in% gene_list,]})
in_timecourse = list.files("data/ase/",pattern = "*tsv", full.names = T)
#gal_genes = file_in[file_in$tpm>= 1 ,]

all_genes$order = 1:nrow(all_genes)
aaa = merge(all_genes, as.data.frame(yeast_verified),by.y="ID", by.x=2)
aaa = merge(all_genes, as.data.frame(yeast_verified),by.y="ID", by.x=2)
aaa = aaa[order(aaa$order),]
aaa =cbind(aaa[,c(1:7)] ,aaa$gene)
#aaa[order(aaa$)]
#write.csv(aaa,file="figures/table7a.csv")

genes = c("gal2a_cbs","gal2b_cbs","YBR018C_cbs_cbs","YBR020W_cbs_cbs","YBR019C_cbs_cbs","YBR018C","YBR020W","YBR019C","YLR081W")

all_gene_list = data.frame()
all_gene_list_overtime = data.frame()
all_gene_errors = data.frame()
pgm1_time_points = data.frame()


out_df = data.frame()
i = 1
for(file in in_timecourse){
  file_in = read.delim(file, sep="\t", header=T,stringsAsFactors = F)
  file_verified = file_in[file_in$target_id %in% yeast_verified$Name | file_in$target_id %in% genes,]
  file_in =file_in[file_in$target_id %in% genes,]
  file_in$time_point = i
  file_in$gal2_cbs = file_in$tpm[1] + file_in$tpm[2]
  file_in  = rbind(file_in,file_in[1,])
  #file_in$tpm[10] = file_in$tpm[1] + file_in$tpm[2]
  #f#ile_in$target_id[10] = "gal2_cbs"
  time_point = i 
  #gene = "YBR018C"
  
  
  gal2a = log2(file_in$est_counts[1]) -log2(file_in$est_counts[9])
  gal2b = log2(file_in$est_counts[2]) - log2(file_in$est_counts[9])  
  gal7 = log2(file_in$est_counts[3]) - log2(file_in$est_counts[6])
  gal10 =  log2(file_in$est_counts[4]) - log2(file_in$est_counts[8])
  gal1 = log2(file_in$est_counts[5]) - log2(file_in$est_counts[7])
  pgm1_fold_change = log2(ab[[i]]$alt_c[6]) - log2(ab[[i]]$ref_c[6])
  
  
  gal2a_p = binom.test(as.integer(file_in$est_counts[1]),as.integer(file_in$est_counts[1]) + as.integer(file_in$est_counts[9]))$p.value
  gal2b_p= binom.test(as.integer(file_in$est_counts[2]), as.integer(file_in$est_counts[9]) + as.integer(file_in$est_counts[2]))$p.value
  gal7_p = binom.test(as.integer(file_in$est_counts[3]), as.integer(file_in$est_counts[3]) + as.integer(file_in$est_counts[6]))$p.value
  gal10_p = binom.test(as.integer(file_in$est_counts[4]), as.integer(file_in$est_counts[8]) + as.integer(file_in$est_counts[4]))$p.value
  gal1_p = binom.test(as.integer(file_in$est_counts[5]), as.integer(file_in$est_counts[5]) + as.integer(file_in$est_counts[7]))$p.value
  pgm1_p = binom.test(ab[[i]]$alt_c[6], ab[[i]]$ref_c[6] + ab[[i]]$alt_c[6])$p.value
  
  tmp_row =c(i, "GAL1",file_in$est_counts[5],file_in$est_counts[7],gal1,gal1_p)
  out_df= rbind(out_df, tmp_row)
  
  tmp_row =c(i, "GAL10",file_in$est_counts[4],file_in$est_counts[8],gal10,gal10_p)
  out_df= rbind(out_df, tmp_row)
  
  tmp_row =c(i, "GAL7",file_in$est_counts[3],file_in$est_counts[6],gal7,gal7_p)
  out_df= rbind(out_df, tmp_row)
  
  tmp_row =c(i, "GAL2a",file_in$est_counts[1],file_in$est_counts[9],gal2a,gal2a_p)
  out_df= rbind(out_df, tmp_row)
  
  tmp_row =c(i, "GAL2b",file_in$est_counts[2],file_in$est_counts[9],gal2b,gal2b_p)
  out_df= rbind(out_df, tmp_row)
  
  
  i = i + 1
}
#write.csv(out_df,"figures/s5b.csv")
i = 1
for(file in in_timecourse){
  file_in = read.delim(file, sep="\t", header=T,stringsAsFactors = F)
  file_verified = file_in[file_in$target_id %in% yeast_verified$Name | file_in$target_id %in% genes,]
  file_in =file_in[file_in$target_id %in% genes,]
  file_in$time_point = i
  file_in$gal2_cbs = file_in$tpm[1] + file_in$tpm[2]
  file_in  = rbind(file_in,file_in[1,])
  file_in$tpm[10] = file_in$tpm[1] + file_in$tpm[2]
  file_in$target_id[10] = "gal2_cbs"
  
  gal2a = log2(file_in$est_counts[1]) -log2(file_in$est_counts[9])
  gal2b = log2(file_in$est_counts[2]) - log2(file_in$est_counts[9])  
  gal7 = log2(file_in$est_counts[3]) - log2(file_in$est_counts[6])
  gal10 =  log2(file_in$est_counts[4]) - log2(file_in$est_counts[8])
  gal1 = log2(file_in$est_counts[5]) - log2(file_in$est_counts[7])
  pgm1_fold_change = log2(ab[[i]]$alt_c[6]) - log2(ab[[i]]$ref_c[6])
  
  
  
  gal2b_p = binom.test(as.integer(file_in$est_counts[1]),as.integer(file_in$est_counts[1]) + as.integer(file_in$est_counts[9]))$p.value
  gal2b_p= binom.test(as.integer(file_in$est_counts[2]), as.integer(file_in$est_counts[9]) + as.integer(file_in$est_counts[2]))$p.value
  gal7_p = binom.test(as.integer(file_in$est_counts[3]), as.integer(file_in$est_counts[3]) + as.integer(file_in$est_counts[6]))$p.value
  gal10_p = binom.test(as.integer(file_in$est_counts[4]), as.integer(file_in$est_counts[8]) + as.integer(file_in$est_counts[4]))$p.value
  gal1_p = binom.test(as.integer(file_in$est_counts[5]), as.integer(file_in$est_counts[5]) + as.integer(file_in$est_counts[7]))$p.value
  pgm1_p = binom.test(ab[[i]]$alt_c[6], ab[[i]]$ref_c[6] + ab[[i]]$alt_c[6])$p.value
  
  gal10 =  log2(file_in$est_counts[4]) - log2(file_in$est_counts[8])
  gal1 = log2(file_in$est_counts[5]) - log2(file_in$est_counts[7])
  pgm1_fold_change = log2(ab[[i]]$alt_c[6]) - log2(ab[[i]]$ref_c[6])
  
  
  
  # sigma^2_f \approx f^2[(sigma^2A/A) + \sigma^2_b/b + 2 (cov(ab)/E(AB))]
  #https://en.wikipedia.org/wiki/Propagation_of_uncertainty
  gal2a_se = sqrt((file_in$est_counts[1]/file_in$est_counts[9])^2 *(file_in$est_counts[1]/file_in$est_counts[1]^2 +
                                                                      file_in$est_counts[9]/file_in$est_counts[9]^2))/ (1/log(2)*file_in$est_counts[1]/file_in$est_counts[9])
  gal2b_se = sqrt((file_in$est_counts[2]/file_in$est_counts[9])^2 *(file_in$est_counts[2]/file_in$est_counts[2]^2 + file_in$est_counts[9]/file_in$est_counts[9]^2))/(1/log(2)*file_in$est_counts[2]/file_in$est_counts[9])
  gal7_se = sqrt((file_in$est_counts[3]/file_in$est_counts[6])^2 *(file_in$est_counts[3]/file_in$est_counts[3]^2 + file_in$est_counts[6]/file_in$est_counts[6]^2))/(1/log(2)*file_in$est_counts[3]/file_in$est_counts[6])
  gal10_se = sqrt((file_in$est_counts[4]/file_in$est_counts[8])^2 *(file_in$est_counts[4]/file_in$est_counts[4]^2 + file_in$est_counts[8]/file_in$est_counts[8]^2))/(1/log(2)*file_in$est_counts[4]/file_in$est_counts[8])
  gal1_se = sqrt((file_in$est_counts[5]/file_in$est_counts[7])^2 *(file_in$est_counts[5]/file_in$est_counts[5]^2 + file_in$est_counts[8]/file_in$est_counts[8]^2))/(1/log(2)*file_in$est_counts[5]/file_in$est_counts[7])
  #print(pgm1_se)
  # gal2=  log2(file_in$est_counts[1] + file_in$est_counts[2]) -log2(file_in$est_counts[9]
  # gal2_se = 
  
  pgm1_se = sqrt((ab[[i]]$alt_c[6]/ab[[i]]$ref_c[6])^2 *(1/(ab[[i]]$alt_c[6]) + 1/(ab[[i]]$ref_c[6])))/(1/log(2)*ab[[i]]$alt_c[6]/ab[[i]]$ref_c[6])
  #print((ab[[i]]$alt_c[6]/ab[[i]]$ref_c[6])^2)
  #print((1/log(2)*ab[[i]]$alt_c[6]/ab[[i]]$ref_c[6]))
  #file_
  #print(pgm1_se)
  tmp_row = data.frame(time_point=i, gal2a=gal2a,gal2b=gal2b, gal7=gal7,gal10=gal10,gal1=gal1,pgm1=pgm1_fold_change)
  errors =  data.frame(time_point=i, gal2a=gal2a_se,gal2b=gal2b_se, gal7=gal7_se,gal10=gal10_se,gal1=gal1_se,pgm1=pgm1_se)
  # counts_row_one = data.frame()
  
  
  pgm1_time_points_row = data.frame(time_point=i, pgm1_alt=ab[[i]]$alt_c[6],pgm1_ref=ab[[i]]$ref_c[6])
  pgm1_time_points = rbind(pgm1_time_points,pgm1_time_points_row)
  all_gene_errors = rbind(all_gene_errors, errors)
  all_gene_list = rbind(all_gene_list, tmp_row)
  all_gene_list_overtime=rbind(all_gene_list_overtime, file_in)
  i = i +1
}

all_gene_list = melt(all_gene_list,id.vars = c("time_point"))
all_gene_errors = melt(all_gene_errors,id.vars=c("time_point"))
all_gene_list = merge(all_gene_list,all_gene_errors,by=c(1,2))
#all_gene_list_se = all_gene_list  %>% filter(grepl("_se",variable))
#all_gene_list_dat = all_gene_list  %>% filter(!(grepl("_se",variable)))


library(ggplot2)
library(dplyr)

#(file_in$est_counts[1]/file_in$tpm[9])^2 *(file_in$tpm[1]/file_in$tpm[1]^2 + file_in$tpm[9]/file_in$tpm[9]^2)
### table s5b
#all_gene_list_overtime

write.table(file="figures/tableS5b.csv")

#all_gene_list %>% ggplot(aes(y=(variable),x=time_point)) + geom_tile(aes(fill=value.x)) + xlab("Time point") + ylab("Gene name") +  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-8,8), space = "Lab", name="Log2FoldChange") 

all_gene_list$time_point[all_gene_list$time_point == 1] = "Mid-log"
all_gene_list$time_point[all_gene_list$time_point == 2] = "0.5 hour"
all_gene_list$time_point[all_gene_list$time_point == 3] = "1 hour"
all_gene_list$time_point[all_gene_list$time_point == 4] = "2 hours"
all_gene_list$time_point[all_gene_list$time_point == 5] = "4 hours"
all_gene_list$time_point[all_gene_list$time_point == 6] = "5.5 hours"
all_gene_list$glucose = ifelse(all_gene_list$time_point == "Mid-log","Glucose","Galactose")
all_gene_list$glucose = factor(all_gene_list$glucose,levels=c("Glucose","Galactose"))
all_gene_list$time_point = factor(all_gene_list$time_point,levels=c("Mid-log","0.5 hour","1 hour","2 hours","4 hours","5.5 hours"))

all_gene_list_overtime$by = "CBS2888"
all_gene_list_overtime$by[grepl("Y",all_gene_list_overtime$target_id) & !grepl("cbs",all_gene_list_overtime$target_id)] = "BY"
#all_gene_list_overtime$by[all_gene_list_overtime$by] #
p2 = all_gene_list %>% mutate(variable = toupper(variable)) %>%  ggplot(aes(y=(variable),x=time_point)) + geom_tile(aes(fill=value.x)) + 
  xlab("Time point") + ylab("Gene name") +  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-8,8), space = "Lab", name="Log2FoldChange")  + 
  facet_grid(~glucose,space="free",scales="free_x") + theme(text=element_text(size=20))
#all_gene_list%>% mutate(variable = toupper(variable)) %>%  ggplot(aes(y=(value.x),x=time_point,group=variable)) + geom_errorbar(aes(ymin=value.x - 2 *value.y, ymax=value.x + 2*value.y),width=0.2) + geom_line() + geom_point() + xlab("Time point") + ylab("Log2FC") + theme_bw()  + expand_limits(y=0) + facet_wrap(~variable)

p2 = all_gene_list %>% filter(glucose == "Glucose") %>% mutate(variable = toupper(variable)) %>%  ggplot(aes(y=value.x,x=variable)) + geom_bar(stat="identity")  + theme_bw() +
  xlab("Gene") + ylab("Log2 Fold-change") +geom_errorbar(aes(ymin=value.x - 2 *value.y, ymax=value.x + 2*value.y),width=0.2)  + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(text=element_text(size=20))
svg("figures/4c.svg",height=6,width=4)
p2
dev.off()
library(cowplot)
# Setup the axes
all_gene_list_overtime$time_point[all_gene_list_overtime$time_point == 1] = "Mid-log"
all_gene_list_overtime$time_point[all_gene_list_overtime$time_point == 2] = "0.5 hour"
all_gene_list_overtime$time_point[all_gene_list_overtime$time_point == 3] = "1 hour"
all_gene_list_overtime$time_point[all_gene_list_overtime$time_point == 4] = "2 hours"
all_gene_list_overtime$time_point[all_gene_list_overtime$time_point == 5] = "4 hours"
all_gene_list_overtime$time_point[all_gene_list_overtime$time_point == 6] = "5.5 hours"
all_gene_list_overtime$time_point= factor(all_gene_list_overtime$time_point,levels=c("Mid-log","0.5 hour","1 hour","2 hours","4 hours","5.5 hours"))
all_gene_list_overtime$glucose = ifelse(all_gene_list_overtime$time_point == "Mid-log","Glucose","Galactose")
all_gene_list_overtime$glucose = factor(all_gene_list_overtime$glucose,levels=c("Glucose","Galactose"))

write.csv(all_gene_list_overtime,"figures/s5b.csv")
write.csv(all_gene_list, "figures/s5d.csv")

#all_gene_list_overtime$by[grepl("Y",all_gene_list_overtime$target_id) & !grepl("cbs",all_gene_list_overtime$target_id)] = "BY"
all_gene_list_overtime$by = "CBS"
all_gene_list_overtime$by[grepl("Y",all_gene_list_overtime$target_id) & !grepl("cbs",all_gene_list_overtime$target_id)] = "BY"
#all_gene_list_overtime$by = "CBS"

#all_gene_list_overtime$se = sqrt(all_gene_list_overtime$tpm)
all_gene_list_overtime$se = (sqrt(all_gene_list_overtime$tpm)/all_gene_list_overtime$tpm)
p1=all_gene_list_overtime %>% filter( all_gene_list_overtime$target_id == "YBR018C" | all_gene_list_overtime$target_id == "YBR018C_cbs_cbs" ) %>%  mutate(allele = ifelse(by == "CBS","CBS2888","BY"))  %>% ggplot(aes(y=log2(tpm +1),x=as.factor(time_point),color=allele,group=allele)) + geom_line() + geom_point() + xlab("Time point") + ylab("Log2(TPM+1)") + theme_bw() +expand_limits(y=0) + ggtitle("GAL1") + facet_grid(~ glucose,scales = "free",space="free_x")   + scale_color_brewer(name="Strain",palette = "Set1")+  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p2=all_gene_list_overtime %>% filter( all_gene_list_overtime$target_id == "YBR019C" | all_gene_list_overtime$target_id == "YBR019C_cbs_cbs" ) %>% mutate(allele = ifelse(by == "CBS","CBS2888","BY"))  %>% ggplot(aes(y=log2(tpm + 1),x=time_point,color=allele, group=allele)) + geom_line() + geom_point() + xlab("Time point") + ylab("Log2(TPM+1)") + theme_bw() + expand_limits(y=0) + ggtitle("GAL10") + scale_color_brewer(name="Strain",palette = "Set1")  + facet_grid(~glucose,scales = "free",space="free_x") +  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p3=all_gene_list_overtime %>% filter( all_gene_list_overtime$target_id == "YBR020W" | all_gene_list_overtime$target_id == "YBR020W_cbs_cbs" ) %>% mutate(allele = ifelse(by == "CBS","CBS2888","BY"))  %>% ggplot(aes(y=log2(tpm + 1),x=time_point,color=allele, group=allele)) + geom_line() + geom_point() + xlab("Time point") + ylab("Log2(TPM+1)") + theme_bw() + expand_limits(y=0) + ggtitle("GAL1")+ scale_color_viridis_d() +  facet_grid(~glucose,scales = "free",space="free_x") + scale_color_brewer(name="Strain",palette = "Set1") +  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p4=all_gene_list_overtime %>% filter( all_gene_list_overtime$target_id == "gal2a_cbs" | all_gene_list_overtime$target_id == "YLR081W" )%>% mutate(allele = ifelse(by == "CBS","CBS2888","BY")) %>% ggplot(aes(y=log2(tpm + 1),x=time_point,color=allele, group=allele)) + geom_line() + geom_point() + xlab("Time point") + ylab("Log2(TPM+1)") + theme_bw() + expand_limits(y=0) + ggtitle("GAL2a")+ scale_color_viridis_d() +  facet_grid(~glucose,scales = "free",space="free_x") + scale_color_brewer(name="Strain",palette = "Set1") +  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p5=all_gene_list_overtime %>% filter( all_gene_list_overtime$target_id == "gal2b_cbs" | all_gene_list_overtime$target_id == "YLR081W" )%>%mutate(allele = ifelse(by == "CBS","CBS2888","BY"))   %>% ggplot(aes(y=log2(tpm + 1),x=time_point,color=allele, group=by)) + geom_line() + geom_point() + xlab("Time point") + ylab("Log2(TPM+1)") + theme_bw() + expand_limits(y=0) + ggtitle("Gal2b")+ scale_color_viridis_d() + facet_grid(~glucose,scales = "free",space="free_x") + scale_color_brewer(name="Strain",palette = "Set1") +  theme(axis.text.x = element_text(angle = 90, hjust = 1))
#p6=all_gene_list_overtime %>% filter( all_gene_list_overtime$target_id == "gal2_cbs" | all_gene_list_overtime$target_id == "YLR081W" )%>% ggplot(aes(y=log2(tpm + 1),x=time_point,color=by)) + geom_line() + geom_point() + xlab("Time point") + ylab("Log2(TPM+1)") + theme_bw() + expand_limits(y=0) + ggtitle("GAL2")+ scale_color_viridis_d() + geom_errorbar(aes(ymin=(log2(tpm + 1 - se)),ymax=(log2(tpm + 1 + se))),color="black",alpha=0.5,width=0.5)


pgm1_time_points = melt(pgm1_time_points,id.vars = c("time_point"))


pgm1_time_points$time_point[pgm1_time_points$time_point == 1] = "Mid-log"
pgm1_time_points$time_point[pgm1_time_points$time_point == 2] = "0.5 hour"
pgm1_time_points$time_point[pgm1_time_points$time_point == 3] = "1 hour"
pgm1_time_points$time_point[pgm1_time_points$time_point == 4] = "2 hours"
pgm1_time_points$time_point[pgm1_time_points$time_point == 5] = "4 hours"
pgm1_time_points$time_point[pgm1_time_points$time_point == 6] = "5.5 hours"
pgm1_time_points$time_point= factor(pgm1_time_points$time_point,levels=c("Mid-log","0.5 hour","1 hour","2 hours","4 hours","5.5 hours"))
pgm1_time_points$glucose = ifelse(pgm1_time_points$time_point == "Mid-log","Glucose","Galactose")
pgm1_time_points$glucose = factor(pgm1_time_points$glucose,levels=c("Glucose","Galactose"))
pgm1_time_points$by = "DIV"
pgm1_time_points$by[grepl("ref",pgm1_time_points$variable)] = "REF"
#pgm1_time_points$by = "CBS"

#pgm1_time_points$se = sqrt(pgm1_time_points$tpm)
pgm1_time_points$se = (sqrt(pgm1_time_points$value)/pgm1_time_points$value)

p6 = pgm1_time_points %>% mutate(allele = ifelse(by == "DIV","CBS2888","BY")) %>% ggplot(aes(y=log2(value + 1),x=time_point,color=allele, group=by)) + geom_line() + geom_point(size=4) + xlab("Time point") + ylab("Log2(allele count)") +
  theme_bw() + expand_limits(y=0) + theme(text=element_text(size=20)) + ggtitle("PGM1")+ scale_color_viridis_d() +
  geom_errorbar(aes(ymin=(pmax(log2(value + 1) - 2 * se,0)),ymax=(log2(value + 1) + 2 *se)),color="black",alpha=0.5,width=0.5)+ facet_grid(~glucose,scales = "free",space="free_x") + scale_color_brewer(name="Allele",palette = "Set1") +  theme(axis.text.x = element_text(angle = 90, hjust = 1))


p6 = pgm1_time_points %>% mutate(allele = ifelse(by == "DIV","CBS2888","BY")) %>% ggplot(aes(y=log2(value + 1),x=time_point,color=allele, group=by)) + geom_line() + geom_point() + xlab("Time point") + ylab("Log2(allele count)") +
  theme_bw() + expand_limits(y=0)  + ggtitle("PGM1")+ scale_color_viridis_d() +
  facet_grid(~glucose,scales = "free",space="free_x") + scale_color_brewer(name="Allele",palette = "Set1") +  theme(axis.text.x = element_text(angle = 90, hjust = 1))
svg("figures/S6.svg",width=8, height=6)
p6
dev.off()

svg("figures/S23.svg",width=12,height=8)
plot_grid(p1,p2,p3,p4,p5,p6, nrow=3)
dev.off()
