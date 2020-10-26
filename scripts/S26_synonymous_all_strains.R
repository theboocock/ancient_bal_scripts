####
#### @author James Boocock
#### @date 24 Sept 2019
####
####


strain_annot_gal = read.csv("data/balancing_selection//suptable_annotations2.csv", header=T)

dn_ds_f = list.files("data/balancing_selection//gall_dn_ds_flanks_all/",pattern="*dn_ds", full.names=T)
gal_genes = c("YBR018C","YBR019C","YBR020W","YLR081W")
sample_names_c = c()
gal_genes_list = list()
i = 1
for(file in dn_ds_f){
  print(file)
  a = read.table(file, header=F,stringsAsFactors = F)
  #print(head(a))
  sample_names = str_split(str_replace_all(basename(file),"_[0-9]",""),"\\.")[[1]][1]
  sample_names_c = c(sample_names_c, sample_names)
  overall = a %>% filter(V1 %in% gal_genes)
  annot_gal = strain_annot_gal %>% filter(strain_annot_gal$Standardized.name ==sample_names)
  if(nrow(overall) == 0){
    print(sample_names)
    next
  }
  overall$sample_name = sample_names
  overall$order = 1:nrow(overall)
  
  gal_genes_df = merge(overall, annot_gal, by.x="sample_name",by.y=1)
  gal_genes_df = gal_genes_df[order(gal_genes_df$order),]
  #if(sum(g))
  #gal_genes_df = data.frame()
  gal_genes_list[[i]] = gal_genes_df
  i =i + 1
}

gal_gene_all = bind_rows(gal_genes_list)

gal_gene_all %>% filter(V2 == "OVERALL" & V1 == "YBR018C") %>% group_by(V1, GAL7 ) %>% 
  summarise(mean=mean(V4),lo=quantile(V4,probs=c(0.025),hi=quantile(V4,probs=0.975)),min=min(V4),max=max(V4), n = length(V4))
gal_gene_all %>% filter(V2 == "OVERALL" & V1 == "YBR019C") %>% group_by(V1, GAL10 ) %>% 
  summarise(mean=mean(V4),lo=quantile(V4,probs=c(0.025),hi=quantile(V4,probs=0.975)),min=min(V4),max=max(V4), n = length(V4))
gal_gene_all %>% filter(V2 == "OVERALL" & V1 == "YBR020W") %>% group_by(V1, GAL1 ) %>% 
  summarise(mean=mean(V4),lo=quantile(V4,probs=c(0.025),hi=quantile(V4,probs=0.975)),min=min(V4),max=max(V4), n= length(V4))



p1 =gal_gene_all %>% filter(V2 == "OVERALL"  & V1 == "YBR020W") %>% ggplot(aes(y=V4,x=GAL1)) + geom_jitter() + ylab("Synonymous substitutions per site") + xlab("Allele")  + theme_bw()+ ggtitle("GAL1") + theme(text=element_text(size=20))
p2 =gal_gene_all %>% filter(V2 == "OVERALL"  & V1 == "YBR019C") %>% ggplot(aes(y=V4,x=GAL10)) + geom_jitter() + ylab("Synonymous substitutions per site") + xlab("Allele")+ theme_bw() +ggtitle("GAL10")   + theme(text=element_text(size=20))
p3 =gal_gene_all %>% filter(V2 == "OVERALL"  & V1 == "YBR018C") %>% ggplot(aes(y=V4,x=GAL7)) + geom_jitter() + ylab("Synonymous substitutions per site") + xlab("Allele") +theme_bw() + ggtitle("GAL7") + theme(text=element_text(size=20))

svg("figures/S26.svg")
cowplot::plot_grid(p1,p2,p3,labels=c("a)","b)","c)"),ncol=3,label_size = 20)
dev.off()
p4 =gal_gene_all %>% filter(V2 == "OVERALL"  & V1 == "YBR020W") %>% ggplot(aes(y=V3,x=GAL1)) + geom_jitter() + ylab("Synonymous substitutions per site") + xlab("Allele")  + theme_bw()+ ggtitle("GAL1")
p5 =gal_gene_all %>% filter(V2 == "OVERALL"  & V1 == "YBR019C") %>% ggplot(aes(y=V3,x=GAL10)) + geom_jitter() + ylab("Synonymous substitutions per site") + xlab("Allele")+ theme_bw() +ggtitle("GAL10")
p6 =gal_gene_all %>% filter(V2 == "OVERALL"  & V1 == "YBR018C") %>% ggplot(aes(y=V3,x=GAL7)) + geom_jitter() + ylab("Synonymous substitutions per site") + xlab("Allele") +theme_bw() + ggtitle("GAL7")
cowplot::plot_grid(p4,p5,p6,labels=c("a)","b)","c)"),ncol=3)
#gal_gene_all %>% filter(V2 == "OVERALL" & V1 == "YLR081W") %>% group_by(V1, GAL2 ) %>% summarise(mean=mean(V4),lo=quantile(V4,probs=c(0.025),hi=quantile(V4,probs=0.975)),min=min(V4),max=max(V4))
#gal_gene_all %>% filter(V2 == "OVERALL" & V1 == "YBR019C") %>% group_by(V1, GAL10 ) %>% summarise(mean=mean(V4),lo=quantile(V4,probs=c(0.025),hi=quantile(V4,probs=0.975)),min=min(V4),max=max(V4))
#gal_gene_all %>% filter(V2 == "OVERALL" & V1 == "YBR020W") %>% group_by(V1, GAL1 ) %>% summarise(mean=mean(V4),lo=quantile(V4,probs=c(0.025),hi=quantile(V4,probs=0.975)),min=min(V4),max=max(V4))



