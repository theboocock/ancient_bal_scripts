## Load all the information about each of the crosses.
## Contains all the information for figure 1 incuding the supplementary information. 
### Makes figures 1a and 1c, S1,S2, and S5.
library(data.table)
#library()
library(ggplot2)
library(dplyr)
library(reshape2)
source("scripts/utils.R")
gal_info_list = readRDS("data/gal_pheno_list.RDs")
load('data/galData_segs.RData')
load('data/jointPeakEffects_JS_variants.RData')
### Phenotypes from the 3043 cross. CBS and CLIB.
pheno_3043 = gal_info_list$pheno_3043
c_y = gal_info_list$segs_3043
gal_markers_3043 = gal_info_list$gal_markers_3043
peak_markers = jointPeakEffects$`Galactose;;1`$`3043`
y = scale(pheno_3043[,grep("Gal",colnames(pheno_3043))])
genotype_names = unlist(lapply(strsplit(rownames(c_y),"_"), function(x) { paste(x[1],x[2],x[3],x[4],sep="_")}))
genotypes_idx_3043 = which(genotype_names %in% gal_markers_3043)
markers_residuals_3043 = peak_markers$peaks[!(peak_markers$peaks %in% gal_markers_3043)]
gt_norm = cor(c_y)
gt = apply(c_y[which(genotype_names %in% markers_residuals_3043),],1,scale)
y_resid = residuals(lm((y ~ gt)))
c_y_out_unscaled = c_y[genotypes_idx_3043,]
c_y_out = data.frame((apply(c_y_out_unscaled,1,scale)))
colnames(c_y_out) = c(gal_markers_3043)
chisq.test(table(paste(c_y_out[,1],c_y_out[,2], c_y_out[,3])))
m2= (lm(y ~ c_y_out[,1]+c_y_out[,2]+c_y_out[,3]))
m1= (lm(y ~ -1 +  c_y_out[,1]*c_y_out[,2]*c_y_out[,3]))
summary((lm(y_resid ~ chrII_238609_C_A*chrXI_201367_T_C*chrXII_299925_A_C, data=c_y_out)))
#### 
colnames(galData$data_3043$geno)
genotypes = t(c_y_out_unscaled)
data_3043 = cbind(y, genotypes)
data_3043 = as.data.frame(data_3043)
### Get lod scores ### 

#cor(pheno_3004,(c_y_3004))
##### #####
# 0 is CBS for the first SNP.
# 1 is CBS for the second SNP.
# 2 is CBS for the 3rd SNP.
data_3043[,2][data_3043[,2] == 0] = 2
data_3043[,2][data_3043[,2] == 1] = 0
data_3043[,2][data_3043[,2] == 2] = 1
summary(lm(scale(data_3043[,1]) ~ scale(data_3043[,2])*scale(data_3043[,3])*scale(data_3043[,4])))
colnames(data_3043) = c("Pheno","GAL1/10/7","PGM1","GAL2")
data_3043$epistatic = as.character(apply(data_3043, 1,function(x){x[x==0] = "REF"; x[x==1] = "ALT"; paste0(x[2:4],collapse = " / ")}))
p_3043_unresidualized =data_3043 %>% ggplot(aes(y=Pheno,x=epistatic)) + theme_bw() + geom_boxplot() + ylab("Normalized Phenotype") +
  xlab("Genotypes at the GAL1/10/7, PGM1, and GAL2 loci") + theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
  theme(axis.text = element_text(size=15)) + theme(axis.title= element_text(size=15))
data_3043$pheno_normalized = y_resid
p_3043_resid = data_3043 %>% ggplot(aes(y=pheno_normalized,x=epistatic)) + theme_bw() + geom_boxplot() + ylab("Normalized Phenotype") +
  xlab("Genotypes at the GAL1/10/7, PGM1, and GAL2 loci") + theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
  theme(axis.text = element_text(size=15)) + theme(axis.title= element_text(size=15))
pheno_3004 = gal_info_list$pheno_3004
c_y_3004 = gal_info_list$segs_3004
gal_markers_3004 = gal_info_list$gal_markers_3004
cor_markers  = !duplicated(c_y_3004,MARGIN=1)
c_y_subset=c_y_3004[cor_markers,]
c_y_subset = apply(c_y_subset,1, function(x){scale(x)})

genotype_names = unlist(lapply(strsplit(rownames(c_y_3004),"_"), function(x) { paste(x[1],x[2],x[3],x[4],sep="_")}))
chrom = unlist(lapply(strsplit(colnames(c_y_subset),"_"), function(x){x[1]}))
pos = as.numeric(unlist(lapply(strsplit(colnames(c_y_subset),"_"), function(x){x[2]})))
y = scale(pheno_3004[,grep("Gal",colnames(pheno_3004))])
r = cor(y,c_y_subset)
lods = get_lod_from_r(r[1,],n=length(y))


mapping_df = data.frame(chrom=chrom,pos=pos,lod=lods)  

mapping_df %>% ggplot(aes(y=lod,x=pos)) +  facet_wrap(~chrom,scales="free_x",nrow=1) + geom_point() + theme_bw() + xlab("Position") + ylab("LOD Score") + theme(text=element_text(size=20),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + geom_hline(yintercept = 4,color="red")


c_y_subset_all = c_y_3004[which(genotype_names %in% gal_markers_3004),]
chrom = unlist(lapply(strsplit((gal_markers_3004),"_"), function(x){x[1]}))
pos = as.numeric(unlist(lapply(strsplit((gal_markers_3004),"_"), function(x){x[2]})))
qtl_plots = data.frame(y=y, qtl1=c_y_subset_all[1,], qtl2=c_y_subset_all[2,],qtl3=c_y_subset_all[3,])
qtl_plots$int = paste0(qtl_plots[,2],qtl_plots[,3],qtl_plots[,4])
colnames(qtl_plots) = c("y","chrII","chrXI","chrXII")
qtl_plots = melt(qtl_plots,id.vars = c("y"))
p_bot = qtl_plots %>% ggplot(aes(y=y,x=factor(value))) + geom_boxplot() + facet_wrap(~variable) + ylab("Normalized growth") + xlab("Genotype calls")
y_unscaled = pheno_3004[,grep("Gal",colnames(pheno_3004))]
genotype_names = unlist(lapply(strsplit(rownames(c_y_3004),"_"), function(x) { paste(x[1],x[2],x[3],x[4],sep="_")}))
genotypes_idx = which(genotype_names %in% gal_markers_3004)
c_y_out_3004_unormalized = c_y_3004[genotypes_idx,]
c_y_out_3004 = t(apply(c_y_out_3004_unormalized,1,scale))
peak_markers = jointPeakEffects$`Galactose;;1`$`3004`
markers_residuals_3004 = peak_markers$peaks[!(peak_markers$peaks %in% gal_markers_3004)]
gt = apply(c_y_3004[which(genotype_names %in% markers_residuals_3004),],1,scale)
y_resid = residuals(lm((y ~ gt)))
### Test that the loci are unlinked.
chisq.test(table(paste(c_y_out_3004[1,],c_y_out_3004[2,], c_y_out_3004[3,])))

m2 = (lm(y~ -1 + c_y_out_3004[1,]*c_y_out_3004[2,]*c_y_out_3004[3,]))
m1 = (lm(y~ -1 + c_y_out_3004[1,]+c_y_out_3004[2,]+c_y_out_3004[3,]))
m3 = (lm(y~ -1 + c_y_out_3004[1,]+c_y_out_3004[2,]+c_y_out_3004[3,] +  c_y_out_3004[1,]:c_y_out_3004[2,]+  c_y_out_3004[1,]:c_y_out_3004[3,]  + c_y_out_3004[2,]:c_y_out_3004[3,]))
# Two way interaction model
summary(m3)
# Three-way interaction model
summary(m2)
# Additive model
summary(m1)
genotypes = t(c_y_out_3004_unormalized)
data_3004 = cbind(y, genotypes)
colnames(data_3004) = c("Pheno","GAL1/10/7","PGM1","GAL2")
data_3004 = as.data.frame(data_3004)
data_3004$epistatic = as.character(apply(data_3004, 1,function(x){ x[x==0] = "REF"; x[x==1] = "ALT"; paste0(x[2:4],collapse = " / ")}))
### Figure S1
y_3004 = scale(pheno_3004[,grep("Gal",colnames(pheno_3004))])
cor_3004 =cor(y_3004,t(gal_info_list$segs_3004))
library(stringr)
splitstr = colnames(cor_3004)
chrom = unlist(lapply(str_split(splitstr,"_"), function(x){x[1]}))
pos =  as.numeric(unlist(lapply(str_split(splitstr,"_"), function(x){x[2]})))
r_df =get_lod_from_r(cor_3004,n = 867)
df_all= data.frame(chrom=chrom,pos=pos,r=r_df[1,])
df_all$cross = "YJM981xCBS2888"
#df_all %>% filter(chrom ==) %>% ggplot(aes(y=r,x=pos))  + geom_point()+ facet_wrap(~chrom)
y_3043 = scale(pheno_3043[,grep("Gal",colnames(pheno_3043))])
cor_3043 = cor(y_3043,t(gal_info_list$segs_3043))
splitstr = colnames(cor_3043)
chrom = unlist(lapply(str_split(splitstr,"_"), function(x){x[1]}))
pos =  as.numeric(unlist(lapply(str_split(splitstr,"_"), function(x){x[2]})))
r_df = get_lod_from_r(cor_3043,n=943)
df_two = data.frame(chrom=chrom,pos=pos,r=r_df[1,])
df_two$cross = "CBS2888xCLIB219"
df_all= bind_rows(df_all,df_two)
data_in = data.frame(chrom=c("chrII","chrXI","chrXII"),start=c(273703,201771,289250),end=c(281443,205253,292394))
p1 = df_all %>% filter(cross == "YJM981xCBS2888") %>% ggplot(aes(y=r,x=pos)) + geom_point()  + facet_wrap(~chrom,nrow=2,scales="free_x") + theme_bw() + ylim(c(0,45))    + theme(text=element_text(size = 16)) + ylab("LOD") + xlab("Position")
p1 = p1 + geom_vline(data = data_in, mapping = aes(xintercept=start), color="red") + geom_vline(data = data_in, mapping = aes(xintercept=end),color="red") + theme(axis.text.x=element_text(angle=90, hjust=1)) + ggtitle("YJM981xCBS2888")
p2 = df_all %>% filter(cross == "CBS2888xCLIB219") %>% ggplot(aes(y=r,x=pos)) + geom_point()  + facet_wrap(~chrom,nrow=2,scales="free_x") + theme_bw() + ylim(c(0,45))    + theme(text=element_text(size = 16)) + ylab("LOD") + xlab("Position")
p2 = p2 + geom_vline(data = data_in, mapping = aes(xintercept=start), color="red") + geom_vline(data = data_in, mapping = aes(xintercept=end),color="red") + theme(axis.text.x=element_text(angle=90, hjust=1)) + ggtitle("CBS2888xCLIB219")
p5= df_all %>% filter(cross == "YJM981xCBS2888")  %>% filter(chrom == "chrII" | chrom == "chrXI" | chrom == "chrXII") %>% ggplot(aes(y=r,x=pos)) + geom_point()  + facet_wrap(~chrom,nrow=1,scales="free")  + theme(text=element_text(size = 16))+ theme_bw() + ylab("LOD") + xlab("Position")
p5 = p5+ geom_vline(data = data_in, mapping = aes(xintercept=start), color="red") + geom_vline(data = data_in, mapping = aes(xintercept=end),color="red")  + theme(axis.text.x=element_text(angle=90, hjust=1)) + ggtitle("YJM981xCBS2888")
p6 = df_all %>% filter(cross == "CBS2888xCLIB219")  %>% filter(chrom == "chrII" | chrom == "chrXI" | chrom == "chrXII") %>% ggplot(aes(y=r,x=pos)) + geom_point()  + facet_wrap(~chrom,nrow=1,scales="free")  + theme(text=element_text(size = 16))+ theme_bw()  + ylab("LOD") + xlab("Position")
p6 = p6+ geom_vline(data = data_in, mapping = aes(xintercept=start), color="red") + geom_vline(data = data_in, mapping = aes(xintercept=end),color="red")  + theme(axis.text.x=element_text(angle=90, hjust=1))+ ggtitle("CBS2888xCLIB219")
png("figures/S1.png",width=30*100, height=15*100,res=150)
cowplot::plot_grid(p1,p5,p2,p6,labels=c("A","C","B","D"))
dev.off()
### Figure S2 ### 
p_3004 =data_3004 %>% ggplot(aes(y=Pheno,x=epistatic)) + theme_bw() + geom_boxplot(width=0.3,size=0.5,outlier.shape = NA)  + geom_point(size=0.5,alpha=0.5, position = position_jitter(w=0.2,h=0))+ ylab("Normalized Phenotype") + ylab("Normalized Phenotype") +
  xlab("Genotypes at the GAL1/10/7, PGM1, and GAL2 loci") + theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
  theme(axis.text = element_text(size=15)) + theme(axis.title= element_text(size=15)) + theme(axis.title.x=element_blank(),
                                                                                              axis.text.x=element_blank(),
                                                                                              axis.ticks.x=element_blank())
data_3004$pheno_normalized = scale(y_resid)
p_3004_resid =data_3004 %>% ggplot(aes(y=pheno_normalized,x=epistatic)) + theme_bw() + geom_boxplot()  +  ylab("Normalized Phenotype") +
  xlab("Genotypes at the GAL1/10/7, PGM1, and GAL2 loci") + theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
  theme(axis.text = element_text(size=15)) + theme(axis.title= element_text(size=15))
data_3043$cross = "CLIB419xCBS2888"
data_3004$cross = "YJM981xCBS2888"
merged_crosses = rbind(data_3043,data_3004)
merged_crosses_ggplot2 = reshape2::melt(merged_crosses,measure.vars=c("Pheno","pheno_normalized"),id.vars=c("cross","epistatic"))
levels(merged_crosses_ggplot2$variable) = c("Scaled phenotype","Scaled phenotype (Other QTL regressed out)")
p1 = merged_crosses_ggplot2 %>% filter(merged_crosses_ggplot2$variable != "Scaled phenotype (Other QTL regressed out)")%>% ggplot(aes(y=value,x=epistatic)) + theme_bw() + geom_boxplot(width=0.3,size=0.5,outlier.shape = NA) + ylab("Growth (arbitrary units)") +
  xlab("Genotypes at the GAL1/10/7, PGM1, and GAL2 loci") + theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
  theme(axis.text = element_text(size=15)) + theme(axis.title= element_text(size=15)) + facet_wrap(~cross ,nrow=2)
svg("figures/S2.svg", width=12, height=8)
p1
dev.off()
### Figure 1 A and C ###
source("scripts//plate_read_utils.R")
library(platetools)
library(segmented)
# Load the allele replacement experiments. 
plate = "data/plate_reader/mutants_final_freezer//gal_mutants_gal_feb2019.csv"
conditions = "data/plate_reader/mutants_final_freezer/gal_mutants_glu_feb2019_conditions.csv"
samples = "data/plate_reader/mutants_final_freezer/gal_mutants_glu_feb2019_samples.csv"
mutant_pheno=  "data/plate_reader/mutants_final_freezer/gal_mutant_pheno.csv"
galactose_growth = process_phenotype_mutant_plates(plate, conditions, samples, mutant_pheno )
make_genotype_boxplot(galactose_growth$df_m%>% filter(ID <= 16),variable="doubling")
parents = galactose_growth$df_m%>% filter(ID <= 16)
# Rename to a numeric vector
parents$epistatic_rename = str_replace_all(str_replace_all(parents$epistatic,"C2888","ALT"),"WT","REF")
parents$epistatic_numeric= as.numeric(factor((parents$epistatic_rename)))
## Ok this works finally I should label the biological replicates
p_validation=  parents %>% ggplot(aes(y=doubling,x=epistatic_numeric, group=epistatic_numeric)) + theme_bw() +   
  geom_point(aes(x=epistatic_numeric), size=4,position = position_jitter(width=0.2))+ ylab("Growth (doublings per hour)") + ylim(c(0.2,0.41)) + 
  xlab("Genotypes at the GAL1/10/7, PGM1, and GAL2 loci") + theme(plot.title = element_text(hjust = 0.5,size=22),axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text = element_text(size=22)) + theme(axis.title= element_text(size=22)) +ggtitle("Allele replacements") +
  scale_x_continuous(breaks=1:8,labels=levels(factor(parents$epistatic_rename))) + scale_color_manual(values=c("#984ea3","#ff7f00")) + theme(legend.position = c(.9,0.2),legend.background = element_rect(color="white"),legend.box.background = element_rect(size=1,color="black"))
p_3004_no_x = data_3004 %>% ggplot(aes(y=Pheno,x=epistatic)) + theme_bw() + geom_boxplot(width=0.5,size=0.5)  +  ylab("Normalized growth rate") +
  xlab("Genotypes at the GAL1/10/7, PGM1, and GAL2 loci") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text = element_text(size=22))  + ggtitle("QTL mapping")+ theme(plot.title = element_text(hjust = 0.5,size=22), axis.title= element_text(size=22)) +   theme(axis.title.x=element_blank(),
                                                                                                                                                                               axis.text.x=element_blank(),
                                                                                                                                                                               axis.ticks.x=element_blank())
p_validation =  parents %>% ggplot(aes(y=doubling,x=epistatic_rename, group=epistatic_rename)) + theme_bw() +   
  geom_point(size=4,position = position_jitter(width=0.2))+ ylab("Growth (doublings per hour)") + ylim(c(0.2,0.41)) + 
  xlab("Genotypes at the GAL1/10/7, PGM1, and GAL2 loci") + theme(plot.title = element_text(hjust = 0.5,size=22),axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text = element_text(size=22)) + theme(axis.title= element_text(size=22)) +ggtitle("Allele replacements") +
  scale_color_manual(values=c("#984ea3","#ff7f00")) + theme(legend.position = c(.9,0.2),legend.background = element_rect(color="white"),legend.box.background = element_rect(size=1,color="black"))

pdf("figures/1.pdf",width = 16, height = 10)
cowplot::plot_grid(p_3004_no_x + theme(legend.position = "none"),p_validation,labels="AUTO",ncol=1, rel_heights = c(1,1.5), align = "v")
dev.off()

parents$GAL7 = 1 -(as.numeric(factor(parents$GAL7))-1)
parents$PGM1 = 1 - (as.numeric(factor(parents$PGM1)) -1)
parents$GAL2 = 1 - (as.numeric(factor(parents$GAL2))-1)

paste(parents$GAL7, parents$PGM1, parents$GAL2)
# Figure S4
gal_strains_lm = lm(scale(parents$doubling) ~  (parents$GAL7 *  parents$PGM1  * parents$GAL2))
gal_strains_lm_two_way = lm(scale(parents$doubling) ~  (parents$GAL7  +  parents$PGM1   + parents$GAL2) + parents$GAL7:parents$PGM1 + parents$GAL7:parents$GAL2 + parents$PGM1:parents$GAL2)
gal_strains_lm_add = lm(scale(parents$doubling) ~  (parents$GAL7  +  parents$PGM1   + parents$GAL2))
parents$epistatic_rename = str_replace_all(str_replace_all(parents$epistatic,"C2888","ALT"),"WT","REF")
# Compare different strains with the galactose pathway.
p1 = parents %>% filter(epistatic_rename == "REF / REF / REF" | epistatic_rename == "ALT / ALT / ALT")
t.test(p1$doubling ~ p1$epistatic_rename)
p1 = parents %>% filter(epistatic_rename == "ALT / REF / ALT" | epistatic_rename == "ALT / ALT / ALT")
t.test(p1$doubling ~ p1$epistatic_rename)
galactose_growth$df_long_m$epistatic_rename = str_replace_all(str_replace_all(galactose_growth$df_long_m$epistatic,"C2888","ALT"),"WT","REF")
svg("figures/S4.svg",width=8,height=6)
galactose_growth$df_long_m %>%  filter(newtime >= 10 & newtime <=50) %>% mutate(value2=value - min(value)) %>% ggplot(aes(y=(value2),x=newtime/4,color=epistatic_rename,group=well))  +
  geom_line()  + theme_bw()  + ylim(c(0,1.1)) +  xlim(c(10/4,50/4)) + scale_color_brewer(palette = "Set1") + xlab("Time in hours") + ylab("Growth (arbitrary units)")
dev.off()

