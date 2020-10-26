#### ###### 
library(tidyverse)
library(ggpubr)
source("scripts//plate_read_utils.R")

plate = "data/plate_reader//mutants_final_freezer/gal_mutants_one_feb2019.csv"
conditions = "data/plate_reader/mutants_final_freezer/gal_mutants_glu_feb2019_conditions.csv"
samples = "data/plate_reader//mutants_final_freezer/gal_mutants_glu_feb2019_samples.csv"
mutant_pheno=  "data/plate_reader//mutants_final_freezer/gal_mutant_pheno.csv"
mixed_growth = process_phenotype_mutant_plates(plate, conditions, samples, mutant_pheno = mutant_pheno)
mixed_growth$df_long_m$epistatic_rename = str_replace_all(str_replace_all(mixed_growth$df_long_m$epistatic,"C2888","DIV"),"WT","REF")


mixed_growth_ar = mixed_growth$df_long_m %>% filter(epistatic_rename == "REF / REF / REF" | epistatic_rename == "DIV / DIV / DIV" )%>% filter(ID <= 16)

mixed_growth_ar$parents = "Allele replacements"

#mixed_growth$df_parents = 
mixed_growth_parents = mixed_growth$df_parents
mixed_growth_parents$epistatic_rename = str_replace_all(str_replace_all(mixed_growth_parents$epistatic,"CBS","DIV"),"BY","REF")


mixed_growth_parents$parents = "Parental strains"


mixed_growth_df = rbind(mixed_growth_ar,mixed_growth_parents)
svg("figures/S22a.svg",width = 12, height = 10)
mixed_growth$df_long_m%>% filter(ID <= 16) %>% ggplot(aes(y=value, x=newtime/4)) +  geom_point(aes(color=epistatic_rename)) + xlim(c(0,75/4))  +
  theme_bw()  +theme(plot.title = element_text(hjust = 0.5,size=26),text=element_text(size=20), axis.title= element_text(size=26),strip.text = element_text(size=26))  + scale_color_brewer(palette = "Set1")+ facet_wrap(~epistatic_rename,ncol=4,nrow=2)  + ylab("OD600") + xlab("Hours") 
dev.off()



mixed_growth$df_m$epistatic_rename = str_replace_all(str_replace_all(mixed_growth$df_m$epistatic,"C2888","DIV"),"WT","REF")
svg("figures/S22b.svg",width=12,height=10)
mixed_growth$df_m %>% ggplot(aes(y=switching,x=epistatic_rename)) + theme_bw() +   
  geom_point(size=4,position = position_jitter(width=0.2))+ ylab("Growth (doublings per hour)") + 
  xlab("Genotypes at the GAL1/10/7, PGM1, and GAL2 loci") + theme(plot.title = element_text(hjust = 0.5,size=22),axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text = element_text(size=22)) + theme(axis.title= element_text(size=22)) +ggtitle("Allele replacements") +
  theme(legend.position = c(.9,0.2),legend.background = element_rect(color="white"),legend.box.background = element_rect(size=1,color="black"))
dev.off()