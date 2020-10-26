### script makes S11. 4a 4b, and 4d.
library(tidyverse)
library(ggpubr)
source("scripts//plate_read_utils.R")

dropbox_in= "~/Dropbox/"
library(glue)
plate = glue(dropbox_in,"/PHDTHESIS/projects/introgressed_gal_cbs288a/data/mutants_final_freezer_gal/gal_mutants_glu_feb2019.csv")
conditions = glue(dropbox_in,"/PHDTHESIS/projects/introgressed_gal_cbs288a/data/mutants_final_freezer_gal/gal_mutants_glu_feb2019_conditions.csv")
#conditions = "/media/theboocock/data/Dropbox/PHDTHESIS/projects/introgressed_gal_cbs288a/data/mutants_final_freezer_gal/gal_mutants_glu_feb2019_conditions.csv"
samples = glue(dropbox_in,"//PHDTHESIS/projects/introgressed_gal_cbs288a/data/mutants_final_freezer_gal/gal_mutants_glu_feb2019_samples.csv")
#samples = "/media/theboocock/data/Dropbox/PHDTHESIS/projects/introgressed_gal_cbs288a/data/mutants_final_freezer_gal/gal_mutants_glu_feb2019_samples.csv"
mutant_pheno = glue(dropbox_in, "//PHDTHESIS/projects/introgressed_gal_cbs288a/data/mutants_final_freezer_gal/gal_mutant_pheno.csv")
#mutant_pheno=  "/media/theboocock/data/Dropbox/PHDTHESIS/projects/introgressed_gal_cbs288a/data/mutants_final_freezer_gal/gal_mutant_pheno.csv"
glu_growth = process_phenotype_mutant_plates(plate, conditions, samples, mutant_pheno )
glu_growth$df_m$epistatic_rename = str_replace_all(str_replace_all(glu_growth$df_m$epistatic,"C2888","ALT"),"WT","REF")
glu_growth$df_m$epistatic_numeric= as.numeric(factor((glu_growth$df_m$epistatic_rename)))
plate = "data/plate_reader/mutants_final_freezer_gal/gal_mutants_gal_feb2019.csv"
conditions = "data/plate_reader/mutants_final_freezer_gal/gal_mutants_glu_feb2019_conditions.csv"
samples = "data/plate_reader/mutants_final_freezer_gal/gal_mutants_glu_feb2019_samples.csv"
mutant_pheno=  "data/plate_reader/mutants_final_freezer_gal/gal_mutant_pheno.csv"
galactose_growth = process_phenotype_mutant_plates(plate, conditions, samples, mutant_pheno )
galactose_growth$df_m$epistatic_rename = str_replace_all(str_replace_all(galactose_growth$df_m$epistatic,"C2888","ALT"),"WT","REF")
galactose_growth$df_m$epistatic_numeric= as.numeric(factor((galactose_growth$df_m$epistatic_rename)))

svg("figures/S11.svg")
glu_growth$df_m %>% ggplot(aes(y=doubling,x=epistatic_numeric, group=epistatic_numeric)) + theme_bw() +   
  geom_point(aes(x=epistatic_numeric), size=4,position = position_jitter(width=0.2))+ ylab("Growth (doublings per hour)")  + 
  xlab("Genotypes at the GAL1/10/7, PGM1, and GAL2 loci") + theme(plot.title = element_text(hjust = 0.5,size=22),axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text = element_text(size=22)) + theme(axis.title= element_text(size=22)) +ggtitle("Allele replacements") +
  scale_x_continuous(breaks=1:8,labels=levels(factor(parents$epistatic_rename))) + scale_color_manual(values=c("#984ea3","#ff7f00")) + theme(legend.position = c(.9,0.2),legend.background = element_rect(color="white"),legend.box.background = element_rect(size=1,color="black"))
dev.off()
#glu_growth$df_m$epistatic_rename == "ALT / ALT / ALT" | 
svg("figures/4a.svg")
galactose_growth$df_m %>% filter(epistatic_rename == "ALT / ALT / ALT" | epistatic_rename == "REF / REF / REF") %>% ggplot(aes(y=doubling,x=epistatic_rename, group=epistatic_numeric)) + theme_bw() +   
  geom_point(size=4,position = position_jitter(width=0.2)) + geom_boxplot(width=0.1,outlier.shape = NA)+ ylab("Growth (doublings per hour)")  + 
  xlab("Genotypes at the GAL1/10/7, PGM1, and GAL2 loci") + theme(plot.title = element_text(hjust = 0.5,size=22),axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text = element_text(size=22)) + theme(axis.title= element_text(size=22)) +ggtitle("Allele replacements")
dev.off()

svg("figures/4b.svg")
glu_growth$df_m %>% filter(epistatic_rename == "ALT / ALT / ALT" | epistatic_rename == "REF / REF / REF") %>% ggplot(aes(y=doubling,x=epistatic_rename, group=epistatic_numeric)) + theme_bw() +   
  geom_point(size=4,position = position_jitter(width=0.2)) + geom_boxplot(width=0.1,outlier.shape = NA)+ ylab("Growth (doublings per hour)")  + 
  xlab("Genotypes at the GAL1/10/7, PGM1, and GAL2 loci") + theme(plot.title = element_text(hjust = 0.5,size=22),axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text = element_text(size=22)) + theme(axis.title= element_text(size=22)) +ggtitle("Allele replacements")
dev.off()


plate = "data/plate_reader/mutants_final_freezer/gal_mutants_one_feb2019.csv"
conditions = "data/plate_reader/mutants_final_freezer/gal_mutants_glu_feb2019_conditions.csv"
samples = "data/plate_reader/mutants_final_freezer/gal_mutants_glu_feb2019_samples.csv"
mutant_pheno=  "data/plate_reader/mutants_final_freezer/gal_mutant_pheno.csv"
mixed_growth = process_phenotype_mutant_plates(plate, conditions, samples, mutant_pheno )
mixed_growth$df_long_m$epistatic_rename = str_replace_all(str_replace_all(mixed_growth$df_long_m$epistatic,"C2888","ALT"),"WT","REF")#mixed_growth$$epistatic_rename = str_replace_all(str_replace_all(galactose_growth$df_m$epistatic,"C2888","ALT"),"WT","REF")
#galactose_growth$df_m$epistatic_numeric= as.numeric(factor((galactose_growth$df_m$epistatic_rename)))
svg("figures/4d.svg")
mixed_growth$df_long_m %>% filter(epistatic == "WT / WT / WT" | epistatic == "C2888 / C2888 / C2888") %>%
  ggplot(aes(y=value,x=newtime/4,color=epistatic_rename)) + geom_point()  + xlim(c(0,75/4)) + scale_color_manual(values = c("#e41a1c","grey50")) + theme_bw()
dev.off()
