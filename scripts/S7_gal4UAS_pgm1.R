source("scripts//plate_read_utils.R")
# Figure S7
plate = "data/plate_reader/pgm1_mutant_growth/growth_pinning_glu_gal_2nd.csv"
conditions = "data/plate_reader/pgm1_mutant_growth/pgm1_mutant_conditions.csv"
samples = "data/plate_reader/pgm1_mutant_growth/pgm1_mutant_samples.csv"
mutant_pheno=  "data/plate_reader/pgm1_mutant_growth/pgm1_mutant_annotation.csv"

galactose_growth_binding_site_ko = process_phenotype_mutant_plates(plate, conditions, samples, mutant_pheno)
#galactose_growth_binding_site_ko$df_long_m
svg("figures/S7.svg")
galactose_growth_binding_site_ko$df_m_all %>% filter(cond == "GAL" & well != "E11" & epistatic != "DIV / DIV_MUT / DIV") %>% ggplot(aes(y=(doubling),x=epistatic))  +
  geom_point(size=4,position = position_jitter(width=0.2))  + theme_bw() + theme(axis.text = element_text(size=22))  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
