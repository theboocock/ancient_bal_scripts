rds_alleles = readRDS("/media/theboocock/data/Dropbox/PHDTHESIS/projects/gal_final_github/data/permutations_ld.RDs")
abline = 0.5916645
df = data.frame(x=rds_alleles$triple)
svg("figures/S8.svg",width=12,height=3)
df %>% ggplot(aes(x=df$x)) + geom_histogram() + scale_x_continuous(limits = c(0,0.7), oob = scales::squish) + geom_vline(xintercept = abline,color="red") +
  xlab("Background distribution of epsilon") + ylab("Count") +  theme(text = element_text(size=30)) + theme_bw()
dev.off()
