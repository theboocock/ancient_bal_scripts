### ###
library(tidyverse)
library(stringr)
in_factors = list.files("data/forward_sims/factor_test/",pattern="factor*",full.names=T)
site_frequency_list = list()
ds_list = list()
summarise_slopes_list = list()
summarise_mutation_list  = list()
i = 1
factor = unlist(lapply(str_split(unlist(lapply(str_split(basename(in_factors),"\\."), function(x){x[[1]]})),"_"), function(x){x[2]}))
for (file_in in in_factors){
  print(file_in)
  x = readRDS(file_in)
  site_frequency_list[[i]] = data.frame(freq=x$site_frequency_spectrum,factor=as.numeric(factor[i]))
  x$pi_between_df$factor = as.numeric(factor[i])
  ds_list[[i]]= x$pi_between_df
  summarise_mutation_list[[i]] = x$mutation_df
  i =  i+ 1
}
names(site_frequency_list) = basename(in_factors)
names(ds_list) = basename(in_factors)

p1 = bind_rows(site_frequency_list,.id = "name") %>% mutate(real_freq = freq * factor) %>% ggplot(aes(x=real_freq)) +  geom_histogram(aes(x=real_freq)) + facet_wrap(~factor) + theme_bw() + 
  theme(text=element_text(size=20)) + xlab("Number of derived alleles") + ylab("Count") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p2 = bind_rows(ds_list, .id="name") %>% ggplot(aes(y=pi_between,x=start)) + geom_point() + facet_wrap(~factor) + theme_bw() + theme(text=element_text(size=20))  + 
  ylab("Pi_between") + xlab("Chromosome position")
cowplot::plot_grid(p1,p2,labels = "AUTO",label_size = 24)

ggsave(filename = "figures/S16.png",width=16,height=12, dpi=150)
#bind_rows(site_frequency_list,.id = "name") %>% mutate(real_freq = freq * factor) %>% filter(factor == 50000)


