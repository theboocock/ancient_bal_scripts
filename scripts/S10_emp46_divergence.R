### Figure 2B
### emp36 analysis
library(tidyverse)
emp46 = read.table("data/balancing_selection/emp46/df_emp46.txt", header=F)
emp46i = read.table("data/balancing_selection//emp46/df_emp_indels.txt", header=F)
emp46$window = seq(1,1142,by=100) + 150
emp46$indels = emp46i$V1
emp_p = pivot_longer(emp46,-window) 
svg("figures/S10.svg",width=4,height = 3)
emp_p %>% ggplot(aes(y=100 - value * 100,x=window,color=name)) + geom_point() + ylab("Sequence dissimilarity (%)") + xlab("Position")+ theme_bw() + xlim(c(0,1342))+
  scale_color_brewer(palette = "Set1",name="Indel treatment", labels=c("Mismatch","Missing"))
dev.off()

