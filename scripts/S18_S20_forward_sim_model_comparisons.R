library(tidyverse)
library(abc)
source("scripts/S17_real_data_conversion.R")
### Date: 
source("summarise_results_and_plot_fx.R")

#sumstats$df %>% pivot_longer(cols=c("alt_advantage","ref_advantage","cloning_rate","mitotic_recomb_rate","pulse_percentage","pulse_generation"),names_to = c("parameter"),values_to = c("val")) %>%
#  ggplot(aes(x=val)) + geom_histogram()+ facet_wrap(~parameter*sim_type,scales="free",ncol=6)
#sumstats$df  %>% mutate(TOTAL_GENERATIONS/)

#sumstats$df = sumstats$df %>% mutate(generations_ago=(TOTAL_GENERATIONS-pulse_generation) * factor)


#ggsave("figures/S3.png",width=8,height=6)
sumstats$df %>% mutate(fixed=!is.na(sumstats$sumstats$X1)) %>% 
  filter(sim_type == "introgression")  %>%  mutate(fixed=ifelse(fixed,"All six alleles present","Atleast one allele lost"))%>%
  ggplot(aes(x=generations_ago)) + geom_histogram(binwidth = 1e6*2) + scale_x_reverse() + facet_wrap(~fixed,nrow=2) + theme_bw()+ theme(text=element_text(size=20)) + xlab("Generations ago") + ylab("Count")
ggsave("figures/S18.png",width=8,height=6)
sumstats$df %>% mutate(fixed=is.na(sumstats$sumstats$X1)) %>% filter(sim_type == "introgression") %>% filter(!fixed) %>% summarise(med=TOTAL_GENERATIONS*factor - mean(pulse_generation)*factor,
                                                                                                                                   low=TOTAL_GENERATIONS*factor -  quantile(pulse_generation,prob=0.025)*factor,
                                                                                                                                   high=TOTAL_GENERATIONS*factor - quantile(pulse_generation,prob=0.975)*factor)
params = c(epsilon, df_out$V1[!is.na(df_out$V1)])
### ###
ab = abc(c(epsilon, df_out$V1[!is.na(df_out$V1)]), sumstats$df, sumstats$sumstats[,c(1,12:ncol(sumstats$sumstats))],tol=0.01,method="rejection")
keep_names_df = read.table("keep_rows.txt")
keep_names = keep_names_df$name
postpr_list = list()
j = 1
sim_type_f = factor(sim_type[!is.na(sumstats$sumstats$X1)],levels=c("balancing_selection","introgression","bal_intro_recent","bal_intro_almost_recent","bal_intro_ancient"))
for(tol in c(0.01,seq(0.05,0.5,by=.05))){
  print(tol)
  #  a = postpr(aa$x[nrow(aa$x),],(sim_type[!is.na(sumstats$sumstats$X1)]),aa$x[1:(nrow(aa$x)-1),],tol=tol,method="rejection")
  
  # ab = postpr(c(epsilon,0.03897673, out_df$max_ds), sim_type, sumstats$sumstats[,c("eld","slopes","max1","max2","max3","max4","max5","max6")] ,tol=0.2,method="mnlogistic")
  
  a=postpr(c(epsilon,df_out$V1[!is.na(df_out$V1)]),
           (sim_type[rownames(sumstats$sumstats) %in% keep_names]),
           sumstats$sumstats[rownames(sumstats$sumstats) %in% keep_names,c(1,12:ncol(sumstats$sumstats))],tol=tol,method = "rejection")
  #s#a = postpr(c(epsilon,df_out$V1[!is.na(df_out$V1)]),(sim_type[!is.na(sumstats$sumstats$X1)]),sumstats$sumstats[!is.na(sumstats$sumstats$X1),c(1,12:ncol(sumstats$sumstats))],tol=tol,method = "rejection")
  postpr_list[[j]] = a
  j = j+1
}



#summary(a)
names(postpr_list) = c(0.01,seq(0.05,0.5,by=.05))


bfs = (lapply(postpr_list, function(x){c(summary(x)$BayesF[4,])}))
bfs = t(bind_rows(bfs))
p_intro = as.data.frame(bfs) %>% mutate(tol=c(0.01,seq(0.05,0.5,by=.05))) %>% 
  dplyr::rename(bal_intro_almost_recent=V1, bal_intro_recent=V3, ancient_balancing_selection=V4,bal_intro_ancient=V2,introgression=V5) %>%
  pivot_longer(cols=c("bal_intro_almost_recent","bal_intro_ancient","bal_intro_recent","ancient_balancing_selection","introgression")) %>%
  filter(name != "ancient_balancing_selection") %>% filter(name == "introgression") %>% mutate(label="Neutral introgression") %>% 
  ggplot(aes(x=(tol),y=(value),group=label)) + geom_point() + geom_line() +  expand_limits(y=0)+  facet_wrap(~label, scales="free_y",nrow=3) + xlab("Tolerance") + ylab(TeX("Bayes factor")) + theme_bw()  +
  theme(text=element_text(size=24)) + expand_limits(y=0) + geom_hline(yintercept = 1)
p_other_models = as.data.frame(bfs) %>% mutate(tol=c(0.01,seq(0.05,0.5,by=.05))) %>% 
  dplyr::rename(M2=V1, M1=V3, ancient_balancing_selection=V4,M3=V2,introgression=V5) %>%
  pivot_longer(cols=c("M1","M2","M3","ancient_balancing_selection","introgression")) %>%
  filter(name != "ancient_balancing_selection") %>% filter(name != "introgression") %>%
  ggplot(aes(x=(tol),y=(value),group=name)) + geom_point() + geom_line() +  expand_limits(y=0)+  facet_wrap(~name, scales="free_y",nrow=3) + xlab("Tolerance") + ylab(TeX("Bayes factor")) + theme_bw()  +
  theme(text=element_text(size=24)) + expand_limits(y=0) + geom_hline(yintercept = 1)


sumstats$df %>% pivot_longer(cols = c("alt_advantage","cloning_rate","mitotic_recomb_rate","pulse_percentage","pulse_generation","ref_advantage")) %>% ggplot(aes(x=value)) + geom_histogram()  + facet_wrap(~name,scales="free")
ab_list = list()
ab_list2 = list()
show_sims = c()
i  = 1
pc_mat = aa$x[1:(nrow(aa$x)-1),]
gal_pc  = aa$x[nrow(aa$x),]
failed_intro = !is.na(sumstats$sumstats$X1)
num_sims = 2
for(type in unique(sim_type)){
  print(type)
  ab2 = abc(c(epsilon, df_out$V1[!is.na(df_out$V1)]), sumstats$df[sim_type == type,c(1:6)], sumstats$sumstats[sim_type==type,c(1,12:ncol(sumstats$sumstats))],tol=0.01,method="rejection")
  
  show_sims = c(show_sims, (sort(ab2$dist)[1:num_sims]))
  ab_list[[i]] =  ab2
  i = i + 1
}

names(ab_list) = unique(sim_type)


df_prior_posterior = list()
j= 1

d1 =bind_rows(all$pi_between_dfs,.id="NAME")
#d2 = bind_rows(all$output_mutation_list,.id="NAME")
d2 %>%  filter(NAME %in% rownames(ab_list[[4]]$ss)) 
d1$real_data = rep(df_out$V1,length(unique(d1$NAME)))
d1sub = d1 %>% filter(NAME %in% names(show_sims))

library(latex2exp)

good_start= df_out$start2[!is.na(df_out$V1)]
elds = lapply(all$output_mutation_list, function(x){data.frame(eld=unique(x$eld))})
elds =bind_rows(elds,.id = "NAME")


#sumstats$df$name = rownames(sumstats$sumstats)

p2 = d1sub %>% inner_join(sumstats$df,by=c("NAME"="name"))  %>% inner_join(elds,by=c("NAME"="NAME"))%>%  filter(start %in% good_start) %>% filter(sim_type2 == "balancing_selection" | sim_type2 == "introgression") %>% ggplot(aes(y=(pi_between),x=start + 300,color=sim_type2)) +  geom_point(aes(y=real_data,color="real_data")) + geom_point(aes(group=NAME))  + 
  facet_wrap(~paste("eLD =" ,format(eld,digits=3),sep=" "),nrow=2,scales="free_y") + theme_bw() + theme(text=element_text(size=20)) +scale_color_viridis_d(name="Simulation model",labels=c("Ancient balancing selection","Neutral introgression","Observed data")) + ylab(TeX("$\\pi_{between}$")) + xlab("Position") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p_other_m123 =  d1sub %>% inner_join(sumstats$df,by=c("NAME"="name"))  %>% inner_join(elds,by=c("NAME"="NAME"))%>%  filter(start %in% good_start) %>% filter(sim_type != "introgression") %>%
  ggplot(aes(y=(pi_between),x=start + 300,color=sim_type)) +  geom_point(aes(y=real_data,color="real_data")) + geom_point(aes(group=NAME))  + 
  facet_wrap(~paste("eLD =" ,format(eld,digits=3),sep=" "),nrow=2,scales="free_y") + theme_bw() + theme(text=element_text(size=20)) +scale_color_brewer(palette = "Set2", name="Simulation model",labels=c("M2","M1","M3","Ancient balancing selection","Observed data")) + ylab(TeX("$\\pi_{between}$")) + xlab("Position") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



cowplot::plot_grid(p_intro,p2,labels="AUTO",label_size =20)
ggsave("figures/S19.png",width=16,height=8)
cowplot::plot_grid(p_other_models,p_other_m123,labels="AUTO",label_size = 20,rel_widths = c(.5,1))
ggsave("figures/S20.png",width=16,height=8)


d1 = bind_rows(all$pi_between_dfs,.id="NAME") %>% inner_join(sumstats$df,by=c("NAME"="name")) %>% dplyr::mutate(sim_type_three=sim_type) %>% filter(NAME %in% res$name) %>% 
  select(NAME,start,pi_between,sim_type)
d1$real_data = rep(df_out$V1,length(unique(d1$NAME)))

keep_rows = sumstats$df %>% group_by(sim_type2) %>% dplyr::sample_n(20000)
keep_names = keep_rows$name

sumstats$sumstats[rownames(sumstats$sumstats) %in% keep_names,]
a = postpr(c(epsilon,df_out$V1[!is.na(df_out$V1)]),
           (sim_type[rownames(sumstats$sumstats) %in% keep_names]),
           sumstats$sumstats[rownames(sumstats$sumstats) %in% keep_names,c(1,12:ncol(sumstats$sumstats))],tol=0.05,method = "rejection")
j  = 1
library(boa)
for(type in unique(sim_type)){
  
  c = sumstats$df %>% filter(sim_type == type) %>% select(alt_advantage,ref_advantage,cloning_rate,mitotic_recomb_rate,pulse_percentage, pulse_generation)
  c$prior = "Prior"
  d = data.frame(ab_list[[which(names(ab_list) == type)]]$unadj.values )
  d$prior = "Posterior"
  
  print(bind_rows(c,d) %>% pivot_longer(cols=c("alt_advantage","ref_advantage","cloning_rate","mitotic_recomb_rate","pulse_percentage","pulse_generation")) %>%
          ggplot(aes(x=value,fill=prior)) + geom_density(alpha=.5) + facet_wrap(~name,scales="free")  + ggtitle(type))
  j = j + 1
}


p_last_figure = (d) %>% mutate(generations_ago = TOTAL_GENERATIONS*factor - (pulse_generation)*factor) %>%
  pivot_longer(cols=c("alt_advantage","ref_advantage","cloning_rate","mitotic_recomb_rate","pulse_percentage","pulse_generation","generations_ago")) %>% 
  filter(name != "pulse_generation" & name != "generations_ago") %>% mutate(distribution=prior) %>%
  ggplot(aes(x=value)) + geom_density(alpha=.5,fill="grey80") + facet_wrap(~name,scales="free")  + xlab("Parameter value") + ylab("Density")  + theme_bw()+ theme(text=element_text(size=14)) +  theme(axis.text.x=element_text(angle=90, hjust=1))

p_generations_ago = d %>% mutate(generations_ago = TOTAL_GENERATIONS*factor - (pulse_generation)*factor) %>%
  pivot_longer(cols=c("alt_advantage","ref_advantage","cloning_rate","mitotic_recomb_rate","pulse_percentage","pulse_generation","generations_ago")) %>% 
  filter(name == "generations_ago") %>% mutate(distribution=prior) %>%
  ggplot(aes(x=value)) + geom_density(alpha=.5,fill="grey80")  + facet_wrap(~name,scales="free")  + scale_x_reverse() + geom_vline(xintercept = 108e6) +theme_bw() + xlab("Parameter value") + ylab("Density") + theme(text=element_text(size=14))
#sumstats$df %>% 

cowplot::plot_grid(p_generations_ago+theme(legend.position = "none"),p_last_figure,labels="AUTO",label_size = 20)
ggsave("figures/S30.png",width=16,height=8)
#ab2
