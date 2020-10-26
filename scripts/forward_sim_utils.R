library(tidyverse)
library(abc)
source("scripts/S17_real_data_conversion.R")

#epsilon = 0.5916645

get_sumstats = function(test){
  elds = unlist(lapply(test$output_mutation_list, function(x) { x$eld[1]}))
  idx_bad = (unlist(lapply(test$output_mutation_list, function(x) {nrow(x)})) < 6)
  idx_good = (unlist(lapply(test$output_mutation_list, function(x) {nrow(x)})) == 6)
  
  
  #change_points  = unlist(lapply(intro$pi_between_dfs[idx_good],function(x){a = cpt.mean(x$pi_between,penalty="AIC",method="PELT");length(a@param.est$mean) - 1 }))
  
  
  ds = rep(NA,length(test$output_mutation_list))
  slopes = rep(NA,length(test$output_mutation_list))
  predict  =  rep(NA,length(test$output_mutation_list))
  var_lag =  rep(NA,length(test$output_mutation_list))
  pi_median = rep(NA,length(test$output_mutation_list))
  
  maxes =  data.frame(max1=rep(NA,length(test$output_mutation_list)),max2=rep(NA,length(test$output_mutation_list)),max3=rep(NA,length(test$output_mutation_list)),
                      max4=rep(NA,length(test$output_mutation_list)),max5=rep(NA,length(test$output_mutation_list)),max6=rep(NA,length(test$output_mutation_list)))
  idx3 = which(!is.na(df_out$V1))
  maxes_df = data.frame(matrix(nrow=length(test$output_mutation_list),ncol=length(idx3)))
  ds[idx_bad] = 0
  slopes[idx_bad] = 0
  predict[idx_bad] = 0
  var_lag[idx_bad] = 0
  pi_median[idx_bad] =0 
  #maxes = rep(NA,length(test$output_mutation_list))
  idx2 = which(df_out$V1 %in% out_df$max_ds)
  idx2 = idx2[-4]
  
  maxes2 = (lapply(test$pi_between_dfs[idx_good], function(x) {c((x$pi_between[idx2]))}))
  maxes2 = t(bind_rows(maxes2))
  # which(unlist(lapply(test$pi_between_dfs[idx_good], function(x){length(x$pi_between)})) != 882)
  maxes3=  (lapply(test$pi_between_dfs[idx_good], function(x) {c((x$pi_all[idx3]))}))
  g2 = t(bind_rows(maxes3))
  
  pi_median2 = (unlist(lapply(test$pi_between_dfs[idx_good], function(x) {median(x$pi_between)})))
  ds2 = unlist(lapply(test$slopes_list[idx_good], function(x) { mean(x$max_ds)}))
  slopes2 = unlist(lapply(test$slopes_list[idx_good], function(x) { mean(abs(x$slope))}))
  predict2 = unlist(lapply(test$slopes_list[idx_good], function(x) { mean((x$predict))}))
  ds[idx_good] = ds2
  slopes[idx_good] = slopes2
  predict[idx_good] = predict2
  pi_median[idx_good] = pi_median2
  maxes[idx_good,] = maxes2
  print(slopes2)
  var_lag2 = unlist(lapply(test$pi_between_dfs[idx_good], function(x) { var(x$pi_between - lag(x$pi_between),na.rm=T)}))
  maxes_df[idx_good,] = g2
  var_lag[idx_good] = var_lag2
  
  sumstats = data.frame(eld=elds,slopes=slopes,predict=ds,var_lag=var_lag,pi_median=pi_median,maxes,maxes_df)
  #sumstats = sumstats[which(!(1:nrow(sumstats) %in% idx_rm)),]
  df = bind_rows(test$output_parameters,.id="sim") %>% pivot_wider(names_from=parameter,values_from=value) %>% 
    dplyr::select(alt_advantage, ref_advantage, cloning_rate,mitotic_recomb_rate,pulse_percentage,pulse_generation,sim_type) %>% mutate(alt_advantage=as.numeric(alt_advantage),
                                                                                                                                        ref_advantage=as.numeric(ref_advantage),
                                                                                                                                        cloning_rate=as.numeric(cloning_rate),
                                                                                                                                        mitotic_recomb_rate=as.numeric(mitotic_recomb_rate),
                                                                                                                                        pulse_percentage=as.numeric(pulse_percentage),
                                                                                                                                        pulse_generation=as.numeric(pulse_generation),sim_type=sim_type)
  
  
  return(list(sumstats=sumstats,df=df))
} 



#aa = readRDS("~/HOFFMAN/slim/gal_bal/scripts/sim_output/test.out.sim.rds")
arun = readRDS("data/arun2.rds")
all = readRDS("data/intro_new.rds")
intro = readRDS("data/intro5.rds")
sumstats2 = get_sumstats(intro)

### Get idx of introgressions ### 

intro_get_all = unlist(lapply(intro$output_parameters, function(x) {x$value[x$parameter == "sim_type"]}))
keep_idx = which(intro_get_all == "introgression")
names_old_intro = paste(names(intro$output_mutation_list),"intro",sep="_")
names(intro$output_mutation_list) = names_old_intro
names(intro$output_parameters) = names_old_intro
names(intro$slopes_list) = names_old_intro
names(intro$pi_between_dfs) = names_old_intro

intro$output_mutation_list = intro$output_mutation_list[keep_idx]
intro$output_parameters = intro$output_parameters[keep_idx]
intro$slopes_list = intro$slopes_list[keep_idx]
intro$pi_between_dfs = intro$pi_between_dfs[keep_idx]



names_new_arun = paste(names(arun$output_mutation_list),"arun",sep="_")
names(arun$output_mutation_list) = names_new_arun
names(arun$output_parameters) = names_new_arun
names(arun$slopes_list) = names_new_arun
names(arun$pi_between_dfs) = names_new_arun




all$output_mutation_list = c(all$output_mutation_list,arun$output_mutation_list,intro$output_mutation_list)
all$output_parameters = c(all$output_parameters,arun$output_parameters,intro$output_parameters)
all$slopes_list = c(all$slopes_list,arun$slopes_list, intro$slopes_list)
all$pi_between_dfs = c(all$pi_between_dfs,arun$pi_between_dfs,intro$pi_between_dfs)
sumstats = get_sumstats(all)
factor = 50000
TOTAL_GENERATIONS=3.2e9/factor + 8 * 10e6/factor
max_generations = c()
#intro$output_mutation_list
sumstats$df$sim_type[grep("arun",rownames(sumstats$sumstats) )] = "balancing_selection_introgression_ancient"

old_sim_type = sumstats$df$sim_type

sumstats$df$name = rownames(sumstats$sumstats)

ranges_recent = c(TOTAL_GENERATIONS - 1, TOTAL_GENERATIONS - 25e6/5e4)
ranges_almost_recent = c(TOTAL_GENERATIONS - 25e6/5e4 - 1,TOTAL_GENERATIONS - 50e6/5e4)
ranges_old = c(TOTAL_GENERATIONS - 50e6/5e4 -1,TOTAL_GENERATIONS - 100e7/factor)
ranges_list = list(ranges_recent=ranges_recent, ranges_almost_recent=ranges_almost_recent, ranges_old=ranges_old)
sumstats$df = sumstats$df %>% 
  mutate(sim_type2 = ifelse(sim_type == "balancing_selection_introgression" & pulse_generation <= ranges_recent[1] & pulse_generation >= ranges_recent[2],"bal_intro_recent",sim_type)) %>% 
  mutate(sim_type2 = ifelse(sim_type == "balancing_selection_introgression" & pulse_generation <= ranges_almost_recent[1] & pulse_generation >= ranges_almost_recent[2],"bal_intro_almost_recent",sim_type2))  %>% 
  mutate(sim_type2 = ifelse(sim_type == "balancing_selection_introgression" & pulse_generation <= ranges_almost_recent[1] & pulse_generation >= ranges_almost_recent[2],"bal_intro_almost_recent",sim_type2)) %>%
  mutate(sim_type2= ifelse(sim_type == "balancing_selection_introgression_ancient" & pulse_generation <= ranges_old[1] & pulse_generation >= ranges_old[2],"bal_intro_ancient",sim_type2)) %>% 
  mutate(sim_type2 = ifelse(sim_type == "balancing_selection_introgression_ancient" & pulse_generation <= ranges_recent[1] & pulse_generation >= ranges_recent[2],"bal_intro_recent",sim_type2)) %>% 
  mutate(sim_type2 = ifelse(sim_type == "balancing_selection_introgression_ancient" & pulse_generation <= ranges_almost_recent[1] & pulse_generation >= ranges_almost_recent[2],"bal_intro_almost_recent",sim_type2)) 
sim_type = sumstats$df$sim_type2 
sim_type2 = sumstats$df
sumstats$df$sim_type = sim_type
sumstats$df$name = rownames(sumstats$sumstats)
#idx1 = sumstats$df$pulse_generation < 65100 | sumstats$df$sim_type != "balancing_selection_introgression"
library(abc)
#sumstats$sumstats = sumstats$sumstats[idx1,]
#sumstats$df = sumstats$df[idx1,]
sim_type = sumstats$df$sim_type
count_seg = (unlist(lapply(all$output_mutation_list, function(x) {nrow(x)})))
sumstats$df %>% pivot_longer(cols=c("alt_advantage","ref_advantage","cloning_rate","mitotic_recomb_rate","pulse_percentage","pulse_generation"),names_to = c("parameter"),values_to = c("val")) %>%
  ggplot(aes(x=val)) + geom_histogram()+ facet_wrap(~parameter*sim_type,scales="free",ncol=6)
params = c(epsilon, df_out$V1[!is.na(df_out$V1)])