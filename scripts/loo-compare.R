## clearing the environment
rm(list = ls())  
gc()    

require(rstan)
require(loo)
require(tidyverse)
####################################################################################
## setting working dirctory
setwd("/opt/mesh/eigg/sanket/FM_GC_ageCorrected")

#### population of interest
arg1 <- "SPGC"

## Setting all the directories for opeartions
projectDir <- getwd()
scriptDir <- file.path(projectDir, "scripts", "YFP_informed_models", arg1)
modelDir <- file.path(projectDir, "models", "YFP_informed_models",  arg1)
dataDir <- file.path(projectDir, "data")
toolsDir <- file.path(projectDir, "scripts", "tools")
outputDir <- file.path(projectDir, "output", "YFP_informed_models",  arg1)

## getting the list of model and stan files saved 
modelnames_list <- list.files(outputDir)

func_for_looic <- function(modelname){
  
  saved_stan_file <- list.files(file.path(outputDir, modelname), pattern = "\\.rds")
  
  path_for_stanfile <- file.path(outputDir, modelname, saved_stan_file)
  
  object_list <- 
  for (i in 1:length(saved_stan_file)){
    
    object_list[i] <- read_rds(path_for_stanfile[i])
  }
  
  
  fit_obj <- sflist2stanfit(list(object_list))
  
  # calculating PSIS-L00-CV for the fit
  counts_loglik <- extract_log_lik(fit_obj, parameter_name = "log_lik1", merge_chains = TRUE)
  fd_loglik <- extract_log_lik(fit_obj, parameter_name = "log_lik2", merge_chains = TRUE)
  hostki_loglik <- extract_log_lik(fit_obj, parameter_name = "log_lik3", merge_chains = TRUE)
  donorki_loglik <- extract_log_lik(fit_obj, parameter_name = "log_lik4", merge_chains = TRUE)
  
  
  #combined_loglik <- extract_log_lik(fit, parameter_name = "log_lik", merge_chains = TRUE)
  log_lik_comb <- cbind(counts_loglik, fd_loglik, donorki_loglik, hostki_loglik)
  
  # optional but recommended
  ll_array <- extract_log_lik(fit_obj,parameter_name = "log_lik1", merge_chains = FALSE)
  r_eff <- relative_eff(exp(ll_array))
  
  # loo-ic values
  loo_loglik <- loo(log_lik_comb, save_psis = FALSE, cores = 8)
  
  return(loo_loglik$estimates[3])
  
}

modelname <- modelnames_list[8]

func_for_looic(m)

loo_ic_values <- c()
for (i in 1:3){
  
  loo_ic_values[i] <- func_for_looic(modelnames_list[i])
}

## function for calculating akaike weight
model_wt <- function(x) { exp(-0.5 * x) }

## list of looic values
loo_list <- loo_ic_values

### applying the function model_wt to list of looic values
model_list <- cbind(lapply(loo_list, model_wt))

### summing the wights of all models to find the relative avearge for each model 
tot_sum <- 0
model_avg <- function(x){
  for (i in 1:length(loo_list)){
    
    tot_sum <- model_list[[i]] + tot_sum}
  
  ans = exp(-0.5 * x) * 100 /tot_sum
  return(ans)}

### list of relative model weights
wt_list <- cbind(lapply(loo_list, model_avg))

write.csv(wt_list, file = file.path(outputDir))

### confirmng that all model weights sum upto 100
tot2_sum = 0
for (i in 1:length(loo_list)){
  
  tot2_sum <- wt_list[[i]] + tot2_sum}

tot2_sum


