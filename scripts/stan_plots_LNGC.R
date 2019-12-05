## clearing the environment
rm(list = ls())  
gc()    

require(rstan)
require(loo)
require(tidyverse)
####################################################################################

## model specific details that needs to be change for every run
modelName <- "ki67_KHM_LNGC"
data_derived1 <- "FM_data.csv"    # name of the file for precursor pop
data_derived2 <- paste("counts_LNGC.csv", sep="")
data_derived3 <- paste("Nfd_LNGC.csv", sep="")
data_derived4 <- paste("ki67_LNGC.csv", sep="") 

## model files (Stan)
# setting data input according to the precurosr pop used
if (grepl("FM", data_derived1) == TRUE){
  #mod_sour <- paste("_P0_", modelName, ".rds", sep = "")     # name of the file for stan model
  mod_sour <- paste("_P0_", modelName, "_", substr(data_derived1, 1,2), ".rds", sep = "")     # name of the file for stan model
} else {
  mod_sour <- paste("_P0_", modelName, ".rds", sep = "")     # name of the file for stan model
}

## setting working dirctory
setwd("/opt/mesh/eigg/sanket/FM_GC_ageCorrected")

## Setting all the directories for opeartions
projectDir <- getwd()
scriptDir <- file.path(projectDir, "scripts", "YFP_informed_models", "LNGC")
modelDir <- file.path(projectDir, "models", "YFP_informed_models", "LNGC")
dataDir <- file.path(projectDir, "data")
toolsDir <- file.path(projectDir, "scripts", "tools")
outputDir <- file.path(projectDir, "output", "YFP_informed_models", "LNGC")
outDir <- file.path(modelDir, paste(modelName, "_", substr(data_derived1, 1,2), sep=""))
saveDir <- file.path(outputDir, paste(modelName, "_", substr(data_derived1, 1,2), sep=""))
figDir <- file.path(projectDir, "deliv", "figures", "YFP_informed_models", "LNGC", paste(modelName, "_", substr(data_derived1, 1,2), sep=""))
tabDir <- file.path(projectDir, "deliv", "tables", "YFP_informed_models", "LNGC", paste(modelName, "_", substr(data_derived1, 1,2), sep=""))


# loadiong the scr# loadiong the script that contains functions for plotting stan parameters
source(file.path(toolsDir, "stanTools.R"))                # save results in new folder

# compiling multiple stan objects together that ran on different nodes
stanfit5 <- readRDS(file.path(saveDir, "35870_P0_ki67_KHM_LNGC_FM.rds"))
stanfit6 <- readRDS(file.path(saveDir, "35871_P0_ki67_KHM_LNGC_FM.rds"))
stanfit7 <- readRDS(file.path(saveDir, "35872_P0_ki67_KHM_LNGC_FM.rds"))
stanfit8 <- readRDS(file.path(saveDir, "35873_P0_ki67_KHM_LNGC_FM.rds"))
stanfit9 <- readRDS(file.path(saveDir, "35874_P0_ki67_KHM_LNGC_FM.rds"))
stanfit10 <- readRDS(file.path(saveDir, "35875_P0_ki67_KHM_LNGC_FM.rds"))
stanfit11 <- readRDS(file.path(saveDir, "35876_P0_ki67_KHM_LNGC_FM.rds"))
stanfit12 <- readRDS(file.path(saveDir, "35877_P0_ki67_KHM_LNGC_FM.rds"))
stanfit13 <- readRDS(file.path(saveDir, "35878_P0_ki67_KHM_LNGC_FM.rds"))
stanfit14 <- readRDS(file.path(saveDir, "35879_P0_ki67_KHM_LNGC_FM.rds"))

fit <- sflist2stanfit(list(stanfit5, stanfit6, stanfit7, stanfit8, stanfit9, stanfit10, stanfit11, stanfit12, stanfit13, stanfit14))

# setting data input according to the precurosr pop used
if (grepl("T2", data_derived1) == TRUE){
  source_list = data.frame("nu" = -5.53e-05, "theta0" = 14.13, "chiEst" = 0.73, "qEst" = 0.06, "eps_donor" = 0.82, "eps_host" = 0.73)
} else if (grepl("T1", data_derived1) == TRUE) {
  source_list = data.frame("nu" = 1e-04, "theta0" = 13.81, "chiEst" = 0.78, "qEst" = 0.07, "eps_donor" = 0.97, "eps_host" = 0.94)
} else if (grepl("FM", data_derived1) == TRUE){
  source_list = data.frame("nu" = -1.1e-03 , "theta0" = 17.20, "chiEst" = 0.74, "qEst" = 0.023, "eps_donor" = 0.31, "eps_host" = 0.14)
} else {
  source_list = data.frame("nu" = 3.30e-04, "theta0" = 14.83, "chiEst" = 0.75, "qEst" = 0.066, "eps_donor" = 0.88, "eps_host" = 0.81)
}

#posterior drwas of parameters
matrix_of_draws <-  as.data.frame(fit) 

source_counts <- function(t){
  exp(source_list$theta0) * exp(-source_list$nu * (t -40))
}

#posterior drwas of parameters
matrix_of_draws <-  as.data.frame(fit) 

mean((matrix_of_draws$psi * source_counts(75)  * (1-matrix_of_draws$f_fast))/Y1pred1$median[Y1pred1$timeseries == 75])
mean((matrix_of_draws$psi * source_counts(300) * (matrix_of_draws$f_fast))/Y1pred1$median[Y1pred1$timeseries == 300])


print("psiF")
mean((matrix_of_draws$psi * source_counts(75)  * (0.2255149))/Y1pred1$median[Y1pred1$timeseries == 75]) * 100
quantile((matrix_of_draws$psi * source_counts(75)  * (0.05931079))/Y1pred1$median[Y1pred1$timeseries == 75], probs = 0.025) * 100
quantile((matrix_of_draws$psi * source_counts(75)  * (0.3889095))/Y1pred1$median[Y1pred1$timeseries == 75], probs = 0.975) * 100

mean((matrix_of_draws$psi * source_counts(300) * (0.1648974))/Y1pred1$median[Y1pred1$timeseries == 300]) * 100
quantile((matrix_of_draws$psi * source_counts(300)  * (0.02671057))/Y1pred1$median[Y1pred1$timeseries == 300], probs = 0.025) * 100
quantile((matrix_of_draws$psi * source_counts(300)  * (0.3691268))/Y1pred1$median[Y1pred1$timeseries == 300], probs = 0.975) * 100

print("psiS")
mean((matrix_of_draws$psi * source_counts(75)  * (0.7744851))/Y1pred1$median[Y1pred1$timeseries == 75]) * 100
quantile((matrix_of_draws$psi * source_counts(75)  * (0.9406892))/Y1pred1$median[Y1pred1$timeseries == 75], probs = 0.025) * 100
quantile((matrix_of_draws$psi * source_counts(75)  * (0.6110905))/Y1pred1$median[Y1pred1$timeseries == 75], probs = 0.975) * 100

mean((matrix_of_draws$psi * source_counts(300)  * (0.8351026))/Y1pred1$median[Y1pred1$timeseries == 300]) * 100
quantile((matrix_of_draws$psi * source_counts(300)  * (0.9732894))/Y1pred1$median[Y1pred1$timeseries == 300], probs = 0.025) * 100
quantile((matrix_of_draws$psi * source_counts(300)  * (0.6308732))/Y1pred1$median[Y1pred1$timeseries == 300], probs = 0.975) * 100

print("alphaF")
mean(matrix_of_draws$alpha_calc)
quantile(matrix_of_draws$alpha_calc, probs = c(0.025, 0.975))

print("alphaS")
mean(1-matrix_of_draws$alpha_calc)
quantile(1-matrix_of_draws$alpha_calc, probs = c(0.025, 0.975))

print("Beta")
mean(matrix_of_draws$Beta)
quantile(matrix_of_draws$Beta, probs = c(0.025, 0.975))

print("rhoF")
mean(matrix_of_draws$rhoFast_inv)
quantile(matrix_of_draws$rhoFast_inv, probs = c(0.025, 0.975))

print("lambdaF")
mean(matrix_of_draws$lambdaFast_inv)
quantile(matrix_of_draws$lambdaFast_inv, probs = c(0.025, 0.975))

print("deltaF")
mean(matrix_of_draws$deltaFast_inv)
quantile(matrix_of_draws$deltaFast_inv, probs = c(0.025, 0.975))

print("lambdaS")
mean(matrix_of_draws$lambdaSlow_inv)
quantile(matrix_of_draws$lambdaSlow_inv, probs = c(0.025, 0.975))

print("deltaS")
mean(1/(matrix_of_draws$lambdaSlow + matrix_of_draws$rhoSlow))
quantile(matrix_of_draws$deltaSlow_inv, probs = c(0.025, 0.975))

print("rhoS")
mean(1/matrix_of_draws$rhoSlow)
quantile(1/matrix_of_draws$rhoSlow, probs = c(0.025, 0.975))




# finding the parameters used in the model 
# using the last parameter("sigma4") in the array to get the total number of parameters set in the model
num_pars <- which(fit@model_pars %in% "sigma4")      # the variable "sigma4" will change depdending on the data used
parametersToPlot <- fit@model_pars[1:num_pars]

# number of post-burnin samples that are used for plotting 
nPost <- nrow(fit)

################################################################################################
################################################################################################

## loading required datasets for plotting
source_sorted <- read_csv(file.path(dataDir, data_derived1))%>% arrange(age.at.S1K)
counts_sorted <- read_csv(file.path(dataDir, data_derived2))%>% arrange(age.at.S1K)
Nfd_sorted <- read_csv(file.path(dataDir, data_derived3))%>% arrange(age.at.S1K)
ki67_sorted <- read_csv(file.path(dataDir, data_derived4))%>% arrange(age.at.S1K)

#################################################################################################
# calculating PSIS-L00-CV for the fit
counts_loglik <- extract_log_lik(fit, parameter_name = "log_lik1", merge_chains = TRUE)
fd_loglik <- extract_log_lik(fit, parameter_name = "log_lik2", merge_chains = TRUE)
hostki_loglik <- extract_log_lik(fit, parameter_name = "log_lik3", merge_chains = TRUE)
donorki_loglik <- extract_log_lik(fit, parameter_name = "log_lik4", merge_chains = TRUE)

combined_loglik <- extract_log_lik(fit, parameter_name = "log_lik", merge_chains = TRUE)
log_lik_comb <- cbind(counts_loglik, fd_loglik, donorki_loglik, hostki_loglik)


# optional but recommended
#ll_array <- extract_log_lik(fit,parameter_name = "log_lik1", merge_chains = FALSE)
#r_eff <- relative_eff(exp(ll_array))

# loo-ic values
loo_loglik <- loo(log_lik_comb, save_psis = FALSE, cores = 8)

# Widely applicable AIC
AICw_lok <- waic(combined_loglik)

# AIC from LLmax
#AIC_lok <-  -2 * max(combined_loglik)  + 2 * length(parametersToPlot)

ploocv <- rbind("loo-ic"=loo_loglik$estimates[3], "WAIC" = AICw_lok$estimates[3])  #, "AIC" = AIC_lok)
ploocv[1]

### posterior distributions of parameters
ptable <- monitor(as.array(fit, pars = parametersToPlot), warmup = 0, print = FALSE)
ptable[,c(1,4,8)]
################################################################################################
# save results in output directory
if (!file.exists(figDir)){
  dir.create(figDir)
}

if (!file.exists(tabDir)){
  dir.create(tabDir)
}

## writting the parameter values with CIs into a csv file
write.csv(ptable, file = file.path(tabDir, paste(modelName, "_ParameterTable.csv", sep = "")))
#write.csv(parameter_est, file = file.path(tabDir, paste(modelName, "_EstimatesTable.csv", sep = "")))
write.csv(ploocv, file = file.path(tabDir, paste(modelName, "_Informatiion_criteria.csv", sep = "")))



## open graphics device 
# saving  plots for quality control 
#pdf(file = file.path(figDir, paste(modelName,"StanPlots%03d.pdf", sep = "")),
 #   width = 8, height = 5, onefile = FALSE, useDingbats = FALSE)

#pairs(fit, pars = parametersToPlot)
#### monitoring stan sampling process
#### rhat should be close to 1
#rhats <- rhat(fit, pars = parametersToPlot)
#mcmc_rhat(rhats) + yaxis_text() + myTheme

#ratios1 <- neff_ratio(fit, pars = parametersToPlot)
#mcmc_neff(ratios1) + yaxis_text() + myTheme

#### to view divergent conditions during the run
#posterior <- as.array(fit)
#mcmc_acf(posterior, pars = parametersToPlot) + myTheme

#mcmcHistory(fit, pars = parametersToPlot, nParPerPage = 4, myTheme = myTheme)

#### density plots for paramaeters
#mcmc_dens_overlay(posterior, parametersToPlot)
#mcmc_dens(posterior, parametersToPlot) + myTheme

#dev.off()

################################################################################################
################################################################################################
## posterior predictive distributions
options(bayesplot.base_size = 15,
        bayesplot.base_family = "sans")
#color_scheme_set(scheme = "viridis")
myTheme <- theme(text = element_text(size = 10), axis.text = element_text(size = 10), axis.title =  element_text(size = 10, face = "bold"),
                 plot.title = element_text(size=10,  hjust = 0.5, face = "bold"),
                 legend.background = element_blank(), legend.key = element_blank())

# setting ggplot theme for rest fo the plots
theme_set(theme_bw())

fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # remove + after exponent, if exists. E.g.: (e^+2 -> e^2)
  l <- gsub("e\\+","e",l)  
  # turn the 'e' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # convert 1x10^ or 1.000x10^ -> 10^
  l <- gsub("\\'1[\\.0]*\\'\\%\\*\\%", "", l)
  # return this as an expression
  parse(text=l)
}

log10minorbreaks=as.numeric(1:10 %o% 10^(4:8))

# time sequence for predictions specific to age bins within the data
# In this analysis 3 age bins were selected with mean ages for each bin as 48, 72 and 120 respectively.
counts_agebmt_sorted <- counts_sorted %>% arrange(age.at.bmt)
Mean_agebin1 <- round(mean(counts_agebmt_sorted$age.at.bmt[1:11]))
Mean_agebin2 <- round(mean(counts_agebmt_sorted$age.at.bmt[12:16]))
Mean_agebin3 <- round(mean(counts_agebmt_sorted$age.at.bmt[17:23]))

ts_pred1 <- seq(Mean_agebin1, to = 600, by = 1)
tb_pred1 <- rep(Mean_agebin1, length(ts_pred1))
ts_pred2 <- seq(Mean_agebin2, to = 600, by = 1)
tb_pred2 <- rep(Mean_agebin2, length(ts_pred2))
ts_pred3 <- seq(Mean_agebin3, to = 600, by = 1)
tb_pred3 <- rep(Mean_agebin3, length(ts_pred3))

# tb_time_pred -- for solving for initial conditions at each tb (here the mean age of each age bin)
tb_time_pred1 <- c(min(counts_sorted$age.at.bmt), Mean_agebin1)
tb_time_pred2 <- c(min(counts_sorted$age.at.bmt), Mean_agebin2)
tb_time_pred3 <- c(min(counts_sorted$age.at.bmt), Mean_agebin3)

# Total cell counts
counts_binned <- counts_sorted%>%
  mutate(age_bins = ifelse(age.at.bmt <= 56, "Age bin#1: <8Wks",
                           ifelse(age.at.bmt <= 84, "Age bin#2: 8-12Wks", "Age bin#3: >12Wks")))

# normalised donor fractions
Nfd_binned <- Nfd_sorted %>%
  mutate(age_bins = ifelse(age.at.bmt <= 56, "Age bin#1: <8Wks",
                           ifelse(age.at.bmt <= 84, "Age bin#2: 8-12Wks", "Age bin#3: >12Wks")))

## plotting ki67 predictions
ki67_binned <- read_csv(file.path(dataDir, data_derived4))%>%
  gather(- c(days.post.bmt, age.at.S1K, age.at.bmt), value = "prop_ki67hi", key = "subpopulation")%>%
  mutate(age_bins = ifelse(age.at.bmt <= 56, "Age bin#1: <8Wks",
                           ifelse(age.at.bmt <= 84, "Age bin#2: 8-12Wks", "Age bin#3: >12Wks")))

ploocv[1]

# sourceing the file that has handwritten lengthy code for plotting observations
source("scripts/GC_plots.R")

## saving  plots for quality control 
pdf(file = file.path(figDir, paste(modelName,"IndPlots%03d.pdf", sep = "")),
    width = 7, height = 8, onefile = FALSE, useDingbats = FALSE)

cowplot::plot_grid(counts_plot, nfd_plot, ki67_plot, labels = c("A", "B", "C"),  nrow =  3)
cowplot::plot_grid(counts_plot, nfd_plot, ki67_plot2, labels = c("A", "B", "C"),  nrow =  3)

dev.off()


pdf(file = file.path(figDir, paste(modelName,"Plots%03d.pdf", sep = "")),
    width = 7, height = 5, onefile = FALSE, useDingbats = FALSE )

lay1 <- rbind(c(1,1,2,2),
              c(NA,3,3,NA))

gridExtra::grid.arrange(counts_facet, nfd_comb, ki67_comb, layout_matrix = lay1)

toprow <- cowplot::plot_grid(counts_facet, nfd_comb,  labels = c("A", "B"), ncol = 2)
cowplot::plot_grid(toprow, ki67_plot,  labels = c("", "C"), nrow = 2)

dev.off()
