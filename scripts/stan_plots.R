## clearing the environment
rm(list = ls())  
gc()    

require(rstan)
require(loo)
require(tidyverse)
####################################################################################

## model specific details that needs to be change for every run
modelName <- "ki67_SHM_FM"
data_derived1 <- "T1_data.csv"    # name of the file for precursor pop
data_derived2 <- paste("counts_FM.csv", sep="")
data_derived3 <- paste("Nfd_FM.csv", sep="")
data_derived4 <- paste("ki67_FM.csv", sep="") 

## setting working dirctory
setwd("/opt/mesh/eigg/sanket/FM_GC_ageCorrected")

## Setting all the directories for opeartions
projectDir <- getwd()
scriptDir <- file.path(projectDir, "scripts")
modelDir <- file.path(projectDir, "models")
dataDir <- file.path(projectDir, "data")
toolsDir <- file.path(scriptDir, "tools")
outDir <- file.path(modelDir, paste(modelName, "_", substr(data_derived1, 1,2), sep=""))
outputDir <- file.path(projectDir, "output")
saveDir <- file.path(outputDir, paste(modelName, "_", substr(data_derived1, 1,2), sep=""))
figDir <- file.path(projectDir, "deliv", "figures", paste(modelName, "_", substr(data_derived1, 1,2), sep=""))
tabDir <- file.path(projectDir, "deliv", "tables", paste(modelName, "_", substr(data_derived1, 1,2), sep=""))

# loadiong the scr# loadiong the script that contains functions for plotting stan parameters
source(file.path(toolsDir, "stanTools.R"))                # save results in new folder

# compiling multiple stan objects together that ran on different nodes
stanfit1 <- readRDS(file.path(saveDir, "3775_P0_ki67_KHM_FM.rds"))
stanfit2 <- readRDS(file.path(saveDir, "3775_P1_ki67_KHM_FM.rds"))
stanfit3 <- readRDS(file.path(saveDir, "3775_P2_ki67_KHM_FM.rds"))
stanfit4 <- readRDS(file.path(saveDir, "3775_P3_ki67_KHM_FM.rds"))

fit <- sflist2stanfit(list(stanfit1, stanfit2, stanfit3, stanfit4))

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

# ################################################################################################
# calculating PSIS-L00-CV for the fit
combined_loglik <- extract_log_lik(fit, parameter_name = "log_lik", merge_chains = TRUE)

# optional but recommended
ll_array <- extract_log_lik(fit, merge_chains = FALSE)
r_eff <- relative_eff(exp(ll_array))

# loo-ic values
loo_loglik <- loo(combined_loglik, r_eff = r_eff, save_psis = FALSE, cores = 8)

# Widely applicable AIC
AICw_lok <- waic(combined_loglik)

# AIC from LLmax
#AIC_lok <-  -2 * max(combined_loglik)  + 2 * length(parametersToPlot)

ploocv <- rbind("loo-ic"=loo_loglik$estimates[3], "WAIC" = AICw_lok$estimates[3])  #, "AIC" = AIC_lok)
ploocv[1]

### posterior distributions of parameters
ptable <- monitor(as.array(fit, pars = parametersToPlot), warmup = 0, print = FALSE)
1/ptable[,c(1,4,8)]

################################################################################################
dir.create(figDir)
dir.create(tabDir)

## writting the parameter values with CIs into a csv file
write.csv(ptable, file = file.path(tabDir, paste(modelName, "_ParameterTable.csv", sep = "")))
write.csv(ploocv, file = file.path(tabDir, paste(modelName, "_Informatiion_criteria.csv", sep = "")))

## open graphics device 
## saving  plots for quality control 
pdf(file = file.path(figDir, paste(modelName,"StanPlots%03d.pdf", sep = "")),
    width = 8, height = 5, onefile = FALSE)

#pairs(fit, pars = parametersToPlot)

options(bayesplot.base_size = 15,
        bayesplot.base_family = "sans")
color_scheme_set(scheme = "viridis")
myTheme <- theme(text = element_text(size = 11), axis.text = element_text(size = 11), axis.title =  element_text(size = 11, face = "bold"),
                 plot.title = element_text(size=11,  hjust = 0.5, face = "bold"),
                 legend.background = element_blank(), legend.key = element_blank(),
                 legend.text = element_text(size=9), legend.title = element_text(9))

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

#rhats <- rhat(fit, pars = parametersToPlot)
mcmc_rhat(rhats) + yaxis_text() + myTheme

ratios1 <- neff_ratio(fit, pars = parametersToPlot)
mcmc_neff(ratios1) + yaxis_text() + myTheme

posterior <- as.array(fit)
mcmc_acf(posterior, pars = parametersToPlot) + myTheme

mcmcHistory(fit, pars = parametersToPlot, nParPerPage = 4, myTheme = myTheme)

mcmc_dens_overlay(posterior, parametersToPlot)
mcmc_dens(posterior, parametersToPlot) + myTheme

dev.off()

################################################################################################
################################################################################################
## posterior predictive distributions

# time sequence for predictions specific to age bins within the data
# In this analysis 3 age bins were selected with mean ages for each bin as 48, 72 and 120 respectively.
ts_pred1 <- seq(48, to = 600, by = 1)
tb_pred1 <- rep(48, length(ts_pred1))
ts_pred2 <- seq(72, to = 600, by = 1)
tb_pred2 <- rep(72, length(ts_pred2))
ts_pred3 <- seq(120, to = 600, by = 1)
tb_pred3 <- rep(120, length(ts_pred3))

# tb_time_pred -- for solving for initial conditions at each tb (here the mean age of each age bin)
tb_time_pred1 <- c(min(counts_sorted$age.at.bmt), 48)
tb_time_pred2 <- c(min(counts_sorted$age.at.bmt), 72)
tb_time_pred3 <- c(min(counts_sorted$age.at.bmt), 120)

# Total cell counts
counts_binned <- counts_sorted%>%
  mutate(age_bins = ifelse(age.at.bmt <= 56, "Age bin#1: <8Wks",
                           ifelse(age.at.bmt <= 84, "Age bin#2: 8-12Wks", "Age bin#3: >12Wks")))

Cpred1 <- as.data.frame(fit, pars = "countspred_age1") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred1,
            "age_bins" = rep("Age bin#1: <8Wks", length(ts_pred1)))

Y1pred1 <- as.data.frame(fit, pars = "y1_mean_pred_age1") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred1,
            "age_bins" = rep("Age bin#1: <8Wks", length(ts_pred1)))

# Total cell counts
Cpred2 <- as.data.frame(fit, pars = "countspred_age2") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred2,
            "age_bins" = rep("Age bin#2: 8-12Wks", length(ts_pred2)))

Y1pred2 <- as.data.frame(fit, pars = "y1_mean_pred_age2") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred2,
            "age_bins" = rep("Age bin#2: 8-12Wks", length(ts_pred2)))

# Total cell counts
Cpred3 <- as.data.frame(fit, pars = "countspred_age3") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred3 -7,
            "age_bins" = rep("Age bin#3: >12Wks", length(ts_pred3)))

Y1pred3 <- as.data.frame(fit, pars = "y1_mean_pred_age3") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred3-7,
            "age_bins" = rep("Age bin#3: >12Wks", length(ts_pred3)))

#gathering data frames
counts_median_pooled <- rbind(Y1pred1, Y1pred2, Y1pred3)
counts_sd_pooled <- rbind(Cpred1, Cpred2, Cpred3)

#### plots for cell counts
counts_facet <- ggplot() +
  geom_line(data = counts_median_pooled, aes(x = timeseries, y = median, col = age_bins)) +
  geom_ribbon(data = counts_median_pooled, aes(x = timeseries, ymin = lb, ymax = ub, fill = age_bins), alpha = 0.25)+
  geom_point(data = counts_binned, aes(x = age.at.S1K, y = total_counts, col = age_bins)) +
  scale_color_manual(values=c("#356af2", "#28af46", "#ea2745"), name = NULL,
                     labels = c("<8Wks", "8-12Wks", ">12Wks"))+
  scale_fill_manual(values=c("#356af2", "#28af46", "#ea2745"), name = NULL,
                    labels = c("<8Wks", "8-12Wks", ">12Wks"))+
  labs(title=paste('Cell counts: ', substr(modelName, 10, 13)),  y=NULL, x="Host age (days)") + 
  scale_x_continuous(limits = c(60, 600), trans = "log10", breaks = c(10,30,100,300))+
  scale_y_continuous(limits = c(5e6, 2e8), trans="log10", breaks=c(1e4, 1e5, 1e6, 1e7, 1e8), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  guides(color =FALSE, fill = FALSE) + myTheme + theme(legend.position = c(0.5, 0.85), legend.direction = "horizontal")

# normalised donor fractions
Nfd_binned <- Nfd_sorted %>%
  mutate(age_bins = ifelse(age.at.bmt <= 56, "Age bin#1: <8Wks",
                           ifelse(age.at.bmt <= 84, "Age bin#2: 8-12Wks", "Age bin#3: >12Wks")))


fdpred1 <- as.data.frame(fit, pars = "fdpred_age1") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred1,
            "age_bins" = rep("Age bin#1: <8Wks", length(ts_pred1)))%>%
  filter(timeseries <= 575)

Y2pred1 <- as.data.frame(fit, pars = "y2_mean_pred_age1") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred1,
            "age_bins" = rep("Age bin#1: <8Wks", length(ts_pred1)))%>%
  filter(timeseries <= 575)

fdpred2 <- as.data.frame(fit, pars = "fdpred_age2") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred2,
            "age_bins" = rep("Age bin#2: 8-12Wks", length(ts_pred2)))%>%
  filter(timeseries <= 575)

Y2pred2 <- as.data.frame(fit, pars = "y2_mean_pred_age2") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred2,
            "age_bins" = rep("Age bin#2: 8-12Wks", length(ts_pred2)))%>%
  filter(timeseries <= 575)


fdpred3 <- as.data.frame(fit, pars = "fdpred_age3") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred3 - 7,
            "age_bins" = rep("Age bin#3: >12Wks", length(ts_pred3)))%>%
  filter(timeseries <= 575)

Y2pred3 <- as.data.frame(fit, pars = "y2_mean_pred_age3") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred3 - 7,
            "age_bins" = rep("Age bin#3: >12Wks", length(ts_pred3)))%>%
  filter(timeseries <= 575)


#gathering data frames
nfd_median_pooled <- rbind(Y2pred1, Y2pred2, Y2pred3)
nfd_sd_pooled <- rbind(fdpred1, fdpred2, fdpred3)

nfd_facet <- ggplot() +
  geom_line(data = nfd_median_pooled, aes(x = timeseries, y = median, col = age_bins)) +
  geom_ribbon(data = nfd_median_pooled, aes(x = timeseries, ymin = lb, ymax = ub, fill = age_bins), alpha = 0.25)+
  geom_point(data = Nfd_binned, aes(x = age.at.S1K, y = Nfd, col = age_bins)) +
  scale_color_manual(values=c("#356af2", "#28af46", "#ea2745"), name = NULL,
                     labels = c("<8Wks", "8-12Wks", ">12Wks"))+
  scale_fill_manual(values=c("#356af2", "#28af46", "#ea2745"), name = NULL,
                    labels = c("<8Wks", "8-12Wks", ">12Wks"))+
  guides(fill = FALSE) + myTheme + theme(legend.position = c(0.8, 0.3))

# combined plot
nfd_comb <- nfd_facet +
  geom_hline(aes(yintercept = 1), linetype = 2)+
  labs(x = "Host age (days)", y = NULL, title = paste("Normalised Donor fractions: ", substr(modelName, 10, 13))) +
  scale_x_continuous(limits = c(40, 575), breaks = c(0,150,300,450, 600)) +
  scale_y_continuous(limits =c(0, 1.1), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) 


## plotting ki67 predictions
ki67_binned <- read_csv(file.path(dataDir, data_derived4))%>%
  gather(- c(Lamis.ID, days.post.bmt, age.at.S1K, age.at.bmt), value = "prop_ki67hi", key = "subpopulation")%>%
  mutate(age_bins = ifelse(age.at.bmt <= 56, "Age bin#1: <8Wks",
                           ifelse(age.at.bmt <= 84, "Age bin#2: 8-12Wks", "Age bin#3: >12Wks")))

donor_kiPred1 <- as.data.frame(fit, pars = "donor_kiprop_pred1") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred1 +1,
            "age_bins" = rep("Age bin#1: <8Wks", length(ts_pred1)))%>%
  filter(timeseries >= 75)

Y3pred1 <- as.data.frame(fit, pars = "y3_mean_pred1") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred1 +1,
            "age_bins" = rep("Age bin#1: <8Wks", length(ts_pred1)))%>%
  filter(timeseries >= 75)

host_kiPred1 <- as.data.frame(fit, pars = "host_kiprop_pred1") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred1 +1,
            "age_bins" = rep("Age bin#1: <8Wks", length(ts_pred1)))%>%
  filter(timeseries >= 75)

Y4pred1 <- as.data.frame(fit, pars = "y4_mean_pred1") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred1 +1,
            "age_bins" = rep("Age bin#1: <8Wks", length(ts_pred1)))%>%
  filter(timeseries >= 75)

donor_kiPred2 <- as.data.frame(fit, pars = "donor_kiprop_pred2") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred2 -1,
            "age_bins" = rep("Age bin#2: 8-12Wks", length(ts_pred2)))%>%
  filter(timeseries >= 72)

Y3pred2 <- as.data.frame(fit, pars = "y3_mean_pred2") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred2 -1,
            "age_bins" = rep("Age bin#2: 8-12Wks", length(ts_pred2)))%>%
  filter(timeseries >= 72)

host_kiPred2 <- as.data.frame(fit, pars = "host_kiprop_pred2") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred2 -1,
            "age_bins" = rep("Age bin#2: 8-12Wks", length(ts_pred2)))%>%
  filter(timeseries >= 72)

Y4pred2 <- as.data.frame(fit, pars = "y4_mean_pred2") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred2 -1,
         "age_bins" = rep("Age bin#2: 8-12Wks", length(ts_pred2)))%>%
  filter(timeseries >= 72)

donor_kiPred3 <- as.data.frame(fit, pars = "donor_kiprop_pred3") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred3 - 8,
            "age_bins" = rep("Age bin#3: >12Wks", length(ts_pred3)))%>%
  filter(timeseries >= 113)

Y3pred3 <- as.data.frame(fit, pars = "y3_mean_pred3") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred3 - 8,
            "age_bins" = rep("Age bin#3: >12Wks", length(ts_pred3)))%>%
  filter(timeseries >= 113)

host_kiPred3 <- as.data.frame(fit, pars = "host_kiprop_pred3") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred3 - 8,
            "age_bins" = rep("Age bin#3: >12Wks", length(ts_pred3)))%>%
  filter(timeseries >= 113)

Y4pred3 <- as.data.frame(fit, pars = "y4_mean_pred3") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred3 - 8,
            "age_bins" = rep("Age bin#3: >12Wks", length(ts_pred3)))%>%
  filter(timeseries >= 113)

#combined pot for ki67
p3 <- ggplot() +
  geom_line(data = Y3pred1, aes(x = timeseries, y = median), col = "#356af2") +
  geom_line(data = Y4pred1, aes(x = timeseries, y = median), col = '#356af2', linetype = 2) +
  geom_line(data = Y3pred2, aes(x = timeseries, y = median), col = "#28af46") +
  geom_line(data = Y4pred2, aes(x = timeseries, y = median), col = '#28af46', linetype = 2) +
  geom_line(data = Y3pred3, aes(x = timeseries, y = median), col = "#ea2745") +
  geom_line(data = Y4pred3, aes(x = timeseries, y = median), col = '#ea2745', linetype = 2) +
  geom_point(data = ki67_binned, aes(x = age.at.S1K, y = prop_ki67hi, group = interaction(age_bins, subpopulation), color = age_bins, shape = subpopulation)) +
  scale_color_manual(values = c("#356af2", "#28af46", "#ea2745"), name = "", labels = c("<8Wks", "8-12Wks", "<12wks"))+
  scale_shape_manual(values = c(19, 1), name = "",  labels = c("Donor", "Host"))+
  labs(x = "Host age (days)", y = NULL, title = paste("Proportions of ki67Hi cells: ", substr(modelName, 10, 13))) +
  scale_x_continuous(limits = c(60, 600), trans = "log10",  breaks = c(10, 30, 100, 300))+
  scale_y_continuous(limits =c(0.0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) +  myTheme +
  theme(legend.position = c(0.75,0.85), legend.direction = "horizontal") + guides(color = FALSE)

#gathering data frames
donor_ki67_median_pooled <- rbind(Y3pred1, Y3pred2, Y3pred3)%>% select(-contains("key"))
donor_ki67_sd_pooled <- rbind(donor_kiPred1, donor_kiPred2, donor_kiPred3)%>% select(-contains("key"))
host_ki67_median_pooled <- rbind(Y4pred1, Y4pred2, Y4pred3)%>% select(-contains("key"))
host_ki67_sd_pooled <- rbind(host_kiPred1, host_kiPred2, host_kiPred3)%>% select(-contains("key"))


ki67_facet <- ggplot() +
  geom_line(data = donor_ki67_median_pooled, aes(x = timeseries, y = median), color = "#0000ba") +
  geom_ribbon(data = donor_ki67_median_pooled, aes(x = timeseries, ymin = lb, ymax = ub), fill = "#0000ba", alpha = 0.2)+
  geom_line(data = host_ki67_median_pooled, aes(x = timeseries, y = median), color = "#ed3e00") +
  geom_ribbon(data = host_ki67_median_pooled, aes(x = timeseries, ymin = lb, ymax = ub), fill = "#ed3e00", alpha = 0.2)+
  geom_point(data = ki67_binned, aes(x = age.at.S1K, y = prop_ki67hi, color = subpopulation)) +
  scale_alpha_manual(values = c(1, 0.7, 0.3), name = "Age bins", labels = c("<8Wks", "8-12Wks", "<12wks"))+
  scale_color_manual(values=c("#0000ba", "#ed3e00"), name = NULL, labels = c("Donor", "Host"))+ 
  guides(alpha = FALSE) + myTheme + theme(legend.position = c(0.92,0.75), legend.background = element_blank())


## saving  plots for quality control 
pdf(file = file.path(figDir, paste(modelName,"IndPlots%03d.pdf", sep = "")),
    width = 7, height = 8, onefile = FALSE )

# plots for individual age bins
nfd_plot <- nfd_facet +
  geom_ribbon(data = nfd_sd_pooled, aes(x = timeseries, ymin = lb, ymax = ub, fill = age_bins), alpha = 0.1) +
  labs(x = "Host age (days)", y = NULL, title = paste("Normalised Donor fractions: ", substr(modelName, 10, 13))) +
  scale_x_continuous(limits = c(40, 575), breaks = c(0,150,300,450, 600))+
  scale_y_continuous(limits =c(0, 1.15), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) + guides(color = FALSE) +
  facet_wrap(~ age_bins)

# plots for individual age bins
counts_plot <- counts_facet +
  geom_ribbon(data = counts_sd_pooled, aes(x = timeseries, ymin = lb, ymax = ub, fill = age_bins), alpha = 0.1) + guides(color = FALSE) +
  facet_wrap(~ age_bins) 

# plots for individual age bins
ki67_plot <- ki67_facet +
  labs(x = "Host age (days)", y = NULL, title = paste("Normalised Donor fractions: ", substr(modelName, 10, 13))) +
  scale_x_continuous(limits = c(60, 600), trans = "log10",  breaks = c(10, 30, 100, 300))+
  scale_y_continuous(limits =c(0, 1.15), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) + guides(color = FALSE) +
  facet_wrap(~ age_bins)

# plots for individual age bins
ki67_plot2 <- ki67_facet +
  geom_ribbon(data = donor_ki67_sd_pooled, aes(x = timeseries, ymin = lb, ymax = ub), fill = "#0000ba",  alpha = 0.15) +
  geom_ribbon(data = host_ki67_sd_pooled, aes(x = timeseries, ymin = lb, ymax = ub), fill = "#ed3e00", alpha = 0.15) +
  labs(x = "Host age (days)", y = NULL, title = paste("Normalised Donor fractions: ", substr(modelName, 10, 13))) +
  scale_x_continuous(limits = c(60, 600), trans = "log10",  breaks = c(10, 30, 100, 300))+
  scale_y_continuous(limits =c(0, 1.15), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) + guides(color = FALSE) +
  facet_wrap(~ age_bins)

gridExtra::grid.arrange(counts_plot, nfd_plot, ki67_plot,  nrow =  3)
gridExtra::grid.arrange(counts_plot, nfd_plot, ki67_plot2,  nrow =  3)

dev.off()


pdf(file = file.path(figDir, paste(modelName,"Plots%03d.pdf", sep = "")),
    width = 7, height = 5, onefile = FALSE )

lay1 <- rbind(c(1,1,2,2),
             c(NA,3,3,NA))

lay2 <- rbind(c(1,2),
              c(3,3))

gridExtra::grid.arrange(counts_facet, nfd_comb, p3, layout_matrix = lay1)
gridExtra::grid.arrange(counts_facet, nfd_comb, ki67_plot, layout_matrix = lay2)

dev.off()
