## clearing the environment
rm(list = ls())  
gc()    
setwd("/opt/mesh/eigg/sanket/FM_GC_ageCorrected")


library(rstan)
library(bayesplot)
library(tidyverse)
library(parallel)

## Stan-fit specific details

#modelName <- "source_counts"
modelName <- "Source_chi"
data_derived1 <- "T1_data.csv"
data_derived2 <- "T2_data.csv"
data_derived3 <- "Tra_data.csv"
data_derived4 <- "FM_data.csv"

## Relative paths assuming the working directory is the script directory
## containing this script
projectDir <- getwd()
scriptDir <- file.path(projectDir, "scripts")
modelDir <- file.path(projectDir, "models")
outputDir <- file.path(projectDir, "output")
saveDir <- file.path(outputDir, modelName)
figDir <- file.path(projectDir, "deliv", "figures", modelName)

## file names for the saved fit objects
output_filename1 = paste(modelName, "_", substr(data_derived1, 1,2), sep="")
output_filename2 = paste(modelName, "_", substr(data_derived2, 1,2), sep="")
output_filename3 = paste(modelName, "_", substr(data_derived3, 1,2), sep="")
output_filename4 = paste(modelName, "_", substr(data_derived4, 1,2), sep="")

## loading saved stanfit object
stanfit1 <- readRDS(file.path(saveDir, paste(output_filename1, ".rds", sep = "")))
stanfit2 <- readRDS(file.path(saveDir, paste(output_filename2, ".rds", sep = "")))
stanfit3 <- readRDS(file.path(saveDir, paste(output_filename3, ".rds", sep = "")))
stanfit4 <- readRDS(file.path(saveDir, paste(output_filename4, ".rds", sep = "")))


##################################################################################################
## import the data set
data_imp1 <- read_csv(file.path(projectDir, "data", data_derived1))
data_imp2 <- read_csv(file.path(projectDir, "data", data_derived2))
data_imp3 <- read_csv(file.path(projectDir, "data", data_derived3))
data_imp4 <- read_csv(file.path(projectDir, "data", data_derived4))

# time points in the data
ts_pred <- seq(0, 600) + 72

## Specify the variables for which you want history and density plots
parametersToPlot <- c("chiEst","qEst","sigma")

## Additional variables to monitor
otherRVs <- c("ymean_pred","chipred")

parameters <- c(parametersToPlot, otherRVs)

################################################################################################
## posterior distributions of parameters
### plotting options
options(bayesplot.base_size = 15,
        bayesplot.base_family = "sans")
myTheme <- theme(text = element_text(size = 11), axis.text = element_text(size = 11), axis.title =  element_text(size = 11, face = "bold"),
                 plot.title = element_text(size=11,  hjust = 0.5, face = "bold"),
                 legend.background = element_blank(), legend.key = element_blank(),
                 legend.text = element_text(size= 9), legend.title = element_text(9))

# setting ggplot theme for rest fo the plots
theme_set(theme_bw())

####### plotting
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

################################################################################################
## posterior predictive distributions
Ypred1 <- as.data.frame(stanfit1, pars = "ymean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred) %>%  filter(timeseries >=40) %>%  filter(timeseries <= 575)

Ypred2 <- as.data.frame(stanfit2, pars = "ymean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred) %>%  filter(timeseries >=40) %>%  filter(timeseries <= 575)

Ypred3 <- as.data.frame(stanfit3, pars = "ymean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred) %>%  filter(timeseries >=40) %>%  filter(timeseries <= 575)


Ypred4 <- as.data.frame(stanfit4, pars = "ymean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred) %>%  filter(timeseries >=40) %>%  filter(timeseries <= 575)

## plots for counts
#P1 <- ggplot() +
#  geom_hline(aes(yintercept = 1), linetype = 2)+
#  #geom_ribbon(data = fdpred, aes(x = timeseries, ymin = lb, ymax=ub), fill = "#FF0000", alpha = 0.2) +
#  geom_ribbon(data = Ypred1, aes(x = timeseries, ymin = lb, ymax = ub), fill = "#5BBCD6", alpha = 0.2)+
#  geom_line(data = Ypred1, aes(x = timeseries, y = median), col = "#5BBCD6") +
#  geom_point(data = data_imp1, aes(x = age.at.S1K, y = total_counts),  col = "#5BBCD6") +
#  labs(title=paste('Total numbers of', substr(data_derived1, 1, 2), "B cells"),  y=NULL, x="Host age (days)") + 
#  scale_x_continuous(limits = c(60, 600), trans = "log10", breaks = c(75, 150, 300, 600))+
#  scale_y_continuous(limits = c(1e5, 1e7), trans="log10", breaks=c(1e4, 1e5, 1e6, 1e7, 1e8), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
#  guides(color = FALSE)+ myTheme
#
#P2 <- ggplot() +
#  geom_hline(aes(yintercept = 1), linetype = 2)+
#  #geom_ribbon(data = fdpred, aes(x = timeseries, ymin = lb, ymax=ub), fill = "#F2AD00", alpha = 0.2) +
#  geom_ribbon(data = Ypred2, aes(x = timeseries, ymin = lb, ymax = ub), fill = "#F2AD00", alpha = 0.2)+
#  geom_line(data = Ypred2, aes(x = timeseries, y = median), col = "#F2AD00") +
#  geom_point(data = data_imp2, aes(x = age.at.S1K, y = total_counts),  col = "#F2AD00") +
#  labs(title=paste('Total numbers of', substr(data_derived2, 1, 2), "B cells"),  y=NULL, x="Host age (days)") + 
#  scale_x_continuous(limits = c(60, 600), trans = "log10", breaks = c(75, 150, 300, 600))+
#  scale_y_continuous(limits = c(1e5, 1e7), trans="log10", breaks=c(1e4, 1e5, 1e6, 1e7, 1e8), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
#  guides(color = FALSE)+ myTheme
#
#P3 <- ggplot() +
#  geom_hline(aes(yintercept = 1), linetype = 2)+
#  #geom_ribbon(data = fdpred, aes(x = timeseries, ymin = lb, ymax=ub), fill = "#F98400", alpha = 0.2) +
#  geom_ribbon(data = Ypred3, aes(x = timeseries, ymin = lb, ymax = ub), fill = "#00A08A", alpha = 0.2)+
#  geom_line(data = Ypred3, aes(x = timeseries, y = median), col = "#00A08A") +
#  geom_point(data = data_imp3, aes(x = age.at.S1K, y = total_counts),  col = "#00A08A") +
#  labs(title=paste("Total numbers of T1+T2 B cells"),  y=NULL, x="Host age (days)") + 
#  scale_x_continuous(limits = c(60, 600), trans = "log10", breaks = c(75, 150, 300, 600))+
#  scale_y_continuous(limits = c(1e5, 1e7), trans="log10", breaks=c(1e4, 1e5, 1e6, 1e7, 1e8), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
#  guides(color = FALSE)+ myTheme
#
#P4 <- ggplot() +
#  geom_hline(aes(yintercept = 1), linetype = 2)+
#  #geom_ribbon(data = fdpred, aes(x = timeseries, ymin = lb, ymax=ub), fill = "#5BBCD6", alpha = 0.2) +
#  geom_ribbon(data = Ypred4, aes(x = timeseries, ymin = lb, ymax = ub), fill = "#F98400", alpha = 0.2)+
#  geom_line(data = Ypred4, aes(x = timeseries, y = median), col = "#F98400") +
#  geom_point(data = data_imp4, aes(x = age.at.S1K, y = total_counts),  col = "#F98400") +
#  labs(title=paste('Total numbers of', substr(data_derived4, 1, 2), "B cells"),  y=NULL, x="Host age (days)") + 
#  scale_x_continuous(limits = c(60, 600), trans = "log10", breaks = c(75, 150, 300, 600))+
#  scale_y_continuous(limits = c(1e6, 1.2e8), trans="log10", breaks=c(1e4, 1e5, 1e6, 1e7, 1e8), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
#  guides(color = FALSE)+ myTheme

## plots for chimerism
P1 <- ggplot() +
  geom_hline(aes(yintercept = 1), linetype = 2)+
  #geom_ribbon(data = fdpred, aes(x = timeseries, ymin = lb, ymax=ub), fill = "#FF0000", alpha = 0.2) +
  geom_ribbon(data = Ypred1, aes(x = timeseries, ymin = lb, ymax = ub), fill = "#5BBCD6", alpha = 0.2)+
  geom_line(data = Ypred1, aes(x = timeseries, y = median), col = "#5BBCD6") +
  geom_point(data = data_imp1, aes(x = age.at.S1K, y = fd),  col = "#5BBCD6") +
  labs(x = "Host age (days)", y = NULL, title = paste("Chimerism within ", substr(data_derived1, 1,2), " B cells",  sep="")) +
  scale_x_continuous(limits = c(40, 575), breaks = c(50,150, 250, 350, 450, 550))+
  scale_y_continuous(limits =c(0, 1.05), breaks = c(0, 0.25, 0.5, 0.75, 1.0))+ 
  guides(color = FALSE)+ myTheme

P2 <- ggplot() +
  geom_hline(aes(yintercept = 1), linetype = 2)+
  #geom_ribbon(data = fdpred, aes(x = timeseries, ymin = lb, ymax=ub), fill = "#F2AD00", alpha = 0.2) +
  geom_ribbon(data = Ypred2, aes(x = timeseries, ymin = lb, ymax = ub), fill = "#F2AD00", alpha = 0.2)+
  geom_line(data = Ypred2, aes(x = timeseries, y = median), col = "#F2AD00") +
  geom_point(data = data_imp2, aes(x = age.at.S1K, y = fd),  col = "#F2AD00") +
  labs(x = "Host age (days)", y = NULL, title = paste("Chimerism within ", substr(data_derived2, 1,2), " B cells",  sep="")) +
  scale_x_continuous(limits = c(40, 575), breaks = c(50,150, 250, 350, 450, 550))+
  scale_y_continuous(limits =c(0, 1.05), breaks = c(0, 0.25, 0.5, 0.75, 1.0))+ 
  guides(color = FALSE)+ myTheme

P3 <- ggplot() +
  geom_hline(aes(yintercept = 1), linetype = 2)+
  #geom_ribbon(data = fdpred, aes(x = timeseries, ymin = lb, ymax=ub), fill = "#F98400", alpha = 0.2) +
  geom_ribbon(data = Ypred3, aes(x = timeseries, ymin = lb, ymax = ub), fill = "#00A08A", alpha = 0.2)+
  geom_line(data = Ypred3, aes(x = timeseries, y = median), col = "#00A08A") +
  geom_point(data = data_imp3, aes(x = age.at.S1K, y = fd),  col = "#00A08A") +
  labs(x = "Host age (days)", y = NULL, title = paste("Chimerism within ", substr(data_derived1, 1,2), "+", substr(data_derived2, 1,2), " B cells",  sep="")) +
  scale_x_continuous(limits = c(40, 575), breaks = c(50, 150, 250, 350, 450, 550))+
  scale_y_continuous(limits =c(0, 1.05), breaks = c(0, 0.25, 0.5, 0.75, 1.0))+ 
  guides(color = FALSE)+ myTheme

P4 <- ggplot() +
  geom_hline(aes(yintercept = 1), linetype = 2)+
  #geom_ribbon(data = fdpred, aes(x = timeseries, ymin = lb, ymax=ub), fill = "#5BBCD6", alpha = 0.2) +
  geom_ribbon(data = Ypred4, aes(x = timeseries, ymin = lb, ymax = ub), fill = "#F98400", alpha = 0.2)+
  geom_line(data = Ypred4, aes(x = timeseries, y = median), col = "#F98400") +
  geom_point(data = data_imp4, aes(x = age.at.S1K, y = fd),  col = "#F98400") +
  labs(x = "Host age (days)", y = NULL, title = paste("Chimerism within ", substr(data_derived4, 1,2), " B cells",  sep="")) +
  scale_x_continuous(limits = c(40, 575), breaks = c(50, 150, 250, 350, 450, 550))+
  scale_y_continuous(limits =c(0, 1.05), breaks = c(0, 0.25, 0.5, 0.75, 1.0))+ 
  guides(color = FALSE)+ myTheme

if (!file.exists(figDir)){
  dir.create(figDir)
}

## saving  plots for quality control 
pdf(file = file.path(figDir, paste(modelName,"SplinePlots%03d.pdf", sep = "")),
    width = 7, height = 5, onefile = FALSE, useDingbats = FALSE)

cowplot::plot_grid(P1, P2, P3, P4,  labels = c("A", "B", "C", "D"), nrow = 2)


dev.off()
