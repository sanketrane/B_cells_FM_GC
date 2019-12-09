## log linear model -- for counts as spline

## clearing the environment
rm(list = ls())  
gc()    
setwd("~/Desktop/GIt_repos/FM_GC_ageCorrected")


library(rstan)
library(tidyverse)
library(parallel)

modelName <- "T1_Ontogeny_counts"
data_derived <- "T1_spleen.csv"
#data_derived <- "T2_pooled.csv"


## Relative paths assuming the working directory is the script directory
## containing this script
projectDir <- getwd()
scriptDir <- file.path(projectDir, "scripts")
figDir <- file.path(projectDir, "deliv", "figures", modelName)
tabDir <- file.path(projectDir, "deliv", "tables", modelName)
dataDir <- file.path(projectDir, "data", data_derived)
modelDir <- file.path(projectDir, "models")
outputDir <- file.path(projectDir, "output")
saveDir <- file.path(outputDir, modelName)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
# set different seeds for different tasks
i <- as.numeric(Sys.getenv("SLURM_PROCID"))

seed <- 1256 + i
################################################################################################

## import the data set
data_imp <- read_csv(dataDir) %>% arrange(days_post_t0) #%>% filter(set_name == "young") 

ts_pred <- seq(0, 600, 0.5)

## create data set
data <- with(data_imp,
             list(
               numObs = nrow(data_imp),
               Time = data_imp$days_post_t0,
               counts = data_imp$total_counts,
               numPred = length(ts_pred),
               ts_pred = ts_pred
               ))

## create initial estimates
init <- function() list(
  theta0_Log = rnorm(1, 18, 0.1),
  r1 = exp(rnorm(1, log(0.3), 0.01)),
  n1 = exp(rnorm(1, log(2), 0.01)),
  #r3 = exp(rnorm(1, log(0.001), 0.01)),
  #alpha = exp(rnorm(1, log(0.5), 0.01)),
    
  sigma = exp(rnorm(1,log(1.5), 1)))

## Specify the variables for which you want history and density plots
parametersToPlot <- c("theta0_Log", "r1", "sigma", "n1")#,  "r3", "alpha")#  

## Additional variables to monitor
otherRVs <- c("ymean_pred", "countspred")

parameters <- c(parametersToPlot, otherRVs)

################################################################################################
# run Stan

nChains <- 6
nPost <- 1500 ## Number of post-burn-in samples per chain after thinning
nBurn <- 500 ## Number of burn-in samples per chain after thinning
nThin <- 1

nIter <- (nPost + nBurn) * nThin
nBurnin <- nBurn * nThin

fit <- stan(file = file.path(modelDir, paste(modelName, ".stan", sep = "")),
            data = data,
            pars = parameters,
            iter = nIter,
            warmup = nBurnin,
            thin = nThin, 
            init = init,
            control = list(adapt_delta = 0.9),
            chains = nChains)

# save results in output directory
if (!file.exists(saveDir)){
  dir.create(saveDir)
}

# saving output file as 
output_filename = paste(modelName, "_", substr(data_derived, 1,2), sep="")
saveRDS(fit, file = file.path(saveDir, paste(output_filename, ".rds", sep = "")))

################################################################################################
################################################################################################
### posterior distributions of parameters
ptable <- monitor(as.array(fit, pars = parametersToPlot), warmup = 0, print = FALSE)
ptable[,c(1,4,8)]

## posterior predictive distributions
options(bayesplot.base_size = 15,
        bayesplot.base_family = "sans")
color_scheme_set(scheme = "viridis")
my_theme <- theme(axis.text = element_text(size = 12),
                  axis.title =  element_text(size = 12, face = "bold"),
                  plot.title = element_text(size= 12,  hjust = 0.5, face = "bold"),
                  legend.text = element_text(size= 12), legend.background = element_blank(),
                  legend.title = element_text(size = 12, face = "bold"), legend.position = c(0.9, 0.35))

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

# time sequnce used for plotting 
ts_pred = seq(from = 0, to = 600, by = 0.5)

# Total cell counts
Cpred <- as.data.frame(fit, pars = "countspred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred + 10)

Ypred <- as.data.frame(fit, pars = "ymean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred + 10)

Ypred <- Ypred%>%
  filter(timeseries >= 10) %>% filter(timeseries <= 610)

Cpred <- Cpred%>%
  filter(timeseries >= 10)%>% filter(timeseries <= 610)

pdf(file = file.path("deliv/figures/FM_ontogeny_sim", paste("T1_full", "_exp", ".pdf", sep = "")),
    width = 6, height = 4, onefile = F, useDingbats = FALSE)

ggplot() +
  #geom_ribbon(data = Cpred, aes(x = timeseries, ymin = lb, ymax = ub), fill = "#E1AF00",  alpha = 0.2)+
  geom_ribbon(data = Ypred, aes(x = timeseries, ymin = lb, ymax = ub), fill = "#E1AF00", alpha = 0.4)+
  geom_line(data = Ypred, aes(x = timeseries, y = median), col = "#3B9AB2", size=1) +
  geom_point(data = data_imp, aes(x = age.at.S1K, y = total_counts, colour = set_name)) +
  scale_color_manual(values = c("#F21A00", "#3B9AB2"), name = "Data", labels = c("Chimera", "Ontogeny")) +
  labs(title=paste('T1 counts in spleen'),  y=NULL, x= "Host age (days)") + 
  scale_x_continuous(limits = c(10, 610), trans = "log10", breaks=c(10, 30, 100, 300))+
  scale_y_continuous(limits = c(1e5, 5e7), trans="log10", breaks=c(1e4, 1e5, 1e6, 1e7, 1e8), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  theme(axis.text = element_text(size = 11),
        axis.title =  element_text(size = 11, face = "bold"),
        plot.title = element_text(size=11,  hjust = 0.5, face = "bold"),
        legend.text = element_text(size=11), legend.background = element_blank(),
        legend.title = element_text(size = 11, face = "bold"), legend.position = c(0.85, 0.8))

dev.off()
