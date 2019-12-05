## log linear model -- for counts as spline

## clearing the environment
rm(list = ls())  
gc()    
setwd("/opt/mesh/eigg/sanket/FM_GC_ageCorrected")


library(rstan)
library(bayesplot)
library(tidyverse)
library(parallel)

modelName <- "source_counts"
data_derived <- "T1_data.csv"
## Relative paths assuming the working directory is the script directory
## containing this script
projectDir <- getwd()
scriptDir <- file.path(projectDir, "scripts")
figDir <- file.path(projectDir, "deliv", "figures", modelName)
tabDir <- file.path(projectDir, "deliv", "tables", modelName)
dataDir <- file.path(projectDir, "data", data_derived)
modelDir <- file.path(projectDir, "models", "Miscellaneous")
outputDir <- file.path(projectDir, "output")
saveDir <- file.path(outputDir, modelName)

rstan_options(auto_write = TRUE)
#options(mc.cores = parallel::detectCores())
# set different seeds for different tasks
i <- as.numeric(Sys.getenv("SLURM_PROCID"))

seed <- 1256 + i
################################################################################################

## import the data set
data_imp <- read_csv(dataDir)

ts_pred <- seq(0, 600)

## create data set
data <- with(data_imp,
             list(
               numObs = nrow(data_imp),
               Time = data_imp$age.at.S1K,
               counts = data_imp$total_counts,
               numPred = length(ts_pred),
               ts_pred = ts_pred
               ))

## create initial estimates
init <- function() list(
    theta0 = rnorm(1, 15, 1),
    nu = rnorm(1, 0.001, 0.2),
    sigma = exp(rnorm(1,log(1.5), 1)))

## Specify the variables for which you want history and density plots
parametersToPlot <- c("theta0","nu","sigma")

## Additional variables to monitor
otherRVs <- c("ymean_pred", "countspred")

parameters <- c(parametersToPlot, otherRVs)

################################################################################################
# run Stan

nChains <- 12
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
