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
#data_derived4 <- "FM_data.csv"

## Relative paths assuming the working directory is the script directory
## containing this script
projectDir <- getwd()
scriptDir <- file.path(projectDir, "scripts")
modelDir <- file.path(projectDir, "models")
outputDir <- file.path(projectDir, "output")
saveDir <- file.path(outputDir, modelName)
figDir <- file.path(projectDir, "deliv", "figures", modelName)

## Specify the variables for which you want history and density plots
parametersToPlot <- c("chiEst","qEst","sigma")

## file names for the saved fit objects
output_filename1 = paste(modelName, "_", substr(data_derived1, 1,2), sep="")
output_filename2 = paste(modelName, "_", substr(data_derived2, 1,2), sep="")
output_filename3 = paste(modelName, "_", substr(data_derived3, 1,2), sep="")
#output_filename4 = paste(modelName, "_", substr(data_derived4, 1,2), sep="")

## loading saved stanfit object
stanfit1 <- readRDS(file.path(saveDir, paste(output_filename1, ".rds", sep = "")))
stanfit2 <- readRDS(file.path(saveDir, paste(output_filename2, ".rds", sep = "")))
stanfit3 <- readRDS(file.path(saveDir, paste(output_filename3, ".rds", sep = "")))
#stanfit4 <- readRDS(file.path(saveDir, paste(output_filename4, ".rds", sep = "")))

### posterior distributions of parameters
ptable1 <- monitor(as.array(stanfit1, pars = parametersToPlot), warmup = 0, print = FALSE)
ptable1[,c(1,4,8)]
### posterior distributions of parameters
ptable2 <- monitor(as.array(stanfit2, pars = parametersToPlot), warmup = 0, print = FALSE)
ptable2[,c(1,4,8)]
### posterior distributions of parameters
ptable3 <- monitor(as.array(stanfit3, pars = parametersToPlot), warmup = 0, print = FALSE)
ptable3[,c(1,4,8)]


