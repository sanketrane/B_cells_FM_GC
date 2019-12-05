## Simple linear ODE model -- Homogeneous model 
## clearing the environment
rm(list = ls())  
gc()    
setwd("/opt/mesh/eigg/sanket/FM_GC_ageCorrected")
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.txt"
}

####################################################################################
## Installing r-stan pachage on the go:
if(!("rstan" %in% rownames(installed.packages())) ){
  install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies=TRUE)
}

## Installing loo:
if(!("loo" %in% rownames(installed.packages())) ){
  install.packages("loo", repos = "https://cloud.r-project.org/", dependencies=TRUE)
}

## Installing tidyverse:
if(!("tidyverse" %in% rownames(installed.packages())) ){
  install.packages("tidyverse", repos = "https://cloud.r-project.org/", dependencies=TRUE)
}
####################################################################################

# loading libararies
require(rstan)
require(parallel)
require(loo)
require(tidyverse)

# model sepcific details
modelName <- "ki67_SVM_SPGC"        # name of the file for stan model
data_derived1 <- paste(args[1], "_data.csv", sep = "")            # name of the file for precursor pop
data_derived2 <- "counts_SPGC.csv"  # name of the file for counts data of traget pop
data_derived3 <- "Nfd_SPGC.csv"     # name of the file for normalised fd data of target pop
data_derived4 <- "ki67_SPGC.csv"    # name of the file for ki67 data of  pop


## Setting all the directories for opeartions
projectDir <- getwd()
scriptDir <- file.path(projectDir, "scripts", "YFP_informed_models", "SPGC")
modelDir <- file.path(projectDir, "models", "YFP_informed_models", "SPGC")
dataDir <- file.path(projectDir, "data")
toolsDir <- file.path(scriptDir, "tools")
outputDir <- file.path(projectDir, "output", "YFP_informed_models", "SPGC")
saveDir <- file.path(outputDir, paste(modelName, "_", substr(data_derived1, 1,2), sep=""))  ### varies according to the precursor pop used

# parallelisation for R-Stan
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# set different seeds for different tasks
### if the script is called on different nodes across the cluster then it will use different seed so as to optimise the sampling 
i <- as.numeric(Sys.getenv("SLURM_PROCID"))
seed <- 4009 + i                     ### Use your lucky mumber! haha! 

################################################################################################
################################################################################################
# loading data files
source_sorted <- read_csv(file.path(dataDir, data_derived1))%>% arrange(age.at.S1K)
counts_sorted <- read_csv(file.path(dataDir, data_derived2))%>% arrange(age.at.S1K)
Nfd_sorted <- read_csv(file.path(dataDir, data_derived3))%>% arrange(age.at.S1K)
ki67_sorted <- read_csv(file.path(dataDir, data_derived4))%>% arrange(age.at.S1K)

### solving ODEs only for unique timepoints 
### also Stan ODE solver throws an error for repetitive time-points 
# unique time points in data for odes solver
unique_times_df <- counts_sorted %>% distinct(age.at.S1K, .keep_all = TRUE)

data_time <- counts_sorted$age.at.S1K                                          # data and solver time 
solve_time <- unique_times_df$age.at.S1K                                       # unique time points to solve ode
time_index <- purrr::map_dbl(data_time, function(x) which(x == solve_time))    # keeping track of index of time point in relation to solve_time


# delay time for age correction 
# for each host with different age at BMT -- time zero is different which is becomes another variable in the ODEs -- is it delayed differential equation now?
# tb_time -- for solving for initial conditions at each tb --  we use the earlist host age at BMT to calculate the initial conditions at each tb
solve_ageAtBMT <- unique_times_df$age.at.bmt   
tb_time <- solve_ageAtBMT %>% unique() %>% sort()                              # unique age at BMT points to solve ode
tb_index <- purrr::map_dbl(solve_ageAtBMT, function(x) which(x == tb_time))    #keeping track of index of  delay time in relation to tb_time


# time sequence for predictions specific to age bins within the data
# In this analysis 3 age bins were selected with mean ages for each bin calculated as follows:
counts_agebmt_sorted <- counts_sorted %>% arrange(age.at.bmt)
Mean_agebin1 <- round(mean(counts_agebmt_sorted$age.at.bmt[1:12]))
Mean_agebin2 <- round(mean(counts_agebmt_sorted$age.at.bmt[13:19]))
Mean_agebin3 <- round(mean(counts_agebmt_sorted$age.at.bmt[20:24]))

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

# setting data input according to the precurosr pop used
## these values dictate how splines for counts and chimerism change with time.
## obtained by fitting splines separately to counts and chimerism in the respective precurosr pop
if (grepl("T2", data_derived1) == TRUE){
  source_list = data.frame("nu" = -5.53e-05, "theta0" = 14.13, "chiEst" = 0.73, "qEst" = 0.06, "eps_donor" = 0.82, "eps_host" = 0.73)
} else if (grepl("T1", data_derived1) == TRUE) {
  source_list = data.frame("nu" = 1e-04, "theta0" = 13.81, "chiEst" = 0.78, "qEst" = 0.07, "eps_donor" = 0.97, "eps_host" = 0.94)
} else if (grepl("FM", data_derived1) == TRUE){
  source_list = data.frame("nu" = -1.1e-03 , "theta0" = 17.20, "chiEst" = 0.77, "qEst" = 0.023, "eps_donor" = 0.31, "eps_host" = 0.14)
} else {
  source_list = data.frame("nu" = 3.30e-04, "theta0" = 14.83, "chiEst" = 0.75, "qEst" = 0.066, "eps_donor" = 0.88, "eps_host" = 0.81)
}

## data list that feeds into the Stan model
## Use names that exactly match with the Stan file
data <- list(
  numObs = nrow(counts_sorted),              # Total number of observations within data
  Nd_0 = 0,                                  # donor counts at t0
  solve_time = solve_time,                   # unique observations for ode solver
  num_index = length(solve_time),            # to declare dimensions of the array in Stan
  num_tb = length(tb_time),                 # to declare dimensions of the array in Stan
  time_index = time_index,
  tb_index = tb_index,
  ageAtBMT = solve_ageAtBMT,
  tb_time = tb_time,
  dpBMT = counts_sorted$days.post.bmt,       #to be used for the spline of chimerism of the precurosr pop
  counts = counts_sorted$total_counts,       ## data to fit
  Nfd = Nfd_sorted$Nfd,                      ## data to fit
  ki_host = ki67_sorted$host_ki67_GC,        ## data to fit
  ki_donor = ki67_sorted$donor_ki67_GC,      ## data to fit
  numPred1 = length(ts_pred1),                # to declare dimensions of the array in Stan
  numPred2 = length(ts_pred2),                # to declare dimensions of the array in Stan
  numPred3 = length(ts_pred3),                # to declare dimensions of the array in Stan
  ts_pred1 = ts_pred1,
  ts_pred2 = ts_pred2,
  ts_pred3 = ts_pred3,
  tb_pred1 = tb_pred1,
  tb_pred2 = tb_pred2,
  tb_pred3 = tb_pred3,
  tb_time_pred1 = tb_time_pred1,
  tb_time_pred2 = tb_time_pred2,
  tb_time_pred3 =tb_time_pred3,
  theta0 = exp(source_list$theta0),        #size of precurosr pop ar t0
  nu = source_list$nu,                     #rate of change of precurosr pop with time
  chiEst = source_list$chiEst,             #stabilised value of chimerism for precurosr pop
  qEst = source_list$qEst,                 #rate of change of chimerism in precurosr pop
  eps_donor = source_list$eps_donor,       #proportion of ki67+ cells in donor precurosr pop
  eps_host = source_list$eps_host          #proportion of ki67+ cells in host precurosr pop
)


## create initial estimates
init <- function() list(
  p_Influx = exp(rnorm(1, log(0.1), 0.1)),
  rho = exp(rnorm(1, log(0.1), 0.01)),
  deltaYFP = exp(rnorm(1, log(0.4), 0.01)),
  r_psi = exp(rnorm(1, log(0.001), 0.001)),
  
  y0_Log = rnorm(1, 11 , 0.1),
  kappa_0 = exp(rnorm(1, log(0.9), 0.01)),
  Beta =  exp(rnorm(1, log(3.5), 0.1)),
  
  
  y0_Log = rnorm(1, 11 , 0.5),
  kappa_0 = exp(rnorm(1, log(0.9), 0.01)),
  Beta =  exp(rnorm(1, log(3.5), 0.01)),
  
  sigma1 = exp(rnorm(1,log(1.5), 1)),
  sigma2 = exp(rnorm(1,log(1.5), 1)),
  sigma3 = exp(rnorm(1,log(1.5), 1)),
  sigma4 = exp(rnorm(1,log(1.5), 1))
)


## Specify the variables for which you want history and density plots
parametersToPlot <- c("p_Influx", "psi", "rho", "r_psi", "deltaYFP", "delta", "Beta", "y0_Log", "kappa_0", 
                      "sigma1", "sigma2", "sigma3", "sigma4")

## Additional variables to monitor
otherRVs <- c("y1_mean_pred_age1", "countspred_age1", "y1_mean_pred_age2", "countspred_age2", "y1_mean_pred_age3", "countspred_age3",
              "y2_mean_pred_age1", "fdpred_age1", "y2_mean_pred_age2", "fdpred_age2", "y2_mean_pred_age3", "fdpred_age3", 
              "y3_mean_pred1",  "donor_kiprop_pred1", "y4_mean_pred1",  "host_kiprop_pred1",
              "y3_mean_pred2",  "donor_kiprop_pred2", "y4_mean_pred2",  "host_kiprop_pred2",
              "y3_mean_pred3",  "donor_kiprop_pred3", "y4_mean_pred3",  "host_kiprop_pred3", "log_lik",
              "lambda_inv", "rho_inv", "lambda_calc", "delta_inv", "r_psi_inv",
              "log_lik1", "log_lik2", "log_lik3", "log_lik4")

parameters <- c(parametersToPlot, otherRVs)
################################################################################################
## run Stan

# parameters for running fits
nChains <- 2
nPost <- 500 ## Number of post-burn-in samples per chain after thinning
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

################################################################################################
# save results in output directory
if (!file.exists(saveDir)){
  dir.create(saveDir)
}

### distinguishing models runs from each other and gathering information about individual tasks

# saving output file as 
output_filename=paste(Sys.getenv("SLURM_JOB_ID"), "_P", Sys.getenv("SLURM_PROCID"), "_",  modelName, ".rds",  sep="")

# saving the stan fit object for individual runs as 
saveRDS(fit, file = file.path(saveDir, output_filename))

## calculating an output from individual runs for the validation of a sucessful run
loo_loglik <- extract_log_lik(fit, parameter_name = "log_lik", merge_chains = TRUE)
loo_ic <- loo(loo_loglik)

## this value is added to the job specific run_info file which is gets written in the save directory after each job copmpletion

# information about the cluster run
cluster_run_info <- paste("Time = ", timestamp(),
                          "|| Job ID = ", Sys.getenv("SLURM_JOB_ID"),
                          "|| Modelname = ", modelName,
                          "|| crude_loo_ic = ", loo_ic$estimates[3],
                          "|| Host = ", Sys.getenv("SLURM_SUBMIT_HOST"),
                          "|| Nodes = ",  Sys.info()["nodename"],
                          "|| Process_ID = P", Sys.getenv("SLURM_PROCID"))

# writing the cluster run info as a data file in the out dir 
write.table(cluster_run_info, file = file.path(saveDir, paste("job_", Sys.getenv("SLURM_JOB_ID"), "_run_info.txt")), append = TRUE)
