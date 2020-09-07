## Varying net Loss rate with host age in neonates

## clearing the environment
rm(list = ls())  
gc()    
#setwd("Desktop/GIt_repos/FM_GC_ageCorrected")


library(rstan)
library(tidyverse)
library(parallel)
library(loo)

modelName <- "FM_Ontogeny_fit_lambda"
data_derived <- "FM_Ontogeny.csv"
## Relative paths assuming the working directory is the script directory
## containing this script
projectDir <- getwd()
scriptDir <- file.path(projectDir, "scripts")
figDir <- file.path(projectDir, "deliv", "figures", paste(modelName, "_T1", sep = ""))
tabDir <- file.path(projectDir, "deliv", "tables", paste(modelName, "_T1", sep = ""))
dataDir <- file.path(projectDir, "data", data_derived)
modelDir <- file.path(projectDir, "models")
outputDir <- file.path(projectDir, "output")
saveDir <- file.path(outputDir, paste(modelName, "_T1", sep = ""))
toolsDir <- file.path(scriptDir, "tools")

# loadiong the scr# loadiong the script that contains functions for plotting stan parameters
source(file.path(toolsDir, "stanTools.R"))               

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
# set different seeds for different tasks
i <- as.numeric(Sys.getenv("SLURM_PROCID"))

seed <- 1256 + i
################################################################################################

## import the data set
data_imp <- read_csv(dataDir)%>% 
  mutate(days_post_t0 = age.at.S1K - 14)%>% arrange(days_post_t0)

#### FM counts form Ontogeny data
### reding and wrangling the excel file
### counts in spleen
FM_spl <- readxl::read_excel(path = "data/ontogeny_data.xlsx", sheet = 1) %>% 
  select(contains("ID"), contains("age"), starts_with("FM")) %>% arrange(age.at.S1K)

### counts in LN
FM_LN <- readxl::read_excel(path = "data/ontogeny_data.xlsx", sheet = 2) %>% 
  select(contains("ID"), contains("age"), starts_with("FM")) %>% arrange(age.at.S1K)

### pooled counts
FM_young <- right_join(FM_spl, FM_LN, by = c("mouse_ID", "age.at.S1K"), suffix = c("_spl", "_ln"))%>%
  na.omit() %>% mutate(total_counts = FM_spl + FM_ln,
                       days_post_t0 = age.at.S1K - 15,
                       set_name ="Ontogeny") %>%
  select(contains("age"), contains("days"), contains("set"), contains("total"))

## FM counts from chimera data
FM_adults <- read.csv("data/counts_FM.csv") %>% 
  mutate(days_post_t0 = age.at.S1K - 15,
         set_name = "Chimera")%>%
  select(contains("set_name"), contains("t0"), contains("S1K"), contains("total_counts"))

### pooling ontogeny and chimera data together
FM_counts <- rbind(FM_young, FM_adults) %>% arrange(set_name)


# time points in data
data_time <- data_imp$age.at.S1K
solve_time <- c(unique(data_time))   # unique time points to solve ode

#keep track of index of time point in relation to solve_time
time_index <- purrr::map_dbl(data_time, function(x) which(x == solve_time))

# time sequnce used for plotting 
ts_pred = seq(from = 14, to = 200, 0.1)

## create data set
data <- list(
  numObs = nrow(data_imp),
  solve_time = solve_time,
  num_index = length(solve_time),
  time_index = time_index,
  counts = data_imp$total_counts,
  numPred = length(ts_pred),
  ts_pred = ts_pred,
  theta0 = exp(13.9591976),
  r1 = 0.3485853,
  n1 = 2.1701241,
  eps = 0.97,
  psi =  0.7239,
  rho = 0.0059901,
  delta0 = 0.02761,
  Beta = 5.7562,
  y0_Log = 14.08421,
  kappa0 = 0.1352286)

## create initial estimates
init <- function() list(
  r_neo = exp(rnorm(1, log(0.05), 0.01)),
  #y0_Log = rnorm(1, 8, 0.01),
  
  sigma = exp(rnorm(1,log(1.5), 1)))

## Specify the variables for which you want history and density plots
parametersToPlot <- c("r_neo",  "sigma")

## Additional variables to monitor
otherRVs <- c("ymean_pred", "countspred", "log_lik")

parameters <- c(parametersToPlot, otherRVs)

################################################################################################
# run Stan

nChains <- 10
nPost <- 2000 ## Number of post-burn-in samples per chain after thinning
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
            #control = list(adapt_delta = 0.92),
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

## reading the FM ontongeny data fit
#output_filename = paste(modelName, "_", substr(data_derived, 1,2), sep="")
#fit <- readRDS(file = file.path(outputDir, modelName, paste(output_filename, ".rds", sep = "")))

## posterior predictive distributions
options(bayesplot.base_size = 15,
        bayesplot.base_family = "sans")
#color_scheme_set(scheme = "viridis")

# setting ggplot theme for rest fo the plots
theme_set(theme_bw())
my_theme <- theme(axis.text = element_text(size = 12),
                  axis.title =  element_text(size = 12, face = "bold"),
                  plot.title = element_text(size= 12,  hjust = 0.5, face = "bold"),
                  legend.text = element_text(size= 12), legend.background = element_blank(),
                  legend.title = element_text(size = 11, face = "bold"), legend.position = c(0.85, 0.2))

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

#output_filename = paste(modelName, "_", substr(data_derived, 1,2), sep="")
#fit <- readRDS(file = file.path(saveDir, paste(output_filename, ".rds", sep = "")))

# Total cell counts
Cpred <- as.data.frame(fit, pars = "countspred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred)

Ypred <- as.data.frame(fit, pars = "ymean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred)

Ypred <- Ypred%>%
  filter(timeseries >= 14)%>% filter(timeseries <= 155)

Cpred <- Cpred%>%
  filter(timeseries >= 14)%>% filter(timeseries <= 155)

dir.create(figDir)
dir.create(tabDir)
### posterior distributions of parameters
ptable <- monitor(as.array(fit, pars = parametersToPlot), warmup = 0, print = FALSE)
ptable[,c(1,4,8)]


## parameter estimate from ontogeny fit
## variation in the rate of Net loss with host age
r_neo_est = ptable[1, 1]

## writting the parameter values with CIs into a csv file
#write.csv(ptable, file = file.path(tabDir, paste(modelName, "_ParameterTable.csv", sep = "")))
#write.csv(ploocv, file = file.path(tabDir, paste(modelName, "_Informatiion_criteria.csv", sep = "")))

### ----------------------------------------------------------------------##
### plotting the rate Net loss
ts_new <- seq(0, 100, 0.5)  # prediction timeocurse

## number of cells in source compartment
source_influx <- function(t){
  
  theta0 = exp(13.9591976)
  r1 = 0.3485853
  n1 = 2.1701241
  
  theta = theta0 * (1 + t^n1 * exp(-r1 * t));
  return(theta)
}

## rate of infiltration into FM compaertment inferred from best model fit on adult mice data 
phi_simple <- function(t){
  
  psi =  0.7239
  tb = 40
  
  psi * source_influx(t)
}

### varying rate of net loss with host age in neonates
delta_complex <- function(t, r_neo){
  r_d = 0.00085065
  delta0 = 0.02761
  tb = 40
  
  ifelse(t <= tb,
         delta0 * exp(-r_neo * (t-tb)),
         delta0 * exp(- r_d * (t-tb))
  )
}


### rate of loss with ci bounds
r_pred <- as.data.frame(fit, pars = "r_neo") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) 


r_pred_bounds <- data.frame(lapply(1:nrow(r_pred), function(i){
  delta_complex(ts_new, r_pred$value[i])
}))

res_r_pred_bounds <- t(r_pred_bounds)
rownames(res_r_pred_bounds) <- NULL

sort_r_pred_bounds <- apply(res_r_pred_bounds,2,sort,decreasing=F)
delta_mean <- delta_complex(ts_new, r_neo_est)
delta_ub <- sort_r_pred_bounds[nrow(r_pred)*0.975,]
delta_lb <- sort_r_pred_bounds[nrow(r_pred)*0.025,]


##### Estimating cell-age distribution
G_a_delta <- function(a, t){
  
  phi_simple(t-a) * exp(- integrate(delta_complex, lower = t-a, upper = t, r_neo = r_neo_est)$value)
}

# vectorise
G_age_delta <- Vectorize(G_a_delta)


# normaalise for the pool size
Norm_age_dist_delta <- function(a, t){
  G_age_delta(a, t)/
    integrate(G_age_delta, lower = 0, upper = t, t=t)$value
}


#### Deriving the GFP intensity distribution of FM B cells in neonates
G_f_psi <- function(gfp, Time, parms){
  #func_value <- function(f, Time, parms){
  r_gfp = parms[1]
  gfpmax = exp(parms[2])
  
  ## exponential decay
  age = (log(gfpmax) - log(gfp))/r_gfp
  value = G_age_delta(age, Time)/(r_gfp * gfp)
  
  return(value)
}

### GFP overall intensity
GFP_intensity <- function(gfp, Time, parms){
  gfpmax = exp(parms[2])
  return(gfp * G_f_psi(gfp, Time, parms))
}

### Normalised GFP MFI
GFP_mfi <- function(Time, parms){
  gfpmax = exp(parms[2])
  r_gfp = exp(parms[1])
  
  # lowest gfp intensity possible at time 'Time'
  gfptime = gfpmax * exp(-r_gfp * Time)
  if(gfptime >= 0.1){
    gfpmin = gfptime
  } else {
    gfpmin = 0.1
  }
  
  answer = integrate(GFP_intensity, lower = gfpmin, upper = gfpmax, Time=Time, parms=parms)$value/
    integrate(G_f_psi, lower = gfpmin, upper = gfpmax, Time=Time, parms=parms)$value
  return(answer)
}

# Vectorisation
GFP_mfi_vec <- function(Time, parms){mapply(GFP_mfi, Time, MoreArgs = list(parms))}


# GFP positive cell density
GFP_pos <- function(Time, parms){
  gfpmax = exp(parms[2])
  #gfp_gate = exp(parms[3])
  integrate(GFP_intensity, lower = 1000, upper = gfpmax, Time=Time,  parms=parms)$value/
    integrate(G_f_psi, lower = 1000, upper = gfpmax, Time=Time, parms=parms)$value
}

# vectorisation
GFP_pos_vec <- function(Time, parms){mapply(GFP_pos, Time, MoreArgs = list(parms))}


### loading the data
gfp_data <- read.csv('GFP_MFI.csv') %>%
  rename(host_age = X)

gfp_gat <- gfp_data %>%
  gather(-host_age, key = FM_comp, value = MFI)

ggplot()+
  geom_point(data = gfp_gat, aes(x= host_age,  y= MFI, col = FM_comp)) +
  scale_color_manual(values = c("#0C71E0", 1), name = NULL, labels = c('FM GFP_pos', 'FM total')) +
  geom_line(aes(ts_gfp, Gfp_tot), size = 1.2) + 
  geom_line(aes(ts_gfp, Gfp_high), size = 1.2, col = "#0C71E0") + 
  labs(title=paste('GFP MFI assuming higher loss rate in neonates than adults'),  y= "MFI", x= "Host age (days)") +
  scale_x_continuous(trans = "log", limits = c(8, 128), breaks = c(8, 16, 32, 64, 128))+
  scale_y_continuous(trans = "log", limits = c(512, 16384), breaks = c(512, 1024, 2048, 4096, 8198, 16384)) +
  my_theme + theme(legend.position = c(0.85, 0.85))


# Log likelihood function for fitting GFP MFI to observed intensities
LL_fit <- function(params, boot_data) { 
  
  ## transforming data
  observed_data.df <- boot_data%>%
    mutate(log_tot = log(FM_tot),
           log_pos = log(FM_pos)) 
  
  ## log transformed counts from the model
  pred_tot <- log(GFP_mfi_vec(observed_data.df$host_age, params))
  #pred_pos <- log(GFP_pos_vec(observed_data.df$host_age, params))
  
  ## calculating the sun of squared residuals
  ssqres1 <- sum((pred_tot - observed_data.df$log_tot)^2)
  #ssqres2 <- sum((pred_pos - observed_data.df$log_pos)^2)
  
  k  <- length(params)                #number of unknown parameters 
  n1 <- nrow(observed_data.df)            #number of observations in dataset1
  
  #cost function
  #log-likelihood ignoring all the terms dependent only on the number of observations n
  #matrix multipltication of residual and transpose of residuals
  logl <-  - (n1/2)* log(ssqres1)  #- (n1/2)* log(ssqres2)  
  
  return(-logl)     #since optim minimizes the function by default, ML
} 

# initial guess
pars_vec <- c(-2.2, 8.42)
### optimising the LL function to maximise the loglikelihood
optim_fit <- optim(par = pars_vec, fn = LL_fit, boot_data = gfp_data,
                   method=c("Nelder-Mead"), hessian = F,  control = list(trace = 1, maxit = 1000))

params <- optim_fit$par
paste0("gfp half life = ", log(2)/exp(params[1]))
paste0("gfpmax = ",exp(params[2]))
AIC <- (2 * optim_fit$value) + 2 * length(params)
paste0("AIC = ", AIC)

# Time course for predictions
ts_gfp <- seq(14, 160, 0.5)

# predictions
Gfp_tot_pred <- GFP_mfi_vec(ts_gfp, params)
Gfp_high_pred <- GFP_pos_vec(ts_gfp, params)

#Plots
ggplot()+
  geom_point(data = gfp_gat, aes(x= host_age,  y= MFI, col = FM_comp), size =2) +
  scale_color_manual(values = c("#0C71E0", 1), name = NULL, labels = c('GFPpos', 'Total')) +
  geom_line(aes(ts_gfp, Gfp_tot_pred), size = 1.2) + 
  geom_line(aes(ts_gfp, Gfp_high_pred), size = 1.2, col = "#0C71E0") + 
  labs(title=paste('GFP MFI assuming higher loss rate in neonates than adults'),  y= "MFI within FM B cells", x= "Host age (days)") +
  scale_x_continuous(trans = "log", limits = c(12, 160), breaks = c(8, 16, 32, 64, 128))+
  scale_y_continuous(trans = "log", limits = c(512, 4096), breaks = c(512, 1024, 2048, 4096, 8198, 16384)) +
  my_theme + theme(legend.position = c(0.85, 0.91))



ggsave("gfp_lambda_total.pdf", plot = last_plot(), width = 6, height = 4.5, device = 'pdf',  useDingbats = FALSE)


### age sequence for predictions
age_seq <- seq(from = 0, to = 40, 0.01)

## open graphics device 
## saving  plots for quality control 
pdf(file = file.path(figDir, paste(modelName,"Plots%03d.pdf", sep = "")),
    width = 6, height = 4, onefile = F, useDingbats = FALSE)


ggplot() +
 geom_ribbon(data = Cpred, aes(x = timeseries, ymin = lb, ymax = ub), fill = "#3B9AB2", alpha = 0.1)+
 geom_ribbon(data = Ypred, aes(x = timeseries, ymin = lb, ymax = ub), fill = "#3B9AB2", alpha = 0.3)+
 geom_line(data = Ypred, aes(x = timeseries, y = median), col = "#3B9AB2", size=1) +
 geom_point(data = FM_counts, aes(x = age.at.S1K, y = total_counts, col= set_name)) +
 scale_color_manual(values = c("#F21A00", "#3B9AB2"), name = "Data") +
 labs(title=paste('Counts of FM B cells'),  y=NULL, x= "Host age (days)") + 
 scale_x_continuous(limits = c(10, 600), trans="log10", breaks=c(10, 30, 100, 300))+
 scale_y_continuous(limits = c(4e5, 1.2e8), trans="log10", breaks=c(1e4, 1e5, 1e6, 1e7, 1e8), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
 my_theme

ggplot() +
 geom_ribbon(aes(x = ts_new, ymin = delta_lb, ymax = delta_ub), fill = "#3B9AB2", alpha = 0.3)+
 geom_line(aes(x = ts_new, y = delta_mean), col = "#3B9AB2", size=1.5) +
 labs(title=paste('Rate of loss FM cells'),  y=NULL, x= "Host age (days)") + 
 scale_x_continuous(limits = c(0, 100),  breaks=c(0, 25, 50, 75, 100))+
 scale_y_continuous(limits = c(0, 0.7),  breaks=c(0, 0.2, 0.4, 0.6)) +
 my_theme

ggplot()+
 geom_line(aes(age_seq, Norm_age_dist_delta(age_seq, 40)), size = 1.25, col = "darkblue") +  ylim(0, 0.05) +
 labs(title=paste('Normalised Cell age distribution of FM B cells'),  y= "Frequency", x= "cell age (days)") +
   my_theme

#lay1 <- cbind(c(1,1,1,1),
#              c(1,1,1,1),
#              c(2,2,3,3))
#
#
#gridExtra::grid.arrange(p1, p2, p3, layout_matrix = lay1)


dev.off()



