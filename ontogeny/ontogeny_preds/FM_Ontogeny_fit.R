## log linear model -- for counts as spline

## clearing the environment
rm(list = ls())  
gc()    
setwd("~/Desktop/GIt_repos/FM_GC_ageCorrected")


library(rstan)
library(tidyverse)
library(parallel)
library(loo)

modelName <- "FM_Ontogeny_fit"
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

theta_spline <- function(t, theta0, r1, n1) {
  
  theta = theta0 * (1 + t^n1 * exp(-r1 * t))
  
  return( theta)
}

## import the data set
data_imp <- read_csv(dataDir)%>% 
  mutate(days_post_t0 = age.at.S1K - 14)%>% arrange(days_post_t0)

# time points in data
data_time <- data_imp$age.at.S1K
solve_time <- c(unique(data_time))   # unique time points to solve ode

#keep track of index of time point in relation to solve_time
time_index <- purrr::map_dbl(data_time, function(x) which(x == solve_time))

# time sequnce used for plotting 
ts_pred = seq(from = 0, to = 200, 0.5)

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
  lambda0 = 0.02761,
  Beta = 5.7562,
  r_d = 0.00085065,
  y0_Log = 14.08421,
  kappa0 = 0.1352286)

## create initial estimates
init <- function() list(
  r_psi_Log = rnorm(1, -6, 0.1),
  sigma = exp(rnorm(1,log(1.5), 1))
  )

ggplot()+
  geom_point(aes(ts_pred, theta_spline(ts_pred, data$theta0, data$r1, data$n1)))


## Specify the variables for which you want history and density plots
parametersToPlot <- c("r_psi", "sigma")

## Additional variables to monitor
otherRVs <- c("ymean_pred", "countspred", "log_lik")

parameters <- c(parametersToPlot, otherRVs)

################################################################################################
# run Stan

nChains <- 4
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
            #control = list(adapt_delta = 0.9),
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

# Total cell counts
Cpred <- as.data.frame(fit, pars = "countspred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred + 14)

Ypred <- as.data.frame(fit, pars = "ymean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred +14)

Ypred <- Ypred%>%
  filter(timeseries >= 1)%>% filter(timeseries <= 155)

Cpred <- Cpred%>%
  filter(timeseries >= 1)%>% filter(timeseries <= 155)

dir.create(figDir)
dir.create(tabDir)
### posterior distributions of parameters
ptable <- monitor(as.array(fit, pars = parametersToPlot), warmup = 0, print = FALSE)
ptable[,c(1,4,8)]

## writting the parameter values with CIs into a csv file
#write.csv(ptable, file = file.path(tabDir, paste(modelName, "_ParameterTable.csv", sep = "")))
#write.csv(ploocv, file = file.path(tabDir, paste(modelName, "_Informatiion_criteria.csv", sep = "")))

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
                       days_post_t0 = age.at.S1K - 14,
                       set_name ="Ontogeny") %>%
  select(contains("age"), contains("days"), contains("set"), contains("total"))

## FM counts from chimera data
FM_adults <- read.csv("data/counts_FM.csv") %>% 
  mutate(days_post_t0 = age.at.S1K - 14,
         set_name = "Chimera")%>%
  select(contains("set_name"), contains("t0"), contains("S1K"), contains("total_counts"))

### pooling ontogeny and chimera data together
FM_counts <- rbind(FM_young, FM_adults) %>% arrange(set_name)

r_psi_est = ptable[1,1]

### plotting the rate of source influx
psi_fun <- function(t, psi, r_psi){
  
  psi_var = psi * (1 - exp(-r_psi * t))
  
  return(psi_var)
}


psi_pred <- as.data.frame(fit, pars = "r_psi") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) 

ts_new <- seq(0, 100, 0.5)

psi_env <- matrix(nrow = nrow(psi_pred), ncol = length(ts_new))

for (i in 1:nrow(psi_pred)){
  
  psi_env[i, ] = psi_fun(ts_new, data$psi, psi_pred$value[i])
}


rownames(psi_env) <- NULL
psi_env_sort <- apply(psi_env, 2, sort, decreasing=F)

psi_mean <- psi_fun(ts_new, data$psi, r_psi_est)
psi_lb <- psi_env_sort[nrow(psi_pred) * 0.025, ]
psi_ub <- psi_env_sort[nrow(psi_pred) * 0.975, ]


#### Age distribution
source_influx <- function(t){
  
  theta0 = exp(13.9591976)
  r1 = 0.3485853
  n1 = 2.1701241
  
  theta = theta0 * (1 + t^n1 * exp(-r1 * t));
  return(theta)
}


phi <- function(t, r_psi){
  
  psi =  0.7239
  tb = 0
  
  psi_var = psi * (1 - exp(-r_psi * (t-tb)))
  
  psi_var * source_influx(t)
}

ggplot()+
  geom_point(aes(ts_pred, phi(ts_pred, r_psi_est)))

lambda_simple <- function(t){
  r_d = 0.00085065
  lambda0 = 0.02761
  tb = 40
  
  lambda0  * exp(- r_d * (t-tb))
}


G_a_psi <- function(a, t){
  
  phi(t-a, r_psi_est) * exp(- integrate(lambda_simple, lower = t-a, upper = t)$value)
}

G_age_psi <- Vectorize(G_a_psi)

Norm_age_dist_psi <- function(t){
  #G_age_psi(a, t)/
    integrate(G_age_psi, lower = 0, upper = t, t=t)$value
}

Norm_age_dist_psi(0)

age_seq3 <- seq(from = 0, to = 40, 0.01)



## open graphics device 
## saving  plots for quality control 
pdf(file = file.path(figDir, paste(modelName,"Plots%03d.pdf", sep = "")),
    width = 6, height = 4, onefile = F, useDingbats = FALSE)

#pairs(fit, pars = parametersToPlot)

p1 <- ggplot() +
  geom_ribbon(data = Cpred, aes(x = timeseries, ymin = lb, ymax = ub), fill = "#3B9AB2", alpha = 0.1)+
  geom_ribbon(data = Ypred, aes(x = timeseries, ymin = lb, ymax = ub), fill = "#3B9AB2", alpha = 0.3)+
  geom_line(data = Ypred, aes(x = timeseries, y = median), col = "#3B9AB2", size=1) +
  geom_point(data = FM_counts, aes(x = age.at.S1K, y = total_counts, col= set_name)) +
  scale_color_manual(values = c("#F21A00", "#3B9AB2"), name = "Data") +
  labs(title=paste('Counts of FM B cells'),  y=NULL, x= "Host age (days)") + 
  scale_x_continuous(limits = c(10, 600), trans="log10", breaks=c(10, 30, 100, 300))+
  scale_y_continuous(limits = c(4e5, 1.2e8), trans="log10", breaks=c(1e4, 1e5, 1e6, 1e7, 1e8), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  my_theme


p2 <- ggplot() +
  geom_ribbon(aes(x = ts_new, ymin = psi_lb, ymax = psi_ub), fill = "#3B9AB2", alpha = 0.3)+
  geom_line(aes(x = ts_new, y = psi_mean), col = "#3B9AB2", size=1) +
  labs(title=paste('Rate of transition of T1 cells into FM'),  y=NULL, x= "Host age (days)") + 
  scale_y_continuous(limits = c(0, 0.8), breaks=c(0, 0.2, 0.4, 0.6, 0.8)) +
  my_theme

p3 <- ggplot()+
  geom_line(aes(age_seq3, Norm_age_dist_psi(age_seq3, 40)), size = 1.25, col = "darkblue") + ylim(0, 0.05) +
  labs(title=paste('Normalised Cell age distribution of FM B cells'),  y= NULL, x= "cell age (days)") +
  my_theme

p1
p2
p3
lay1 <- rbind(c(1,1),
              c(2,3))


gridExtra::grid.arrange(p1, p2, p3, layout_matrix = lay1)

dev.off()



