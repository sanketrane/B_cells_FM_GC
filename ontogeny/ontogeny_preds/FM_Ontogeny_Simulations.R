rm(list = ls()); gc()

library(rstan)
library(tidyverse)
library(deSolve)

#### working directories
setwd("~/Desktop/GIt_repos/FM_GC_ageCorrected")

### for better looking plots
## function to make plot labels look scientific
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

my_theme <- theme(axis.text = element_text(size = 12),
                  axis.title =  element_text(size = 12, face = "bold"),
                  plot.title = element_text(size= 12,  hjust = 0.5, face = "bold"),
                  legend.text = element_text(size= 12), legend.background = element_blank(),
                  legend.title = element_text(size = 11, face = "bold"), legend.position = c(0.85, 0.2))


################

#### importing functions written in Stan --  
### ODEs that were used to fit the time-dependent turnover model 
expose_stan_functions("models/FM_Ontogeny_sim.stan")

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

## parameter estimates from model fits to chimera data (time-dependent turnover model -- T1 source)
parms_T1 <- c("psi" = 0.7239, "rho" = 0.00599, "delta0" = 0.027615,
           "Beta" = 5.75624, "r_d" = 0.00085065)

## parameter estimates from model fits to chimera data (time-dependent turnover model -- T1 source)
parms_T2 <- c("psi" = 0.760715662, "rho" = 0.001173273, "lambda0" = 0.03036215,
           "Beta" = 4.81317874, "r_d" = 0.000853235)


### input values for splines of source population 
### two different spline functions are used viz exp and piecewise (pw)
### these values are estimated by fitting respective spline functions to counts of T1 or T2 cells in neonatal mice from ontogeny data
rdata_T1_exp <- c("theta0" = 13.9591976, "r1" = 0.3485853, "n1" = 2.1701241, "eps" = 0.97)
rdata_T2_exp <- c("theta0" = 14.1552612, "r1" = 0.2630012, "n1" = 1.8844567, "eps" = 0.97)
#rdata_T1_pw <- c("theta0" = 14.645994216, "r1" = 0.161825360, "r3" = 0.001055222, "alpha" = 0.837282136 , "eps" = 0.97)
#rdata_T2_pw <- c("theta0" = 15.78696337, "r1" = 0.02083275, "r3" = 0.07395444, "alpha" = 0.30970166 , "eps" = 0.97)

### time sequence for predictions
ts_pred <- seq(0, 550, 0.5)
numPred <- length(ts_pred)   #### dimensions of array for the stan function

##### stand model simulation
FM_sim <- function(source_spline, parms, rdata){
  
  ## initial conditions
  y0_init = mean(FM_young$total_counts[1:3])    #### the initial count of FM at t0 used for simulation
  ### kappa0 is the proportion of ki67hi cells at t0
  ### assumed to be equal to the mean proportion of ki67hi cells in adults
  kappa0 = 0.1352286
  
  init_cond <- c(y0_init * kappa0, y0_init * (1 - kappa0))  
  
  #### Depdending on the spline function used the ODE function to solve varies 
  ### simulation were made using different source popualtions (T1 or T2 cells) and different spline functions were tried
  ### possible values for the source_spline argument are -- "T1_exp", "T2_exp", "T1_pw", "T2_pw"
  ### identifying the source population and the spline function used.
  ### for more details of spline functions check the stan file.
  
  if (grepl("exp", source_spline) == TRUE) {
    ### for exp spline
    TDT_FM_df <- solve_ode_exp(ts_pred, init_cond, parms, rdata, numPred)
    spline_fun <- "exp"
  } else if(grepl("pw", source_spline) == TRUE){
    ### for pw spline
    TDT_FM_df <- solve_ode_pw(ts_pred, init_cond, parms, rdata, numPred)
    spline_fun <- "pw"
  }
  
  stan_pred_df <- data.frame("host_age" = ts_pred + 14,
                             "y_pred" = matrix(unlist(TDT_FM_df), nrow = length(TDT_FM_df), byrow = TRUE))%>%
    mutate(total_counts = y_pred.1 + y_pred.2) %>%
    select(contains("age"), contains("total"))
  
  print(spline_fun)
  return(stan_pred_df)
}

### simulations using T1 or T2 as source and exp spline
FM_sim_exp <- data.frame("T1_source" = FM_sim("T1_exp", parms_T1, rdata_T1_exp))%>%#,
                        # "T2_source" = FM_sim("T2_exp", parms_T2, rdata_T2_exp)) %>%
  mutate(host_age = T1_source.host_age)%>% select(-contains("source.host_age")) %>%
  gather(-host_age, key = source_pop, value = total_counts)

### simulations using T1 or T2 as source and pw spline
#FM_sim_pw <- data.frame("T1_source" = FM_sim("T1_pw", parms_T1, rdata_T1_pw),
#                         "T2_source" = FM_sim("T2_pw", parms_T2, rdata_T2_pw)) %>%
#  mutate(host_age = T1_source.host_age)%>% select(-contains("source.host_age")) %>%
#  gather(-host_age, key = source_pop, value = total_counts)


### plots
pdf(file = file.path("deliv/figures/FM_ontogeny_sim", paste("FM_", "exp", ".pdf", sep = "")),
    width = 6, height = 4, onefile = F, useDingbats = FALSE)

ggplot() +
  geom_line(data= FM_sim_exp, aes(x = host_age, y = total_counts),  col= "#3B9AB2", size = 1.2) +
  geom_point(data= FM_counts, aes(x = age.at.S1K, y = total_counts, col = set_name))+
  scale_y_log10(limits = c(5e5,1e8), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  scale_x_log10(limits = c(10, 600)) + scale_color_manual(values = c("#F21A00", "#3B9AB2"), name = "Data") +
  #scale_linetype_discrete(labels = c("T1", "T2"), name = "Precursor") +
  labs(x = "Host age (days)",  y = NULL, title = "Simulated counts of FM cells") + theme_bw() + my_theme

dev.off()











