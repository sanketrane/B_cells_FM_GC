rm(list = ls()); gc()

library(tidyverse)

## posterior predictive distributions
options(bayesplot.base_size = 15,
        bayesplot.base_family = "sans")
color_scheme_set(scheme = "viridis")
my_theme <- theme(axis.text = element_text(size = 11),
                  axis.title =  element_text(size = 11, face = "bold"),
                  plot.title = element_text(size=11,  hjust = 0.5, face = "bold"),
                  legend.text = element_text(size=11), legend.background = element_blank(),
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

source_influx <- function(t){
  
  theta0 = exp(13.9591976)
  r1 = 0.3485853
  n1 = 2.1701241

  theta = theta0 * (1 + t^n1 * exp(-r1 * t));
  return(theta)
}


phi_complex <- function(t, r_psi){

  psi =  0.890718776
  tb = 40
  
  gaama = log((psi + r_psi)/r_psi) * (1/tb)
  
  psi_var = (psi + r_psi) * ( 1- exp(-r_psi * t))
  
  psi_var * source_influx(t)
}

ts_seq <- seq(from = 1, to = 400, 1)


ggplot()+
  geom_point(aes(ts_seq, source_influx(ts_seq))) + scale_y_log10()+ scale_x_log10()

ggplot()+
  geom_point(aes(ts_seq, phi(ts_seq, 0.03872))) + scale_y_log10() + scale_x_log10()

lambda_simple <- function(t){
  r_d = 0.00116935
  lambda0 = 0.02803697
  tb = 40
  
  lambda0  * exp(- r_d * (t-tb))
}

r_psi_est = 0.03358

ggplot()+
  geom_point(aes(ts_seq, lambda(ts_seq))) 


G_a_psi <- function(a, t){
  
  phi_complex(t-a, r_psi = r_psi_est) * exp(- integrate(lambda_simple, lower = 0, upper = a)$value)
}

G_age_psi <- Vectorize(G_a_psi)

Norm_age_dist_psi <- function(a, t){
  G_age_psi(a, t)/
    integrate(G_age_psi, lower = 0, upper = t, t=t)$value
}

age_seq1 <- seq(from = 0, to = 20, 0.01)
age_seq2 <- seq(from = 0, to = 30, 0.01)
age_seq3 <- seq(from = 0, to = 40, 0.01)

ggplot()+
  geom_point(aes(age_seq1, Norm_age_dist_psi(age_seq1, 20))) + 
  geom_point(aes(age_seq2, Norm_age_dist_psi(age_seq2, 30)), col = 2) +
  geom_point(aes(age_seq3, Norm_age_dist_psi(age_seq3, 40)), col = 3) 




lambda_neo_est = 0.776824

phi_simple <- function(t){
  
  psi =  0.890718776
  tb = 40
  
  psi * source_influx(t)
}

lambda_complex <- function(t, lambda_neo){
  r_d = 0.00116935
  lambda0 = 0.02803697
  tb = 40
  
  gaama = log(lambda_neo/lambda0) * (1/tb);
  
  ifelse(t <= tb,
         lambda_neo * exp(- gaama * t),
         lambda0 * exp(- r_d * (t-tb))
  )
}

ggplot()+
  geom_point(aes(x =  ts_seq, y =lambda_complex(ts_seq, lambda_neo_est )))

G_a_lambda <- function(a, t){
  
  phi_simple(t-a) * exp(- integrate(lambda_complex, lower = 0, upper = a, lambda_neo = lambda_neo_est)$value)
}

G_age_lambda <- Vectorize(G_a_lambda)

Norm_age_dist_lambda <- function(a, t){
  G_age_lambda(a, t)/
    integrate(G_age_lambda, lower = 0, upper = t, t=t)$value
}

age_seq1 <- seq(from = 0, to = 20, 0.01)
age_seq2 <- seq(from = 0, to = 100, 0.01)
age_seq3 <- seq(from = 0, to = 40, 0.01)

projectDir <- getwd()
figDir <- file.path(projectDir, "deliv", "figures")

ggplot()+
  geom_point(aes(age_seq3, Norm_age_dist_psi(age_seq3, 40)), size = 1, col = "darkblue") + 
  geom_point(aes(age_seq2, Norm_age_dist_psi(age_seq2, 100)), size = 1, col = "darkblue") + 
  labs(title=paste('Normalised Cell age distribution of FM B cells'),  y= "Frequency", x= "cell age (days)") +
  my_theme

ggplot()+
  #geom_point(aes(age_seq3, Norm_age_dist_psi(age_seq3, 40)), size = 1, col = "darkred") + 
  geom_line(aes(age_seq3, Norm_age_dist_lambda(age_seq3, 40)), size = 1.25, col = "darkred") +
  labs(title=paste('Normalised Cell age distribution of FM B cells'),  y= "Frequency", x= "cell age (days)") +
  my_theme






