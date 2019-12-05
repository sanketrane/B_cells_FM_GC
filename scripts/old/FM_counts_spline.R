rm(list = ls())
gc()

# loading the libraries required for smooth processing
require(tidyverse)

data_counts <- read.csv("data/counts_FM_T2.csv")

y_data <- data_counts$total_counts
x_data <- data_counts$age.at.S1K 

FM.fit <- lm((y_data) ~ x_data)
summary(FM.fit)

y0 <- FM.fit$coefficients[1]
nu <- FM.fit$coefficients[2]

counts_pred <- function(t, intc, m){
  intc + (m * t)
}

ts <- seq(0, 600)

ggplot()+
  geom_point(data = data_counts, aes(y=total_counts, x= age.at.S1K), col = "darkblue",  size=3)+
  geom_line(aes(y = (counts_pred(ts, (y0), nu)), x = ts), col = "darkblue", size =1.5)+
  geom_line(aes(y = (counts_pred(ts, (y0), 0)), x = ts), col ="darkred", linetype = 2, size =1.5)+
  labs(title=paste('Cell counts: FM'),  y=NULL, x="Host age (days)") + 
  scale_x_continuous(limits = c(0, 600), breaks = c(0, 100,200,300,400, 500))+
  scale_y_continuous(limits = c(5e6, 2e8), trans="log10", breaks=c(1e4, 1e5, 1e6, 1e7, 1e8)) 



#############################################################################################

## spline for T2
data_counts2 <- read.csv("data/T2_data.csv")

y_data <- data_counts2$total_counts
x_data <- data_counts2$age.at.S1K 

T2.fit <- lm(log(y_data) ~ x_data)
summary(T2.fit)

y0 <- T2.fit$coefficients[1]
nu <- T2.fit$coefficients[2]

counts_pred <- function(t, intc, m){
  intc + (m * t)
}

ts <- seq(0, 600)

ggplot()+
  geom_point(data = data_counts2, aes(y=total_counts, x= age.at.S1K), col = "darkblue",  size=3)+
  geom_line(aes(y = exp(counts_pred(ts, (y0), nu)), x = ts), col = "darkblue", size =1.5)+
  geom_line(aes(y = exp(counts_pred(ts, (y0), 0)), x = ts), col ="darkred", linetype = 2, size =1.5)+
  labs(title=paste('Cell counts: FM'),  y=NULL, x="Host age (days)") + 
  scale_x_continuous(limits = c(0, 600), breaks = c(0, 100,200,300,400, 500))+
  scale_y_continuous(limits = c(1e5, 1e7), trans="log10", breaks=c(1e4, 1e5, 1e6, 1e7, 1e8)) 


#############################################################################################

## spline for T1
data_counts1 <- read.csv("data/T1_data.csv")

y_data <- data_counts1$total_counts
x_data <- data_counts1$age.at.S1K 

T1.fit <- lm(log(y_data) ~ x_data)
summary(T1.fit)

y0 <- T1.fit$coefficients[1]
nu <- T1.fit$coefficients[2]

counts_pred <- function(t, intc, m){
  intc + (m * t)
}

ts <- seq(0, 600)

ggplot()+
  geom_point(data = data_counts1, aes(y=total_counts, x= age.at.S1K), col = "darkblue",  size=3)+
  geom_line(aes(y = exp(counts_pred(ts, (y0), nu)), x = ts), col = "darkblue", size =1.5)+
  geom_line(aes(y = exp(counts_pred(ts, (y0), 0)), x = ts), col ="darkred", linetype = 2, size =1.5)+
  labs(title=paste('Cell counts: FM'),  y=NULL, x="Host age (days)") + 
  scale_x_continuous(limits = c(0, 600), breaks = c(0, 100,200,300,400, 500))+
  scale_y_continuous(limits = c(1e5, 1e7), trans="log10", breaks=c(1e4, 1e5, 1e6, 1e7, 1e8)) 


