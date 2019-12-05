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

### total counts

###### agebin #1

Cpred1 <- as.data.frame(fit, pars = "countspred_age1") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred1,
            "age_bins" = rep("Age bin#1: <8Wks", length(ts_pred1)))%>%
  filter(timeseries >= 60) %>%  filter(timeseries <= 460)

Y1pred1 <- as.data.frame(fit, pars = "y1_mean_pred_age1") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred1,
            "age_bins" = rep("Age bin#1: <8Wks", length(ts_pred1)))%>%
  filter(timeseries >= 60) %>%  filter(timeseries <= 460)

###### agebin #2
Cpred2 <- as.data.frame(fit, pars = "countspred_age2") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred2,
            "age_bins" = rep("Age bin#2: 8-12Wks", length(ts_pred2)))%>%
  filter(timeseries >= 72) %>%  filter(timeseries <= 450)

Y1pred2 <- as.data.frame(fit, pars = "y1_mean_pred_age2") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred2,
            "age_bins" = rep("Age bin#2: 8-12Wks", length(ts_pred2)))%>%
  filter(timeseries >= 72) %>%  filter(timeseries <= 450)

###### agebin #3
Cpred3 <- as.data.frame(fit, pars = "countspred_age3") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred3 ,
            "age_bins" = rep("Age bin#3: >12Wks", length(ts_pred3))) %>%
  filter(timeseries >= 110) %>%  filter(timeseries <= 575)

Y1pred3 <- as.data.frame(fit, pars = "y1_mean_pred_age3") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred3 ,
            "age_bins" = rep("Age bin#3: >12Wks", length(ts_pred3))) %>%
  filter(timeseries >= 110) %>%  filter(timeseries <= 575)

#gathering data frames
counts_median_pooled <- rbind(Y1pred1, Y1pred2, Y1pred3)
counts_sd_pooled <- rbind(Cpred1, Cpred2, Cpred3)


#### nfd plots
fdpred1 <- as.data.frame(fit, pars = "fdpred_age1") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred1,
            "age_bins" = rep("Age bin#1: <8Wks", length(ts_pred1))) %>%
  filter(timeseries >= 41) %>%  filter(timeseries <= 575)

Y2pred1 <- as.data.frame(fit, pars = "y2_mean_pred_age1") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred1,
            "age_bins" = rep("Age bin#1: <8Wks", length(ts_pred1))) %>%
  filter(timeseries >= 41) %>%  filter(timeseries <= 575)

fdpred2 <- as.data.frame(fit, pars = "fdpred_age2") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred2,
            "age_bins" = rep("Age bin#2: 8-12Wks", length(ts_pred2))) %>%
  filter(timeseries >= 61) %>%  filter(timeseries <= 575)

Y2pred2 <- as.data.frame(fit, pars = "y2_mean_pred_age2") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred2,
            "age_bins" = rep("Age bin#2: 8-12Wks", length(ts_pred2))) %>%
  filter(timeseries >= 61) %>%  filter(timeseries <= 575)


fdpred3 <- as.data.frame(fit, pars = "fdpred_age3") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred3 ,
            "age_bins" = rep("Age bin#3: >12Wks", length(ts_pred3))) %>%
  filter(timeseries >= 111) %>%  filter(timeseries <= 575)

Y2pred3 <- as.data.frame(fit, pars = "y2_mean_pred_age3") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred3,
            "age_bins" = rep("Age bin#3: >12Wks", length(ts_pred3))) %>%
  filter(timeseries >= 111) %>%  filter(timeseries <= 575)


#gathering data frames
nfd_median_pooled <- rbind(Y2pred1, Y2pred2, Y2pred3)
nfd_sd_pooled <- rbind(fdpred1, fdpred2, fdpred3)


####ki67 plots
donor_kiPred1 <- as.data.frame(fit, pars = "donor_kiprop_pred1") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred1  +12,
            "age_bins" = rep("Age bin#1: <8Wks", length(ts_pred1)))%>%
  filter(timeseries >= 72) %>%  filter(timeseries <= 475)

Y3pred1 <- as.data.frame(fit, pars = "y3_mean_pred1") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred1  +12,
            "age_bins" = rep("Age bin#1: <8Wks", length(ts_pred1)))%>%
  filter(timeseries >= 72) %>%  filter(timeseries <= 475)

host_kiPred1 <- as.data.frame(fit, pars = "host_kiprop_pred1") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred1  +12,
            "age_bins" = rep("Age bin#1: <8Wks", length(ts_pred1)))%>%
  filter(timeseries >= 72) %>%  filter(timeseries <= 475)

Y4pred1 <- as.data.frame(fit, pars = "y4_mean_pred1") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" =ts_pred1  +12,
            "age_bins" = rep("Age bin#1: <8Wks", length(ts_pred1)))%>%
  filter(timeseries >= 72) %>%  filter(timeseries <= 475)

donor_kiPred2 <- as.data.frame(fit, pars = "donor_kiprop_pred2") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred2,
            "age_bins" = rep("Age bin#2: 8-12Wks", length(ts_pred2)))%>%
  filter(timeseries >= 73) %>%  filter(timeseries <= 475)

Y3pred2 <- as.data.frame(fit, pars = "y3_mean_pred2") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred2,
            "age_bins" = rep("Age bin#2: 8-12Wks", length(ts_pred2)))%>%
  filter(timeseries >= 73) %>%  filter(timeseries <= 475)

host_kiPred2 <- as.data.frame(fit, pars = "host_kiprop_pred2") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred2,
            "age_bins" = rep("Age bin#2: 8-12Wks", length(ts_pred2)))%>%
  filter(timeseries >= 73) %>%  filter(timeseries <= 475)

Y4pred2 <- as.data.frame(fit, pars = "y4_mean_pred2") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred2,
            "age_bins" = rep("Age bin#2: 8-12Wks", length(ts_pred2)))%>%
  filter(timeseries >= 73) %>%  filter(timeseries <= 475)

donor_kiPred3 <- as.data.frame(fit, pars = "donor_kiprop_pred3") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" =  ts_pred3 - 7,
            "age_bins" = rep("Age bin#3: >12Wks", length(ts_pred3)))%>%
  filter(timeseries >= 114) %>%  filter(timeseries <= 575)

Y3pred3 <- as.data.frame(fit, pars = "y3_mean_pred3") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" =  ts_pred3 - 7,
            "age_bins" = rep("Age bin#3: >12Wks", length(ts_pred3)))%>%
  filter(timeseries >= 114) %>%  filter(timeseries <= 575)

host_kiPred3 <- as.data.frame(fit, pars = "host_kiprop_pred3") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" =  ts_pred3 - 7,
            "age_bins" = rep("Age bin#3: >12Wks", length(ts_pred3)))%>%
  filter(timeseries >= 114) %>%  filter(timeseries <= 575)

Y4pred3 <- as.data.frame(fit, pars = "y4_mean_pred3") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" =  ts_pred3 - 7,
            "age_bins" = rep("Age bin#3: >12Wks", length(ts_pred3)))%>%
  filter(timeseries >= 114) %>%  filter(timeseries <= 575)

#gathering data frames
donor_ki67_median_pooled <- rbind(Y3pred1, Y3pred2, Y3pred3) 
donor_ki67_median_pooled$key <- rep("Donor", nrow(donor_ki67_median_pooled))
host_ki67_median_pooled <- rbind(Y4pred1, Y4pred2, Y4pred3) 
host_ki67_median_pooled$key <- rep("host", nrow(host_ki67_median_pooled))

## variance
donor_ki67_sd_pooled <- rbind(donor_kiPred1, donor_kiPred2, donor_kiPred3)%>% select(-contains("key"))
host_ki67_sd_pooled <- rbind(host_kiPred1, host_kiPred2, host_kiPred3)%>% select(-contains("key"))

## combined donor host df
ki67_median_pooled <- rbind(donor_ki67_median_pooled, host_ki67_median_pooled)

#### facet plot for individual age bins
counts_facet <- ggplot() +
  geom_line(data = counts_median_pooled, aes(x = timeseries, y = median, col = age_bins)) +
  geom_ribbon(data = counts_median_pooled, aes(x = timeseries, ymin = lb, ymax = ub, fill = age_bins), alpha = 0.25)+
  geom_point(data = counts_binned, aes(x = age.at.S1K, y = total_counts, col = age_bins)) +
  scale_color_manual(values=c("#3B9AB2", "#E1AF00", "#F21A00"), name = NULL,
                     labels = c("<8Wks", "8-12Wks", ">12Wks"))+
  scale_fill_manual(values=c("#3B9AB2", "#E1AF00", "#F21A00"), name = NULL,
                    labels = c("<8Wks", "8-12Wks", ">12Wks"))+
  labs(title=paste('Total numbers of ', substr(modelName, 10, 13), " B cells"),  y=NULL, x="Host age (days)") + 
  scale_x_continuous(limits = c(60, 600), trans = "log10", breaks = c(75, 150, 300, 600))+
  scale_y_continuous(limits = c(5e6, 2e8), trans="log10", breaks=c(1e4, 1e5, 1e6, 1e7, 1e8), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  guides(color =FALSE, fill = FALSE) + myTheme + theme(legend.position = c(0.5, 0.85), legend.direction = "horizontal")

# facet plot for individual age bins
nfd_facet <- ggplot() +
  geom_line(data = nfd_median_pooled, aes(x = timeseries, y = median, col = age_bins)) +
  geom_ribbon(data = nfd_median_pooled, aes(x = timeseries, ymin = lb, ymax = ub, fill = age_bins), alpha = 0.25)+
  geom_point(data = Nfd_binned, aes(x = age.at.S1K, y = Nfd, col = age_bins)) +
  scale_color_manual(values=c("#3B9AB2", "#E1AF00", "#F21A00"), name = NULL,
                     labels = c("<8Wks", "8-12Wks", ">12Wks"))+
  scale_fill_manual(values=c("#3B9AB2", "#E1AF00", "#F21A00"), name = NULL,
                    labels = c("<8Wks", "8-12Wks", ">12Wks"))+
  guides(fill = FALSE) + myTheme + theme(legend.position = c(0.8, 0.3))

# combined plot
nfd_comb <- nfd_facet +
  geom_hline(aes(yintercept = 1), linetype = 2)+
  labs(x = "Host age (days)", y = NULL, title = paste("Normalised chimerism within ", substr(modelName, 10, 13), " B cells")) +
  scale_x_continuous(limits = c(40, 575), breaks = c(50, 200, 350, 500))+
  scale_y_continuous(limits =c(0, 1.1), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) 




####ki67 plots
## afcet plot for individual age bin
ki67_facet <- ggplot() +
  geom_line(data = donor_ki67_median_pooled, aes(x = timeseries, y = median), color = "#3B9AB2") +
  geom_ribbon(data = donor_ki67_median_pooled, aes(x = timeseries, ymin = lb, ymax = ub), fill = "#3B9AB2", alpha = 0.2)+
  geom_line(data = host_ki67_median_pooled, aes(x = timeseries, y = median), color = "#F21A00") +
  geom_ribbon(data = host_ki67_median_pooled, aes(x = timeseries, ymin = lb, ymax = ub), fill = "#F21A00", alpha = 0.2)+
  geom_point(data = ki67_binned, aes(x = age.at.S1K, y = prop_ki67hi, color = subpopulation)) +
  scale_alpha_manual(values = c(1, 0.7, 0.3), name = "Age bins", labels = c("<8Wks", "8-12Wks", "<12wks"))+
  scale_color_manual(values=c("#3B9AB2", "#F21A00"), name = NULL, labels = c("Donor", "Host"))+ 
  guides(alpha = FALSE) + myTheme + theme(legend.position = c(0.92, 0.9), legend.spacing = unit(0.2, "cm") ,
                                          legend.key.size = unit(0.5, "cm"), strip.background = element_blank())


# combined plot
ki67_comb <- ggplot() +
  geom_line(data = ki67_median_pooled, aes(x = timeseries, y = median, col = age_bins, linetype = key)) +
  #geom_ribbon(data = ki67_median_pooled, aes(x = timeseries, ymin = lb, ymax = ub, fill = age_bins), alpha = 0.25)+
  geom_point(data = ki67_binned, aes(x = age.at.S1K, y = prop_ki67hi, group = interaction(age_bins, subpopulation), color = age_bins, shape = subpopulation))+
  scale_color_manual(values=c("#3B9AB2", "#E1AF00", "#F21A00"), name = NULL,
                     labels = c("<8Wks", "8-12Wks", ">12Wks"))+
  scale_fill_manual(values=c("#3B9AB2", "#E1AF00", "#F21A00"), name = NULL,
                    labels = c("<8Wks", "8-12Wks", ">12Wks"))+
  scale_linetype_manual(values = c(1, 3), name = "",  labels = c("Donor", "Host"))+
  scale_shape_manual(values = c(19, 1), name = "",  labels = c("Donor", "Host"))+
  labs(x = "Host age (days)", y = NULL, title = paste("Proportions of host and donor ", substr(modelName, 10, 13), " B cells expressing Ki67")) +
  scale_x_continuous(limits = c(50, 475), trans = "log10",  breaks = c(75, 150, 300, 600))+
  scale_y_continuous(limits =c(0.0, 1.1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) +  myTheme +
  theme(legend.position = c(0.72, 0.8), legend.box = "horizontal", legend.spacing = unit(0.2, "cm") , legend.key.size = unit(0.5, "cm")) + guides(fill = FALSE) 


# plots for individual age bins
nfd_plot <- nfd_facet +
  geom_ribbon(data = nfd_sd_pooled, aes(x = timeseries, ymin = lb, ymax = ub, fill = age_bins), alpha = 0.1) +
  labs(x = "Host age (days)", y = NULL, title = paste("Normalised Donor fraction: ", substr(modelName, 10, 13))) +
  scale_x_continuous(limits = c(40, 575), breaks = c(50, 200, 350, 500))+
  scale_y_continuous(limits =c(0, 1.15), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) + guides(color = FALSE) +
  facet_wrap(~ age_bins)

counts_plot <- counts_facet +
  geom_ribbon(data = counts_sd_pooled, aes(x = timeseries, ymin = lb, ymax = ub, fill = age_bins), alpha = 0.1) + guides(color = FALSE) +
  facet_wrap(~ age_bins) 

ki67_plot <- ki67_facet +
  labs(x = "Host age (days)", y = NULL, title = paste("Proportions of host and donor ", substr(modelName, 10, 13), " B cells expressing Ki67")) +
  scale_x_continuous(limits = c(60, 600), trans = "log10",  breaks = c(75, 150, 300, 600))+
  scale_y_continuous(limits =c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) + 
  facet_wrap(~ age_bins)

# alternative style with error included
ki67_plot2 <- ki67_facet +
  geom_ribbon(data = donor_ki67_sd_pooled, aes(x = timeseries, ymin = lb, ymax = ub), fill = "#0000ba",  alpha = 0.15) +
  geom_ribbon(data = host_ki67_sd_pooled, aes(x = timeseries, ymin = lb, ymax = ub), fill = "#ed3e00", alpha = 0.15) +
  labs(x = "Host age (days)", y = NULL, title = paste("Proportions of host and donor ", substr(modelName, 10, 13), " B cells expressing Ki67")) +
  scale_x_continuous(limits = c(60, 600), trans = "log10",  breaks = c(75, 150, 300, 600))+
  scale_y_continuous(limits =c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) + guides(color = FALSE) +
  facet_wrap(~ age_bins)
