### data wrangling for Marginal Zone compartmnet
rm(list = ls())
gc()

setwd("/opt/mesh/eigg/sanket/FM_GC_ageCorrected")

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

options(bayesplot.base_size = 15,
        bayesplot.base_family = "sans")
myTheme <- theme(text = element_text(size = 11), axis.text = element_text(size = 11), axis.title =  element_text(size = 11, face = "bold"),
                 plot.title = element_text(size=11,  hjust = 0.5, face = "bold"),
                 legend.background = element_blank(), legend.key = element_blank(),
                 legend.text = element_text(size= 9), legend.title = element_text(9))

# setting ggplot theme for rest fo the plots
theme_set(theme_bw())

library(tidyverse)

# source for FM 
parent_pop <- "Transitional"
target_pop <- "FM"

# total counts and donor fractions for the source poppulation
source_host <- readxl::read_excel(path = "data/Bcells_numbers.xlsx", sheet = 1) %>%
  filter(row_number() %% 2 == 1)%>% filter(is.na(notes)) %>%
  select(contains("Lamis"), contains("days"), contains("age"), contains(parent_pop))%>% 
  mutate(Total_T2 = SP.Transitional2 + LN.Transitional2 ,
         Total_T1 = SP.Transitional1,
         Total_Tra = Total_T1 + Total_T2)%>%  select(-contains("SP"), -contains("LN"))%>% 
  mutate_if(is.character, as.numeric) %>% unique()

source_donor <- readxl::read_excel(path = "data/Bcells_numbers.xlsx", sheet = 2) %>%
  filter(row_number() %% 2 == 1)%>% filter(is.na(notes)) %>%
  select(contains("Lamis"), contains("days"), contains("age"), contains(parent_pop))%>% 
  mutate(Total_T2 = SP.Transitional2 + LN.Transitional2 ,
         Total_T1 = SP.Transitional1,
         Total_Tra = Total_T1 + Total_T2)%>%  select(-contains("SP"), -contains("LN"))%>% 
  mutate_if(is.character, as.numeric) %>% unique()

## Ki67 data in the source compartmenmt
# total counts and donor fractions for the source poppulation
source_host_ki67 <- readxl::read_excel(path = "data/Bcells_numbers.xlsx", sheet = 1) %>%
  filter(row_number() %% 2 == 0)%>% filter(is.na(notes)) %>%
  select(contains("Lamis"), contains("days"), contains("age"), contains(parent_pop))%>%
  mutate(ki67Pos_T2 = SP.Transitional2 + LN.Transitional2 ,
         ki67Pos_T1 = SP.Transitional1,
         ki67Pos_Tra = ki67Pos_T1 + ki67Pos_T2) %>%  select(-contains("SP"), -contains("LN"))%>% 
  mutate_if(is.character, as.numeric) %>% unique()

source_donor_ki67 <- readxl::read_excel(path = "data/Bcells_numbers.xlsx", sheet = 2) %>%
  filter(row_number() %% 2 == 0)%>% filter(is.na(notes)) %>%
  select(contains("Lamis"), contains("days"), contains("age"), contains(parent_pop))%>%
  mutate(ki67Pos_T2 = SP.Transitional2 + LN.Transitional2 ,
         ki67Pos_T1 = SP.Transitional1,
         ki67Pos_Tra = ki67Pos_T1 + ki67Pos_T2) %>%  select(-contains("SP"), -contains("LN"))%>% 
  mutate_if(is.character, as.numeric) %>% unique()


# merging total counts for host and donor compartments
source_counts <- full_join(source_host, source_donor, by = c("Lamis.ID", "days.post.bmt", "age.at.S1K", "age.at.bmt"), suffix= c(".host", ".donor"))

# merging ki67Pos counts from host and donor
source_ki67 <- full_join(source_host_ki67, source_donor_ki67, by = c("Lamis.ID", "days.post.bmt", "age.at.S1K", "age.at.bmt"), suffix= c(".host", ".donor"))

# calculating total counts, donor fractions anf ki67Pos fractions
source_data <- full_join(source_ki67, source_counts, by = c("Lamis.ID", "days.post.bmt", "age.at.S1K", "age.at.bmt"), suffix= c(".host", ".donor"))%>%
  mutate(total_counts = Total_T1.host + Total_T1.donor,
         fd = (Total_T1.donor)/total_counts,
         # ki67 fractions without the immature FM compartment
         host_ki67_Tra = ki67Pos_T1.host/ Total_T1.host,
         donor_ki67_Tra = ki67Pos_T1.donor/ Total_T1.donor)%>%
  select(-contains(".host"), -contains(".donor"))


# total counts, donor fractions and ki67Pos fractions  for the target poppulation (FM)
counts_host_spl <- readxl::read_excel(path = "data/Bcells_numbers.xlsx", sheet = 1)%>%
  #mutate_at(c("immature_FM"), funs(lead), n = 1 )%>%
  filter(is.na(Ki67)) %>% filter(is.na(notes))%>%
  #filter(row_number() %% 2 == 1)%>% filter(is.na(notes)) %>%
  select(contains("Lamis"), contains("days"), contains("age"), contains(target_pop)) %>% 
  mutate(Total_FM = SP.FM)%>%  select(-contains("SP"), -contains("LN"))%>% 
  mutate_if(is.character, as.numeric) %>% unique()

counts_donor_spl <- readxl::read_excel(path = "data/Bcells_numbers.xlsx", sheet = 2)%>%
  #mutate_at(c("immature_FM"), funs(lead), n = 1 )%>%
  #filter(row_number() %% 2 == 1)%>% filter(is.na(notes)) %>%
  filter(is.na(Ki67)) %>% filter(is.na(notes))%>%
  select(contains("Lamis"), contains("days"), contains("age"), contains(target_pop))%>%
  mutate(Total_FM = SP.FM)%>%  select(-contains("SP"), -contains("LN"))%>% 
  mutate_if(is.character, as.numeric) %>% unique()

# ki67 for FM
ki67_host_spl <- readxl::read_excel(path = "data/Bcells_numbers.xlsx", sheet = 1)%>%
  #filter(is.na(Ki67)) %>% filter(is.na(notes))
  filter(row_number() %% 2 == 0)%>% filter(is.na(notes)) %>%
  select(contains("Lamis"), contains("days"), contains("age"), contains(target_pop))%>% unique()%>%
  mutate(ki67Pos_FM = SP.FM)%>%  select(-contains("SP"), -contains("LN"))%>%
  mutate_if(is.character, as.numeric) %>% unique()

ki67_donor_spl <- readxl::read_excel(path = "data/Bcells_numbers.xlsx", sheet = 2)%>%
  #filter(is.na(Ki67)) %>% filter(is.na(notes))
  filter(row_number() %% 2 == 0)%>% filter(is.na(notes)) %>%
  select(contains("Lamis"), contains("days"), contains("age"), contains(target_pop))%>% unique()%>%
  mutate(ki67Pos_FM = SP.FM)%>%  select(-contains("SP"), -contains("LN"))%>%
  mutate_if(is.character, as.numeric) %>% unique()

# total counts, donor fractions and ki67Pos fractions  for the target poppulation (FM)
counts_host_LN <- readxl::read_excel(path = "data/Bcells_numbers.xlsx", sheet = 1)%>%
  #mutate_at(c("immature_FM"), funs(lead), n = 1 )%>%
  filter(is.na(Ki67)) %>% filter(is.na(notes))%>%
  #filter(row_number() %% 2 == 1)%>% filter(is.na(notes)) %>%
  select(contains("Lamis"), contains("days"), contains("age"), contains(target_pop)) %>% 
  mutate(Total_FM = LN.FM)%>%  select(-contains("SP"), -contains("LN"))%>% 
  mutate_if(is.character, as.numeric) %>% unique()

counts_donor_LN <- readxl::read_excel(path = "data/Bcells_numbers.xlsx", sheet = 2)%>%
  #mutate_at(c("immature_FM"), funs(lead), n = 1 )%>%
  #filter(row_number() %% 2 == 1)%>% filter(is.na(notes)) %>%
  filter(is.na(Ki67)) %>% filter(is.na(notes))%>%
  select(contains("Lamis"), contains("days"), contains("age"), contains(target_pop))%>%
  mutate(Total_FM = LN.FM)%>%  select(-contains("SP"), -contains("LN"))%>% 
  mutate_if(is.character, as.numeric) %>% unique()

# ki67 for FM
ki67_host_LN <- readxl::read_excel(path = "data/Bcells_numbers.xlsx", sheet = 1)%>%
  #filter(is.na(Ki67)) %>% filter(is.na(notes))
  filter(row_number() %% 2 == 0)%>% filter(is.na(notes)) %>%
  select(contains("Lamis"), contains("days"), contains("age"), contains(target_pop))%>% unique()%>%
  mutate(ki67Pos_FM = LN.FM)%>%  select(-contains("SP"), -contains("LN"))%>%
  mutate_if(is.character, as.numeric) %>% unique()

ki67_donor_LN <- readxl::read_excel(path = "data/Bcells_numbers.xlsx", sheet = 2)%>%
  #filter(is.na(Ki67)) %>% filter(is.na(notes))
  filter(row_number() %% 2 == 0)%>% filter(is.na(notes)) %>%
  select(contains("Lamis"), contains("days"), contains("age"), contains(target_pop))%>% unique()%>%
  mutate(ki67Pos_FM = LN.FM)%>%  select(-contains("SP"), -contains("LN"))%>%
  mutate_if(is.character, as.numeric) %>% unique()

# merging total counts for host and donor compartments
FM_counts_spl <- full_join(counts_host_spl, counts_donor_spl, by = c("Lamis.ID", "days.post.bmt", "age.at.S1K", "age.at.bmt"), suffix= c(".host", ".donor"))

# merging ki67Pos counts from host and donor
FM_ki67_spl <- full_join(ki67_host_spl, ki67_donor_spl, by = c("Lamis.ID", "days.post.bmt", "age.at.S1K", "age.at.bmt"), suffix= c(".host", ".donor"))

# calculating total counts, donor fractions anf ki67Pos fractions
FM_data_spl <- full_join(FM_ki67_spl, FM_counts_spl, by = c("Lamis.ID", "days.post.bmt", "age.at.S1K", "age.at.bmt"), suffix= c(".host", ".donor"))%>%
  mutate(total_counts = Total_FM.host + Total_FM.donor,
         fd = Total_FM.donor/ total_counts,
         # ki67 fractions without the immature FM compartment
         host_ki67_FM = ki67Pos_FM.host/ Total_FM.host,
         donor_ki67_FM = ki67Pos_FM.donor/ Total_FM.donor)%>%
  select(-contains(".host"), -contains(".donor"))


# merging total counts for host and donor compartments
FM_counts_LN <- full_join(counts_host_LN, counts_donor_LN, by = c("Lamis.ID", "days.post.bmt", "age.at.S1K", "age.at.bmt"), suffix= c(".host", ".donor"))

# merging ki67Pos counts from host and donor
FM_ki67_LN <- full_join(ki67_host_LN, ki67_donor_LN, by = c("Lamis.ID", "days.post.bmt", "age.at.S1K", "age.at.bmt"), suffix= c(".host", ".donor"))

# calculating total counts, donor fractions anf ki67Pos fractions
FM_data_LN <- full_join(FM_ki67_LN, FM_counts_LN, by = c("Lamis.ID", "days.post.bmt", "age.at.S1K", "age.at.bmt"), suffix= c(".host", ".donor"))%>%
  mutate(total_counts = Total_FM.host + Total_FM.donor,
         fd = Total_FM.donor/ total_counts,
         # ki67 fractions without the immature FM compartment
         host_ki67_FM = ki67Pos_FM.host/ Total_FM.host,
         donor_ki67_FM = ki67Pos_FM.donor/ Total_FM.donor)%>%
  select(-contains(".host"), -contains(".donor"))

# normalising donor fraction in FM by dividing with the donor fractions in the source compartment
fd_spl <- FM_data_spl %>%
  select(-contains("ki67")) %>%
  full_join(source_data, by = c("Lamis.ID", "days.post.bmt", "age.at.S1K", "age.at.bmt"), suffix= c(".FM", ".source"))%>%
  mutate(Nfd = fd.FM/ fd.source) %>%
  select(contains("Lamis"), contains("days"), contains("age"), contains("Nfd"))

fd_LN <- FM_data_LN %>%
  select(-contains("ki67")) %>%
  full_join(source_data, by = c("Lamis.ID", "days.post.bmt", "age.at.S1K", "age.at.bmt"), suffix= c(".FM", ".source"))%>%
  mutate(Nfd = fd.FM/ fd.source) %>%
  select(contains("Lamis"), contains("days"), contains("age"), contains("Nfd"))

fd_pooled <- full_join(fd_spl, fd_LN, by = c("Lamis.ID", "days.post.bmt", "age.at.S1K", "age.at.bmt"), suffix = c(".SP", ".LN"))%>%
  gather(-c(Lamis.ID, days.post.bmt, age.at.S1K, age.at.bmt), key = subpopulation, value = Nfd)


# plolt
fd_plot <- ggplot(fd_pooled) +
  geom_point(data = fd_pooled, aes(x = age.at.S1K, y = Nfd, colour = subpopulation)) +
  geom_hline(yintercept = 1.00, linetype = 2) +
  labs(x = "Host age (days)", y = NULL, title = "Normalised chimerism within FM B cells") + 
  scale_x_continuous(limits = c(40, 575), breaks = c(50, 200, 350, 500))+
  scale_y_continuous(limits =c(0, 1.15), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) + 
  scale_color_manual(values = wesanderson::wes_palette(n=2, name = "Darjeeling1"), name = NULL, labels = c("Lymph nodes", "Spleen"))+
  myTheme + theme(legend.position = c(0.8, 0.2))


# total nunbers
counts_spl <- FM_data_spl %>%
  select(contains("Lamis"), contains("days"), contains("age"), contains("counts")) 

counts_LN <- FM_data_LN %>%
  select(contains("Lamis"), contains("days"), contains("age"), contains("counts"))

counts_pooled <- full_join(counts_spl, counts_LN, by = c("Lamis.ID", "days.post.bmt", "age.at.S1K", "age.at.bmt"), suffix = c(".SP", ".LN"))%>%
  gather(-c(Lamis.ID, days.post.bmt, age.at.S1K, age.at.bmt), key = subpopulation, value = total_counts)

# plot
counts_plot <- ggplot(counts_pooled) +
  geom_point(aes(x = age.at.S1K, y = total_counts, colour= subpopulation)) +
  scale_y_log10(limits = c(5e5, 1e8), breaks = c(1e6, 1e7, 1e8),  minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  scale_x_continuous(limits = c(60, 600), trans = "log10", breaks = c(75, 150, 300, 500)) + 
  scale_color_manual(values = wesanderson::wes_palette(n=2, name = "Darjeeling1"), name = NULL, labels = c("Lymph nodes", "Spleen"))+
  labs(x = "Host age (days)", y = NULL, title = "Total numbers of FM B cells") + myTheme + theme(legend.position = c(0.2, 0.2))


# ki67
ki67_spl <-  FM_data_spl %>%
  select(contains("Lamis"), contains("days"), contains("age"), contains("ki67"))%>%
  select(contains("Lamis"), contains("days"), contains("age"),  contains("FM"))%>%
  gather(-c(Lamis.ID, days.post.bmt, age.at.S1K, age.at.bmt), key = subpopulation, value = ki67_fraction)

ki67_LN <-  FM_data_LN %>%
  select(contains("Lamis"), contains("days"), contains("age"), contains("ki67"))%>%
  select(contains("Lamis"), contains("days"), contains("age"),  contains("FM"))%>%
  gather(-c(Lamis.ID, days.post.bmt, age.at.S1K, age.at.bmt), key = subpopulation, value = ki67_fraction)


ki67_pooled <- full_join(ki67_spl, ki67_LN, by = c("Lamis.ID", "days.post.bmt", "age.at.S1K", "age.at.bmt", "subpopulation"), suffix = c(".SP", ".LN"))%>%
  gather(-c(Lamis.ID, days.post.bmt, age.at.S1K, age.at.bmt, subpopulation), key = subfacet, value = ki67_fraction)

#facet_names <- c("ki67_fraction.SP" = "Spleen", "ki67_fraction.LN" = "Lymph nodes")

# plot
ki67_plot <- ggplot(ki67_pooled)+
  geom_point(aes(x = age.at.S1K, y = ki67_fraction, color = subfacet))+
  scale_color_manual(values = wesanderson::wes_palette(n=2, name = "Darjeeling1"), name = NULL, labels = c("Lymph nodes", "Spleen"))+
  #scale_color_manual(values=c("#00A08A", "#F98400"), name = NULL,
   #                  labels = c("Lymph nodes", "Spleen"))+
  ylim(0, 1) + scale_x_log10(limits = c(60, 600), breaks = c(75, 150, 300, 500)) +  theme_bw()+
  labs(x = "Host age (days)", y = NULL, title = "Proportions of FM B cells expressing Ki67") +
  myTheme + theme(legend.position = c(0.9, 0.85), legend.background = element_blank(), strip.background = element_blank(), strip.text = element_blank()) +
  facet_wrap(~ subpopulation) #, labeller = as_labeller(facet_names)) 


pdf(file = file.path("~/Dropbox/LNGC_plots", paste("FM_plots.pdf", sep = "")),
    width = 7, height = 5, onefile = FALSE, useDingbats = FALSE)

toprow <- cowplot::plot_grid(counts_plot, fd_plot,  labels = c("A", "B"), ncol = 2)
cowplot::plot_grid(toprow, ki67_plot,  labels = c("", "C"), nrow = 2)

dev.off()


