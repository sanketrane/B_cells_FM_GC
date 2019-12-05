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
counts_host <- readxl::read_excel(path = "data/Bcells_numbers.xlsx", sheet = 1)%>%
  #mutate_at(c("immature_FM"), funs(lead), n = 1 )%>%
  filter(is.na(Ki67)) %>% filter(is.na(notes))%>%
  #filter(row_number() %% 2 == 1)%>% filter(is.na(notes)) %>%
  select(contains("Lamis"), contains("days"), contains("age"), contains(target_pop)) %>% 
  mutate(Total_FM = SP.FM)%>%  select(-contains("SP"), -contains("LN"))%>% 
  mutate_if(is.character, as.numeric) %>% unique()

counts_donor <- readxl::read_excel(path = "data/Bcells_numbers.xlsx", sheet = 2)%>%
  #mutate_at(c("immature_FM"), funs(lead), n = 1 )%>%
  #filter(row_number() %% 2 == 1)%>% filter(is.na(notes)) %>%
  filter(is.na(Ki67)) %>% filter(is.na(notes))%>%
  select(contains("Lamis"), contains("days"), contains("age"), contains(target_pop))%>%
  mutate(Total_FM = SP.FM)%>%  select(-contains("SP"), -contains("LN"))%>% 
  mutate_if(is.character, as.numeric) %>% unique()

# ki67 for FM
ki67_host <- readxl::read_excel(path = "data/Bcells_numbers.xlsx", sheet = 1)%>%
  #filter(is.na(Ki67)) %>% filter(is.na(notes))
  filter(row_number() %% 2 == 0)%>% filter(is.na(notes)) %>%
  select(contains("Lamis"), contains("days"), contains("age"), contains(target_pop))%>% unique()%>%
  mutate(ki67Pos_FM = SP.FM)%>%  select(-contains("SP"), -contains("LN"))%>%
  mutate_if(is.character, as.numeric) %>% unique()

ki67_donor <- readxl::read_excel(path = "data/Bcells_numbers.xlsx", sheet = 2)%>%
  #filter(is.na(Ki67)) %>% filter(is.na(notes))
  filter(row_number() %% 2 == 0)%>% filter(is.na(notes)) %>%
  select(contains("Lamis"), contains("days"), contains("age"), contains(target_pop))%>% unique()%>%
  mutate(ki67Pos_FM = SP.FM)%>%  select(-contains("SP"), -contains("LN"))%>%
  mutate_if(is.character, as.numeric) %>% unique()

# merging total counts for host and donor compartments
FM_counts <- full_join(counts_host, counts_donor, by = c("Lamis.ID", "days.post.bmt", "age.at.S1K", "age.at.bmt"), suffix= c(".host", ".donor"))

# merging ki67Pos counts from host and donor
FM_ki67 <- full_join(ki67_host, ki67_donor, by = c("Lamis.ID", "days.post.bmt", "age.at.S1K", "age.at.bmt"), suffix= c(".host", ".donor"))

# calculating total counts, donor fractions anf ki67Pos fractions
FM_data <- full_join(FM_ki67, FM_counts, by = c("Lamis.ID", "days.post.bmt", "age.at.S1K", "age.at.bmt"), suffix= c(".host", ".donor"))%>%
  mutate(total_counts = Total_FM.host + Total_FM.donor,
         fd = Total_FM.donor/ total_counts,
         # ki67 fractions without the immature FM compartment
         host_ki67_FM = ki67Pos_FM.host/ Total_FM.host,
         donor_ki67_FM = ki67Pos_FM.donor/ Total_FM.donor)%>%
  select(-contains(".host"), -contains(".donor"))

# normalising donor fraction in FM by dividing with the donor fractions in the source compartment
FM_fd <- FM_data%>%
  select(-contains("ki67")) %>%
  full_join(source_data, by = c("Lamis.ID", "days.post.bmt", "age.at.S1K", "age.at.bmt"), suffix= c(".FM", ".source"))%>%
  mutate(Nfd = fd.FM/ fd.source) %>%
  select(contains("Lamis"), contains("days"), contains("age"), contains("Nfd"))


counts_FM <- FM_data%>%
  select(contains("Lamis"), contains("days"), contains("age"), contains("total_counts"))

write.csv(counts_FM, "~/Desktop/Git_repos/FM_cells_dynamics/data/counts_FM.csv")
write.csv(FM_fd, "~/Desktop/Git_repos/FM_cells_dynamics/data/Nfd_FM.csv")
write.csv(source_data, "~/Desktop/Git_repos/FM_cells_dynamics/data/source_data.csv")

ggplot(FM_data) +
  geom_point(aes(x = age.at.S1K, y = total_counts), col=4, size = 4)+
  scale_y_log10(limits = c(1e6, 1e8), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  scale_x_continuous(limits = c(50, 600), trans = "log10", breaks = c(300, 100, 30)) + 
  theme_bw() + labs(x = "Host age (days)", y = NULL, title = "Total counts: SP FM") +
  theme(axis.text = element_text(size = 22),
        axis.title =  element_text(size = 20, face = "bold"),
        plot.title = element_text(size=20,  hjust = 0.5),
        legend.text = element_text(size=20),
        legend.title = element_text(size = 20, face = "bold"))

ggplot(FM_fd) +
  geom_point(aes(x = days.post.bmt, y = Nfd, color = age.at.bmt),  size = 4)+
  geom_hline(yintercept = 1.00, linetype = 2, size =1.5, col="darkred")+
  ylim(0, 1.1) + theme_bw() + labs(x = "Days post t0", y = NULL, title = "Noirmalised donor fractions: SP FM") +
  theme(axis.text = element_text(size = 22),
        axis.title =  element_text(size = 20, face = "bold"),
        plot.title = element_text(size=20,  hjust = 0.5),
        legend.text = element_text(size=20),
        legend.title = element_text(size = 20, face = "bold"))


ki67_FM <-  FM_data %>%
  select(contains("Lamis"), contains("days"), contains("age"), contains("ki67"))%>%
  select(contains("Lamis"), contains("days"), contains("age"),  contains("FM"))%>%
  gather(-c(Lamis.ID, days.post.bmt, age.at.S1K, age.at.bmt), key = subpopulation, value = ki67_fraction)


ki67_source <-  source_data %>%
  select(contains("days"), contains("ki67"))%>%
  select(contains("days"), contains("Tra"))%>%
  gather(-days.post.bmt, key = subpopulation, value = ki67_fraction)


write.csv(ki67_FM, "~/Desktop/Git_repos/FM_cells_dynamics/data/ki67_FM.csv")

ggplot(ki67_FM)+
  geom_point(aes(x = days.post.bmt, y = ki67_fraction, color = subpopulation), alpha = 0.8, size = 4)+
  scale_color_manual(values = c(4,2), name = "Subset", labels = c("Donor", "Host"))+
  ylim(0, 1) + scale_x_log10(limits = c(10,500), breaks = c(30, 100, 300)) +  theme_bw()+
  labs(x = "Days post t0", y = NULL, title = "Proportions of ki67Hi cells: SP FM") +
  theme(axis.text = element_text(size = 22),
        axis.title =  element_text(size = 20, face = "bold"),
        plot.title = element_text(size=20,  hjust = 0.5),
        legend.text = element_text(size=20),
        legend.title = element_text(size = 20, face = "bold"))


ggplot(source_data) +
  geom_point(aes(x = age.at.S1K, y = total_counts), col=4, size = 4)+
  scale_y_log10(limits = c(1e5, 1e8), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  scale_x_continuous(limits = c(50, 600), trans = "log10", breaks = c(300, 100, 30)) + 
  theme_bw() + labs(x = "Host age (days)", y = NULL, title = "Total counts: Source") +
  theme(axis.text = element_text(size = 22),
        axis.title =  element_text(size = 20, face = "bold"),
        plot.title = element_text(size=20,  hjust = 0.5),
        legend.text = element_text(size=20),
        legend.title = element_text(size = 20, face = "bold"))

ggplot(source_data) +
  geom_point(aes(x = days.post.bmt, y = fd),  size = 4)+
  geom_hline(yintercept = 1.00, linetype = 2, size =1.5, col="darkred")+
  ylim(0, 1.2) + theme_bw() + labs(x = "Days post t0", y = NULL, title = "Noirmalised donor fractions: Source") +
  theme(axis.text = element_text(size = 22),
        axis.title =  element_text(size = 20, face = "bold"),
        plot.title = element_text(size=20,  hjust = 0.5),
        legend.text = element_text(size=20),
        legend.title = element_text(size = 20, face = "bold"))

ggplot(ki67_source)+
  geom_point(aes(x = days.post.bmt, y = ki67_fraction, color = subpopulation), alpha = 0.8, size = 4)+
  scale_color_manual(values = c(4,2), name = "Subset", labels = c("Donor", "Host"))+
  ylim(0, 1) + scale_x_log10(limits = c(10,500), breaks = c(30, 100, 300)) +  theme_bw()+
  labs(x = "Days post t0", y = NULL, title = "Proportions of ki67Hi cells: source") +
  theme(axis.text = element_text(size = 22),
        axis.title =  element_text(size = 20, face = "bold"),
        plot.title = element_text(size=20,  hjust = 0.5),
        legend.text = element_text(size=20),
        legend.title = element_text(size = 20, face = "bold"))


