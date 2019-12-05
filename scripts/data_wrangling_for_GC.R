### data wrangling for Marginal Zone compartmnet
rm(list = ls())
gc()

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
parent_pop <- "Transitional1"
target_pop <- "LN.GCs"

# total counts and donor fractions for the source poppulation
source_host <- readxl::read_excel(path = "data/Bcells_numbers.xlsx", sheet = 1) %>%
  filter(row_number() %% 2 == 1)%>% filter(is.na(notes)) %>%
  select(contains("Lamis"), contains("days"), contains("age"), contains(parent_pop))%>% 
  mutate(Total_counts = SP.Transitional1)%>%  select(-contains("SP"), -contains("LN"))%>% 
  mutate_if(is.character, as.numeric) %>% unique()

source_donor <- readxl::read_excel(path = "data/Bcells_numbers.xlsx", sheet = 2) %>%
  filter(row_number() %% 2 == 1)%>% filter(is.na(notes)) %>%
  select(contains("Lamis"), contains("days"), contains("age"), contains(parent_pop))%>% 
  mutate(Total_counts = SP.Transitional1 + LN.Transitional1)%>%  select(-contains("SP"), -contains("LN"))%>% 
  mutate_if(is.character, as.numeric) %>% unique()

## Ki67 data in the source compartmenmt
# total counts and donor fractions for the source poppulation
source_host_ki67 <- readxl::read_excel(path = "data/Bcells_numbers.xlsx", sheet = 1) %>%
  filter(row_number() %% 2 == 0)%>% filter(is.na(notes)) %>%
  select(contains("Lamis"), contains("days"), contains("age"), contains(parent_pop))%>%
  mutate(ki67Pos_source = SP.Transitional1) %>%  select(-contains("SP"), -contains("LN"))%>% 
  mutate_if(is.character, as.numeric) %>% unique()

source_donor_ki67 <- readxl::read_excel(path = "data/Bcells_numbers.xlsx", sheet = 2) %>%
  filter(row_number() %% 2 == 0)%>% filter(is.na(notes)) %>%
  select(contains("Lamis"), contains("days"), contains("age"), contains(parent_pop))%>%
  mutate(ki67Pos_source = SP.Transitional1) %>%  select(-contains("SP"), -contains("LN"))%>% 
  mutate_if(is.character, as.numeric) %>% unique()

# merging total counts for host and donor compartments
source_counts <- full_join(source_host, source_donor, by = c("Lamis.ID", "days.post.bmt", "age.at.S1K", "age.at.bmt"), suffix= c(".host", ".donor"))

# merging ki67Pos counts from host and donor
source_ki67 <- full_join(source_host_ki67, source_donor_ki67, by = c("Lamis.ID", "days.post.bmt", "age.at.S1K", "age.at.bmt"), suffix= c(".host", ".donor"))

# calculating total counts, donor fractions anf ki67Pos fractions
source_data <- full_join(source_ki67, source_counts, by = c("Lamis.ID", "days.post.bmt", "age.at.S1K", "age.at.bmt"), suffix= c(".host", ".donor"))%>%
  mutate(total_counts = Total_counts.host + Total_counts.donor,
         fd = (Total_counts.donor)/total_counts,
         # ki67 fractions without the immature FM compartment
         host_ki67_source = ki67Pos_source.host/ Total_counts.host,
         donor_ki67_source = ki67Pos_source.donor/ Total_counts.donor)%>%
  select(-contains(".host"), -contains(".donor"))


# total counts, donor fractions and ki67Pos fractions  for the target poppulation (FM)
counts_host <- readxl::read_excel(path = "data/Bcells_numbers.xlsx", sheet = 1)%>%
  #mutate_at(c("immature_FM"), funs(lead), n = 1 )%>%
  filter(is.na(Ki67)) %>% filter(is.na(notes))%>%
  #filter(row_number() %% 2 == 1)%>% filter(is.na(notes)) %>%
  select(contains("Lamis"), contains("days"), contains("age"), contains(target_pop)) %>% 
  mutate(Total_counts = LN.GCs)%>%  select(-contains("SP"), -contains("LN"))%>% 
  mutate_if(is.character, as.numeric) %>% unique()%>% na.omit()

counts_donor <- readxl::read_excel(path = "data/Bcells_numbers.xlsx", sheet = 2)%>%
  #mutate_at(c("immature_FM"), funs(lead), n = 1 )%>%
  #filter(row_number() %% 2 == 1)%>% filter(is.na(notes)) %>%
  filter(is.na(Ki67)) %>% filter(is.na(notes))%>%
  select(contains("Lamis"), contains("days"), contains("age"), contains(target_pop))%>%
  mutate(Total_counts = LN.GCs)%>%  select(-contains("SP"), -contains("LN"))%>% 
  mutate_if(is.character, as.numeric) %>% unique() %>% na.omit()


# ki67 for FM
ki67_host <- readxl::read_excel(path = "data/Bcells_numbers.xlsx", sheet = 1)%>%
  #filter(is.na(Ki67)) %>% filter(is.na(notes))
  filter(row_number() %% 2 == 0)%>% filter(is.na(notes)) %>%
  select(contains("Lamis"), contains("days"), contains("age"), contains(target_pop))%>% unique()%>%
  mutate(ki67Pos_GC = LN.GCs)%>%  select(-contains("SP"), -contains("LN"))%>%
  mutate_if(is.character, as.numeric) %>% unique() %>% na.omit()

ki67_donor <- readxl::read_excel(path = "data/Bcells_numbers.xlsx", sheet = 2)%>%
  #filter(is.na(Ki67)) %>% filter(is.na(notes))
  filter(row_number() %% 2 == 0)%>% filter(is.na(notes)) %>%
  select(contains("Lamis"), contains("days"), contains("age"), contains(target_pop))%>% unique()%>%
  mutate(ki67Pos_GC = LN.GCs)%>%  select(-contains("SP"), -contains("LN"))%>%
  mutate_if(is.character, as.numeric) %>% unique() %>% na.omit()

# merging total counts for host and donor compartments
GC_counts <- full_join(counts_host, counts_donor, by = c("Lamis.ID", "days.post.bmt", "age.at.S1K", "age.at.bmt"), suffix= c(".host", ".donor"))

# merging ki67Pos counts from host and donor
GC_ki67 <- full_join(ki67_host, ki67_donor, by = c("Lamis.ID", "days.post.bmt", "age.at.S1K", "age.at.bmt"), suffix= c(".host", ".donor"))

# calculating total counts, donor fractions anf ki67Pos fractions
GC_data <- full_join(GC_ki67, GC_counts, by = c("Lamis.ID", "days.post.bmt", "age.at.S1K", "age.at.bmt"), suffix= c(".host", ".donor"))%>%
  mutate(total_counts = Total_counts.host + Total_counts.donor,
         fd = Total_counts.donor/ total_counts,
         # ki67 fractions without the immature FM compartment
         host_ki67_GC = ki67Pos_GC.host/ Total_counts.host,
         donor_ki67_GC = ki67Pos_GC.donor/ Total_counts.donor)%>%
  select(-contains(".host"), -contains(".donor"))

countless_GC <- GC_data%>%
  select(contains("Lamis"), contains("days"), contains("age"), contains("total_counts"))%>%
  select(-contains("total_counts"))

counts_source <- left_join(countless_GC, source_data, by = c("Lamis.ID", "days.post.bmt", "age.at.S1K", "age.at.bmt"), suffix= c(".host", ".donor"))


# normalising donor fraction in FM by dividing with the donor fractions in the source compartment
GC_fd <- GC_data%>%
  select(-contains("ki67")) %>%
  full_join(counts_source, by = c("Lamis.ID", "days.post.bmt", "age.at.S1K", "age.at.bmt"), suffix= c(".GC", ".source"))%>%
  mutate(Nfd = fd.GC/ fd.source) %>%
  filter(Nfd < 1.5)%>%
  select(contains("Lamis"), contains("days"), contains("age"), contains("Nfd"))%>% na.omit()

counts_GC <- left_join(GC_fd, GC_data, by = c("Lamis.ID", "days.post.bmt", "age.at.S1K", "age.at.bmt"))%>%
  select(contains("Lamis"), contains("days"), contains("age"), contains("total_counts"))


ggplot(GC_data) +
  geom_point(aes(x = age.at.S1K, y = total_counts), col=4, size = 4)+
  scale_y_log10(limits = c(1e3, 5e7), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  scale_x_continuous(limits = c(50, 600), trans = "log10", breaks = c(300, 100, 30)) + 
  theme_bw() + labs(x = "Host age (days)", y = NULL, title = "Total counts: SP FM") +
  theme(axis.text = element_text(size = 22),
        axis.title =  element_text(size = 20, face = "bold"),
        plot.title = element_text(size=20,  hjust = 0.5),
        legend.text = element_text(size=20),
        legend.title = element_text(size = 20, face = "bold"))

ggplot(GC_fd) +
  geom_point(aes(x = age.at.S1K, y = Nfd),  size = 4)+
  geom_hline(yintercept = 1.00, linetype = 2, size =1.5, col="darkred")+ xlim(0, 500)+
  ylim(0, 1.2) + theme_bw() + labs(x = "Days post t0", y = NULL, title = "Noirmalised donor fractions: SP FM") +
  theme(axis.text = element_text(size = 22),
        axis.title =  element_text(size = 20, face = "bold"),
        plot.title = element_text(size=20,  hjust = 0.5),
        legend.text = element_text(size=20),
        legend.title = element_text(size = 20, face = "bold"))


ki67_data <-  GC_data %>%
  right_join(GC_fd, by = c("Lamis.ID", "days.post.bmt", "age.at.S1K", "age.at.bmt"))%>%
  select(contains("days"), contains("age"), contains("ki67"))%>%
  select(contains("days"), contains("age"),  contains("GC"))

ki67_GC <- ki67_data%>%
  gather(-c(days.post.bmt, age.at.S1K, age.at.bmt), key = subpopulation, value = ki67_fraction)


ggplot(ki67_GC)+
  geom_point(aes(x = days.post.bmt, y = ki67_fraction, color = subpopulation), alpha = 0.8, size = 4)+
  scale_color_manual(values = c(4,2), name = "Subset", labels = c("Donor", "Host"))+
  ylim(0, 1) + scale_x_log10(limits = c(10,500), breaks = c(30, 100, 300)) +  theme_bw()+
  labs(x = "Days post t0", y = NULL, title = "Proportions of ki67Hi cells: FM") +
  theme(axis.text = element_text(size = 22),
        axis.title =  element_text(size = 20, face = "bold"),
        plot.title = element_text(size=20,  hjust = 0.5),
        legend.text = element_text(size=20),
        legend.title = element_text(size = 20, face = "bold"))

write.csv(counts_GC, "data/counts_LNGC.csv")
write.csv(GC_fd, "data/Nfd_LNGC.csv")
write.csv(ki67_data, "data/ki67_LNGC.csv")





ki67_FM <-  source_data %>%
  select(contains("days"), contains("ki67"))%>%
  select(contains("days"), contains("source"))%>%
  gather(-days.post.bmt, key = subpopulation, value = ki67_fraction)



ggplot(source_data) +
  geom_point(aes(x = age.at.S1K, y = total_counts), col=4, size = 4)+
  scale_y_log10(limits = c(5e6, 2e8), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  scale_x_continuous(limits = c(50, 600), trans = "log10", breaks = c(300, 100, 30)) + 
  theme_bw() + labs(x = "Host age (days)", y = NULL, title = "Total counts: source") +
  theme(axis.text = element_text(size = 22),
        axis.title =  element_text(size = 20, face = "bold"),
        plot.title = element_text(size=20,  hjust = 0.5),
        legend.text = element_text(size=20),
        legend.title = element_text(size = 20, face = "bold"))

ggplot(source_data) +
  geom_point(aes(x = days.post.bmt, y = fd),  size = 4)+
  geom_hline(yintercept = 1.00, linetype = 2, size =1.5, col="darkred")+
  ylim(0, 1.2) + theme_bw() + labs(x = "Days post t0", y = NULL, title = "Noirmalised donor fractions: source") +
  theme(axis.text = element_text(size = 22),
        axis.title =  element_text(size = 20, face = "bold"),
        plot.title = element_text(size=20,  hjust = 0.5),
        legend.text = element_text(size=20),
        legend.title = element_text(size = 20, face = "bold"))


ggplot(ki67_FM)+
  geom_point(aes(x = days.post.bmt, y = ki67_fraction, color = subpopulation), alpha = 0.8, size = 4)+
  scale_color_manual(values = c(4,2), name = "Subset", labels = c("Donor", "Host"))+
  ylim(0, 1) + scale_x_log10(limits = c(10,500), breaks = c(30, 100, 300)) +  theme_bw()+
  labs(x = "Days post t0", y = NULL, title = "Proportions of ki67Hi cells: T2") +
  theme(axis.text = element_text(size = 22),
        axis.title =  element_text(size = 20, face = "bold"),
        plot.title = element_text(size=20,  hjust = 0.5),
        legend.text = element_text(size=20),
        legend.title = element_text(size = 20, face = "bold"))





