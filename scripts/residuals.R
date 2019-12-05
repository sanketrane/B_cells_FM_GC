solve_time <- unique(counts_binned$age.at.S1K)
data_time <- counts_binned$age.at.S1K

time_index <-  purrr::map_dbl(data_time, function(x) which(x == solve_time))

y1res <- filter(Y1pred1, timeseries %in% solve_time) 
y2res <- filter(Y2pred1, timeseries %in% solve_time) 
y3res <- filter(Y3pred1, timeseries %in% solve_time) 
y4res <- filter(Y4pred1, timeseries %in% solve_time) 

ind_y1res <- vector()
ind_y2res <- vector()
ind_y3res <- vector()
ind_y4res <- vector()


for (i in 1:length(counts_binned$total_counts)){
  
  ind_y1res[i]  <- y1res$median[time_index[i]] 
  ind_y2res[i]  <- y2res$median[time_index[i]] 
  ind_y3res[i]  <- y3res$median[time_index[i]] 
  ind_y4res[i]  <- y4res$median[time_index[i]] 
}

res_counts  <-  ind_y1res - counts_binned$total_counts
res_fd      <-  ind_y2res - Nfd_binned$Nfd
res_donorki <-  ind_y3res - filter(ki67_binned, subpopulation == "donor_ki67_FM")$prop_ki67hi
res_hostki  <-  ind_y4res - filter(ki67_binned, subpopulation == "host_ki67_FM")$prop_ki67hi

ggplot()+
  geom_point(aes(x = data_time, y = res_counts)) +
  geom_hline(yintercept = 0, linetype =2, col = "darkred")

ggplot()+
  geom_point(aes(x = data_time, y = res_fd)) +
  geom_hline(yintercept = 0, linetype =2, col = "darkred")

ggplot()+
  geom_point(aes(x = data_time, y = res_donorki)) +
  geom_hline(yintercept = 0, linetype =2, col = "darkred")

ggplot()+
  geom_point(aes(x = data_time, y = res_hostki)) +
  geom_hline(yintercept = 0, linetype =2, col = "darkred")