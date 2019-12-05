### model averaging using LOO-CV values

## For FM cells

#TDM
modelName <- "simple_version_TDM_FM2"
outDir <- "/home/sanket/Desktop/Git_repos/Stan_B_cell_modelling/models"
modelDir <- file.path(outDir, modelName)
FM_TDM <- readRDS(file.path(modelDir, paste(modelName, "Fit.rds", sep = "")))


#SHM
modelName <- "simple_model_FM"
outDir <- "/home/sanket/Desktop/Git_repos/Stan_B_cell_modelling/models"
modelDir <- file.path(outDir, modelName)
FM_SHM <- readRDS(file.path(modelDir, paste(modelName, "Fit.rds", sep = "")))

#INC
modelName <- "simple_version_INC_FM"
outDir <- "/home/sanket/Desktop/Git_repos/Stan_B_cell_modelling/models"
modelDir <- file.path(outDir, modelName)
FM_INC <- readRDS(file.path(modelDir, paste(modelName, "Fit.rds", sep = "")))

model_list <- list(FM_TDM, FM_SHM, FM_INC)
log_lik_list <- lapply(model_list, extract_log_lik)

# optional but recommended
r_eff_list <- lapply(model_list, function(x) {
  ll_array <- extract_log_lik(x, merge_chains = FALSE)
  relative_eff(exp(ll_array))
})

# stacking method:
wts1 <- loo_model_weights(
  log_lik_list,
  method = "stacking",
  r_eff_list = r_eff_list,
  optim_control = list(reltol=1e-10)
)
print(wts1)

# can also pass a list of psis_loo objects to avoid recomputing loo
loo_list <- lapply(1:length(log_lik_list), function(j) {
  loo(log_lik_list[[j]], r_eff = r_eff_list[[j]])
})

wts2 <- loo_model_weights(
  loo_list,
  method = "stacking",
  optim_control = list(reltol=1e-10)
)
all.equal(wts1, wts2)


model_wt <- function(x) { exp(-0.5 * x) }
loo_list <- c(130.55, 127.68, 121.5, 118.94, 118.86, 117.83, 124.40, 127.25, 120.35,
              130.57, 131.14, 123.51, 131.77, 130.17, 131.7, 133.4, 128.5) - 117.83    # SPGC
#loo_list <- c(74.86, 77.42, 75.54, 79.24, 75.11, 75.82, 76.79, 81.81, 77.46, 81.11, 74.83, 77.69) - 74.83    ## LNGC
model_list <- cbind(lapply(loo_list, model_wt))

tot_sum <- 0
model_avg <- function(x){
  for (i in 1:length(loo_list)){
    
    tot_sum <- model_list[[i]] + tot_sum}
  
  ans = exp(-0.5 * x) * 100 /tot_sum
  return(ans)}

wt_list <- cbind(lapply(loo_list, model_avg))

tot2_sum = 0
for (i in 1:length(loo_list)){
  
  tot2_sum <- wt_list[[i]] + tot2_sum}

