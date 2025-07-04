########################################################################################################################
# PRIMARY SCHOOL (COUNT) DATA APPLICATION (POSTERIOR PREDICTIVE CHECKS) [USING LATENTNET AS A COMPARISON]
########################################################################################################################

rm(list = ls())

setwd("/Volumes/My Passport/PhD/PhD Project 1/LaPCoM_Code/Applications/Primary_School_Count/")

# ======================================================================================================================
# LOAD DATA
# ======================================================================================================================
set.seed(123)

psc_multi = readRDS("primary_school_adj_matrices_count.rds")

# ----------------------------------------------------------------------------------------------------------------------    
# SET UP
# ----------------------------------------------------------------------------------------------------------------------
output_psc_ppc_latentnet = list()

N = dim(psc_multi)[1]
M = dim(psc_multi)[3]
G_max = 9

# no_cores = M
no_cores = 2

bic_list = vector("list", length = M)

doParallel::registerDoParallel(cores = no_cores)
cl = parallel::makeCluster(no_cores, type = "FORK")

`%dopar%` = foreach::`%dopar%`

tictoc::tic()

output_psc_ppc_latentnet = foreach::foreach(m = c(10, 11)) %dopar% {

  # ````````````````````````````````````````````````````````````````````````````````````````````````````````````````````
  # CREATE A LOG FILE
  # ````````````````````````````````````````````````````````````````````````````````````````````````````````````````````
  log_file = paste0("PSC_PPC_latentnet_Logs/Updated_PSC_PPC_latentnet_", m, ".txt")
  file.create(log_file)
  
  bic_list[[m]] = rep(NA, G_max)
  
  for (g in 1:G_max) {
    res_lpcm = tryCatch(
      latentnet::ergmm(network::as.network(psc_multi[, , m]) ~ euclidean(d = 2, G = g),
                       control = latentnet::ergmm.control(burnin = 5000, sample.size = 2000, interval = 5)),
      error = function(e) {
        write(paste0("g = ", g, ": Error - ", e$message, "\n"), file = log_file, append = TRUE)
        return(NULL)
      } # end error function
    ) # end tryCatch
    
    if (!is.null(res_lpcm)) {
      sum_lpcm = summary(res_lpcm)
      bic_list[[m]][g] = sum_lpcm$bic$overall
      write(paste0("BIC[", g, "] = ", bic_list[[m]][g], ": Done\n"), file = log_file, append = TRUE)
    } else {
      bic_list[[m]][g] = NA
    } # end if else statement
  } # end g for loop
  rm(g)
  
  # ====================================================================================================================
  # FIND BEST FITTING LATENTNET MODEL FOR EACH NETWORK
  # ====================================================================================================================
  K_opt = which.min(bic_list[[m]])
  write(paste0("K_opt obtained: ", K_opt, "\n"), file = log_file, append = TRUE)
  
  # ====================================================================================================================
  # REFIT EACH BEST ONE AND EXTRACT LAST 100 ALPHAS AND ZS
  # ====================================================================================================================
  best_fit = latentnet::ergmm(network::as.network(psc_multi[, , m]) ~ euclidean(d = 2, G = K_opt),
                               control = latentnet::ergmm.control(burnin = 5000, sample.size = 2000, interval = 5))
  write(paste0("Best fit obtained!\n"), file = log_file, append = TRUE)
  
  output = list(bic_list = bic_list[[m]], K_opt = K_opt, best_fit = best_fit)

  saveRDS(output, file = paste0("PSC_PPC_latentnet_Output/Updated_PSC_PPC_latentnet_", m, ".RDS"))
  
} # end foreach loop

parallel::stopCluster(cl)

tictoc::toc()