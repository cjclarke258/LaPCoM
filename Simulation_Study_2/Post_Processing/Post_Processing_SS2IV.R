####################################################################################################
# POST-PROCESSING FOR SS2 SCENARIO IV
####################################################################################################

rm(list = ls())
setwd("/Users/cjclarke/Desktop/PhD/PhD Project 1/LaPCoM_Code")

# ==================================================================================================
# LOAD THE RELEVANT FUNCTION/SCRIPT
# ==================================================================================================
source("Simulation_Study_2/Scripts_Results/Seeds_Study_2.R")

source("LaPCoM_PP_MethodA.R")

# ==================================================================================================
# RUN THE FUNCTION/SCRIPT
# ==================================================================================================
output_combined = vector("list", length = 10)
for (sim in 1:10) {
  # load(paste0("Simulation_Study_2/Scripts_Results/SS2_SIV_CJ/Output_SS2_SIV_CJ_Simulation_", sim, ".Rdata"))
  load(paste0("/Volumes/My Passport/PhD/PhD Project 1/LaPCoM_Code/Simulation_Study_2/Scripts_Results/SS2_SIV_CJ/Output_SS2_SIV_CJ_Simulation_", sim, ".Rdata"))
  output_combined[[sim]] = output
  output_combined[[sim]]$res_LaPCoM = output_combined[[sim]]$res_CJ
  rm(output)
} # end sim for loop

# ------------------------------------------------------------------------------------------------
# CALCULATE AVERAGE ACCEPTANCE RATES ACROSS THE REPLICATIONS
# ------------------------------------------------------------------------------------------------
avg_ARs = vector("list", length = 0)

avg_AR_alpha = avg_AR_e = rep(NA, 10)
avg_AR_Z = avg_AR_w = matrix(NA, nrow = 10, ncol = 12)

for (sim in 1:10) {
  avg_AR_alpha[sim] = mean(output_combined[[sim]]$res_LaPCoM$store_ARs$acc_rate_alpha)
  temp_Z = rowMeans(output_combined[[sim]]$res_LaPCoM$store_ARs$acc_rate_Z)
  avg_AR_Z[sim, ] = c(temp_Z, rep(NA, 12 - length(temp_Z)))
  avg_AR_e[sim] = mean(output_combined[[sim]]$res_LaPCoM$store_ARs$acc_rate_e)
  temp_w = rowMeans(output_combined[[sim]]$res_LaPCoM$store_ARs$acc_rate_w)
  avg_AR_w[sim, ] = c(temp_w, rep(NA, 12 - length(temp_w)))
} # end sim for loop

avg_ARs$alpha = mean(avg_AR_alpha)
avg_ARs$Z = colMeans(avg_AR_Z, na.rm = T)
avg_ARs$e = mean(avg_AR_e)
avg_ARs$w = colMeans(avg_AR_w, na.rm = T)

avg_ARs

# ------------------------------------------------------------------------------------------------
# SET THE SEED 
# ------------------------------------------------------------------------------------------------
seed = seeds_study_2[1]
set.seed(seed)

# ------------------------------------------------------------------------------------------------
# POST-PROCESSING
# ------------------------------------------------------------------------------------------------
pp_res_IV = post_process_LaPCoM(output = output_combined,
                            multi = output_combined[[1]]$sim_multi$multi,
                            plots_path = "Simulation_Study_2/Post_Processing/SS2IV_CJ_PP/PP_Plots_IV/")
save(pp_res_IV, file = "Simulation_Study_2/Post_Processing/SS2IV_CJ_PP/pp_res_IV.Rdata")  
