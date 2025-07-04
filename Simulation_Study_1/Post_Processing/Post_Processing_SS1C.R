####################################################################################################
# POST-PROCESSING FOR SS1 SCENARIO C
####################################################################################################

rm(list = ls())
setwd("/Users/cjclarke/Desktop/PhD/PhD Project 1/LaPCoM_code")

# ==================================================================================================
# LOAD THE RELEVANT FUNCTION/SCRIPT
# ==================================================================================================
source("Simulation_Study_1/Scripts_Results/Seeds_Study_1.R")

source("monoLaPCM_PP_MethodA.R")
source("LaPCoM_PP_MethodA.R")

# ==================================================================================================
# RUN THE FUNCTION/SCRIPT
# ==================================================================================================
output_combined = vector("list", length = 10)
output_combined_monoLaPCM = vector("list", length = 10)
output_combined_LaPCoM = vector("list", length = 10)
for (sim in 1:10) {
  load(paste0("/Volumes/My Passport/PhD/PhD Project 1/LaPCoM_code/Simulation_Study_1/Scripts_Results/SS1C_Output/Output_SS1C_Simulation_", sim, ".Rdata"))
  output_combined[[sim]] = output_combined_monoLaPCM[[sim]] = output_combined_LaPCoM[[sim]] = output
  output_combined_monoLaPCM[[sim]]$res_monoLaPCM = output_combined[[sim]]$res_without_node_clusters
  output_combined_LaPCoM[[sim]]$res_LaPCoM = output_combined[[sim]]$res_with_node_clusters
  rm(output)
} # end sim for loop

# ------------------------------------------------------------------------------------------------
# CALCULATE AVERAGE ACCEPTANCE RATES ACROSS THE REPLICATIONS FOR BOTH MODELS
# ------------------------------------------------------------------------------------------------
avg_ARs_monoLaPCM = avg_ARs_LaPCoM = vector("list", length = 0)

avg_AR_alpha_monoLaPCM = avg_AR_e_monoLaPCM = avg_AR_alpha_LaPCoM = avg_AR_e_LaPCoM = rep(NA, 10)
avg_AR_Z_monoLaPCM = avg_AR_Z_LaPCoM = avg_AR_w_LaPCoM = matrix(NA, nrow = 10, ncol = 12)

for (sim in 1:10) {
  avg_AR_alpha_monoLaPCM[sim] = mean(output_combined_monoLaPCM[[sim]]$res_monoLaPCM$store_ARs$acc_rate_alpha)
  AR_Z_temp = rowMeans(output_combined_monoLaPCM[[sim]]$res_monoLaPCM$store_ARs$acc_rate_Z)
  avg_AR_Z_monoLaPCM[sim, ] = c(AR_Z_temp, rep(NA, 12 - length(AR_Z_temp)))
  avg_AR_e_monoLaPCM[sim] = mean(output_combined_monoLaPCM[[sim]]$res_monoLaPCM$store_ARs$acc_rate_e)
  
  avg_AR_alpha_LaPCoM[sim] = mean(output_combined_LaPCoM[[sim]]$res_LaPCoM$store_ARs$acc_rate_alpha)
  AR_Z_temp = rowMeans(output_combined_LaPCoM[[sim]]$res_LaPCoM$store_ARs$acc_rate_Z)
  avg_AR_Z_LaPCoM[sim, ] = c(AR_Z_temp, rep(NA, 12 - length(AR_Z_temp)))
  avg_AR_e_LaPCoM[sim] = mean(output_combined_LaPCoM[[sim]]$res_LaPCoM$store_ARs$acc_rate_e)
  AR_w_temp = rowMeans(output_combined_LaPCoM[[sim]]$res_LaPCoM$store_ARs$acc_rate_w)
  avg_AR_w_LaPCoM[sim, ] = c(AR_w_temp, rep(NA, 12 - length(AR_w_temp)))
} # end sim for loop

avg_ARs_monoLaPCM$alpha = mean(avg_AR_alpha_monoLaPCM)
avg_ARs_monoLaPCM$Z = colMeans(avg_AR_Z_monoLaPCM, na.rm = T)
avg_ARs_monoLaPCM$e = mean(avg_AR_e_monoLaPCM)

avg_ARs_LaPCoM$alpha = mean(avg_AR_alpha_LaPCoM)
avg_ARs_LaPCoM$Z = colMeans(avg_AR_Z_LaPCoM, na.rm = T)
avg_ARs_LaPCoM$e = mean(avg_AR_e_LaPCoM)
avg_ARs_LaPCoM$w = colMeans(avg_AR_w_LaPCoM, na.rm = T)

avg_ARs_monoLaPCM
avg_ARs_LaPCoM

# ------------------------------------------------------------------------------------------------
# SET THE SEED 
# ------------------------------------------------------------------------------------------------
seed = seeds_study_1[23]
set.seed(seed)

# ------------------------------------------------------------------------------------------------
# POST-PROCESSING
# ------------------------------------------------------------------------------------------------
pp_res_C = vector("list", length = 2)
pp_res_C$without_node_clustering = post_process_monoLaPCM(output = output_combined_monoLaPCM, 
                                                       multi = output_combined[[1]]$sim_multi$multi,
                                                       cols = khroma::color("light")(9),
                                                       plots_path = "Simulation_Study_1/Post_Processing/SS1C_PP/WITHOUT_NC")
cols_index = c(1:DescTools::Mode(pp_res_C$without_node_clustering$posterior_means_modes$G_plus_hat)[1])
pp_res_C$with_node_clustering = post_process_LaPCoM(output = output_combined_LaPCoM, 
                                                multi = output_combined[[1]]$sim_multi$multi,
                                                cols = khroma::color("light")(9)[-cols_index],
                                                plots_path = "Simulation_Study_1/Post_Processing/SS1C_PP/WITH_NC")
save(pp_res_C, file = "Simulation_Study_1/Post_Processing/SS1C_PP/pp_res_C.Rdata")  
