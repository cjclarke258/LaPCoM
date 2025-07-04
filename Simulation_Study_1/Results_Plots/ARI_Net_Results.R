########################################################################################################################
# MEAN ARI WITH SD FOR monoLaPCM AND LaPCoM ACROSS ALL REPETITIONS
########################################################################################################################

rm(list = ls())

cols = khroma::color("light")(9)

scenarios = LETTERS[1:8]

ARI_monoLaPCM = rep(NA, 10)
ARI_LaPCoM = rep(NA, 10)

ARI_df = data.frame(ARI_mean_monoLaPCM = rep(NA, 8),
                    ARI_sd_monoLaPCM = rep(NA, 8),
                    ARI_mean_LaPCoM = rep(NA, 8),
                    ARI_sd_LaPCoM = rep(NA, 8))

for (scen_num in 1:length(scenarios)) {
  
  scenario = scenarios[scen_num]
  
  load(paste0("/Volumes/My Passport/PhD//PhD Project 1/LaPCoM_Code/Simulation_Study_1/Post_Processing/SS1", scenario, 
              "_PP/pp_res_", scenario, ".Rdata"))
  output_monoLaPCM_PP = output_LaPCoM_PP = get(paste0("pp_res_", scenario)); rm(list = paste0("pp_res_", scenario))
  output_monoLaPCM_PP = output_monoLaPCM_PP$without_node_clustering
  output_LaPCoM_PP = output_LaPCoM_PP$with_node_clustering
  
  # ====================================================================================================================
  # monoLaPCM  
  # ====================================================================================================================
  for (sim in 1:10) {
    load(paste0("/Volumes/My Passport/PhD/PhD Project 1/LaPCoM_Code/Simulation_Study_1/Scripts_Results/SS1", scenario, 
                "_Output/Output_SS1", scenario, "_Simulation_", sim, ".Rdata"))
    output_monoLaPCM_sim = get("output"); rm(output)
    true_C = output_monoLaPCM_sim$sim_multi$true_params$C_vector
    est_C = output_monoLaPCM_PP$posterior_means_modes$C[sim, ]
    ARI_monoLaPCM[sim] = mclust::adjustedRandIndex(true_C, est_C)
    rm(output_monoLaPCM_sim)
  } # end sim for loop
  ARI_df$ARI_mean_monoLaPCM[scen_num] = mean(ARI_monoLaPCM)
  ARI_df$ARI_sd_monoLaPCM[scen_num] = sd(ARI_monoLaPCM)
  cat(paste0("Scenario ", scenario, ": monoLaPCM done\n"))
  
  # ====================================================================================================================
  # LaPCoM  
  # ====================================================================================================================
  for (sim in 1:10) {
    load(paste0("/Volumes/My Passport/PhD/PhD Project 1/LaPCoM_Code/Simulation_Study_1/Scripts_Results/SS1", scenario, 
                "_Output/Output_SS1", scenario, "_Simulation_", sim, ".Rdata"))
    output_LaPCoM_sim = get("output"); rm(output)
    true_C = output_LaPCoM_sim$sim_multi$true_params$C_vector
    est_C = output_LaPCoM_PP$posterior_means_modes$C[sim, ]
    ARI_LaPCoM[sim] = mclust::adjustedRandIndex(true_C, est_C)
    rm(output_LaPCoM_sim)
  } # end sim for loop
  ARI_df$ARI_mean_LaPCoM[scen_num] = mean(ARI_LaPCoM)
  ARI_df$ARI_sd_LaPCoM[scen_num] = sd(ARI_LaPCoM)
  cat(paste0("Scenario ", scenario, ": LaPCoM done\n"))
  
  rm(output_LaPCoM_PP)
  
} # end scen_num for loop

# ----------
round(ARI_df, 2)
