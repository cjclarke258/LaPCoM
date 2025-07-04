########################################################################################################################
# CHECK THE RUNTIMES OF SIMULATION STUDY 1
########################################################################################################################
rm(list = ls())

scenarios = LETTERS[1:8]

time_df = data.frame(ARI_mean_monoLaPCM = rep(NA, 8),
                     ARI_sd_monoLaPCM = rep(NA, 8),
                     ARI_mean_LaPCoM = rep(NA, 8),
                     ARI_sd_LaPCoM = rep(NA, 8))

for (scen_num in 1:length(scenarios)) {
  
  scenario = scenarios[scen_num]
  
  # ====================================================================================================================
  # LaPCoM  
  # ====================================================================================================================
  time_monoLaPCM = rep(NA, 10)
  time_LaPCoM = rep(NA, 10)
  for (sim in 1:10) {
    load(paste0("/Users/cjclarke/Desktop/PhD/PhD Project 1/LaPCoM_Code/Simulation_Study_1/Scripts_Results/SS1", scenario, 
    "_Output/Output_SS1", scenario, "_Simulation_", sim, ".Rdata"))
    output_sim = get("output")
    time_monoLaPCM[sim] = as.numeric(output_sim$res_without_node_clusters$time_taken, units = "hours")
    time_LaPCoM[sim] = as.numeric(output_sim$res_with_node_clusters$time_taken, units = "hours")
    rm(output_sim)
  } # end sim for loop
  time_df$ARI_mean_monoLaPCM[scen_num] = mean(time_monoLaPCM)
  time_df$ARI_sd_monoLaPCM[scen_num] = sd(time_monoLaPCM)
  time_df$ARI_mean_LaPCoM[scen_num] = mean(time_LaPCoM)
  time_df$ARI_sd_LaPCoM[scen_num] = sd(time_LaPCoM)
  
} # end scen_num for loop

# 

time_df
round(time_df)