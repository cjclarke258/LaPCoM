########################################################################################################################
# CHECK THE RUNTIMES OF SIMULATION STUDY 2
########################################################################################################################
rm(list = ls())

scenarios = c("I", "II", "III", "IV", "V")

time_df = data.frame(ARI_mean_CJ = rep(NA, 5),
                     ARI_sd_CJ = rep(NA, 5),
                     ARI_mean_Durante = rep(NA, 5),
                     ARI_sd_Durante = rep(NA, 5),
                     ARI_mean_Rebafka = rep(NA, 5),
                     ARI_sd_Rebafka = rep(NA, 5),
                     ARI_mean_Mantziou = rep(NA, 5),
                     ARI_sd_Mantziou = rep(NA, 5),
                     ARI_mean_Signorelli = rep(NA, 5),
                     ARI_sd_Signorelli = rep(NA, 5))

for (scen_num in 1:length(scenarios)) {
  
  scenario = scenarios[scen_num]
  
  # ====================================================================================================================
  # LaPCoM  
  # ====================================================================================================================
  time_LaPCoM = rep(NA, 10)
  for (sim in 1:10) {
    load(paste0("/Users/cjclarke/Desktop/PhD/PhD Project 1/LaPCoM_Code/Simulation_Study_2/Scripts_Results/SS2_S", scenario, 
                "_CJ/Output_SS2_S", scenario, "_CJ_Simulation_", sim, ".Rdata"))
    output_LaPCoM_sim = get("output")
    time_LaPCoM[sim] = as.numeric(output_LaPCoM_sim$res_CJ$time_taken, units = "hours")
    rm(output_LaPCoM_sim)
  } # end sim for loop
  time_df$ARI_mean_CJ[scen_num] = mean(time_LaPCoM)
  time_df$ARI_sd_CJ[scen_num] = sd(time_LaPCoM)

  # ====================================================================================================================
  # Durante - PopNet
  # ====================================================================================================================
  time_Durante = rep(NA, 10)
  for (sim in 1:10) {
    load(paste0("/Users/cjclarke/Desktop/PhD/PhD Project 1/LaPCoM_Code/Simulation_Study_2/Scripts_Results/SS2_S", scenario, 
                "_Durante/Output_SS2_S", scenario, "_Durante_Simulation_", sim, ".Rdata"))
    output_Durante_sim = get("output")
    if (length(output_Durante_sim$res_Durante) == 2) {
      time_Durante[sim] = as.numeric(output_Durante_sim$res_Durante$time_taken, units = "hours")
    } else {
      time_Durante[sim] = NA
    } # end if else statement
    rm(output_Durante_sim)
  } # end sim for loop
  time_df$ARI_mean_Durante[scen_num] = mean(time_Durante)
  time_df$ARI_sd_Durante[scen_num] = sd(time_Durante)
  
  # ====================================================================================================================
  # Rebafka - graphclust
  # ====================================================================================================================
  time_Rebafka = rep(NA, 10)
  for (sim in 1:10) {
    load(paste0("/Users/cjclarke/Desktop/PhD/PhD Project 1/LaPCoM_Code/Simulation_Study_2/Scripts_Results/SS2_S", scenario, 
                "_Rebafka/Output_SS2_S", scenario, "_Rebafka_Simulation_", sim, ".Rdata"))    
    output_Rebafka_sim = get("output")
    time_Rebafka[sim] = as.numeric(output_Rebafka_sim$res_Rebafka$time_taken, units = "hours")
    rm(output_Rebafka_sim)
  } # end sim for loop
  time_df$ARI_mean_Rebafka[scen_num] = mean(time_Rebafka)
  time_df$ARI_sd_Rebafka[scen_num] = sd(time_Rebafka)
  rm(output)
  
  # ====================================================================================================================
  # Mantziou
  # ====================================================================================================================
  time_Mantziou = rep(NA, 10)
  for (sim in 1:10) {
    load(paste0("/Users/cjclarke/Desktop/PhD/PhD Project 1/LaPCoM_Code/Simulation_Study_2/Scripts_Results/SS2_S", scenario,
                "_Mantziou/Output_SS2_S", scenario, "_Mantziou_Simulation_", sim, ".Rdata")) 
    output_Mantziou_sim = get("output")
    time_Mantziou[sim] = as.numeric(output_Mantziou_sim$res_Mantziou$time_taken, units = "hours")
    rm(output_Mantziou_sim)
  } # end sim for loop 
  time_df$ARI_mean_Mantziou[scen_num] = mean(time_Mantziou)
  time_df$ARI_sd_Mantziou[scen_num] = sd(time_Mantziou)
  
  # ====================================================================================================================
  # Signorelli and Wit 
  # ====================================================================================================================
  time_Signorelli = rep(NA, 10)
  for (sim in 1:10) {
    load(paste0("/Users/cjclarke/Desktop/PhD/PhD Project 1/LaPCoM_Code/Simulation_Study_2/Scripts_Results/SS2_S", scenario, 
                "_Signorelli/Output_SS2_S", scenario, "_Signorelli_Simulation_", sim, ".Rdata"))
    output_SignorelliWit_sim = get("output")
    time_Signorelli[sim] = as.numeric(output_SignorelliWit_sim$res_Signorelli$time_taken, units = "hours")
    rm(output_SignorelliWit_sim)
  } # end sim for loop
  time_df$ARI_mean_Signorelli[scen_num] = mean(time_Signorelli)
  time_df$ARI_sd_Signorelli[scen_num] = sd(time_Signorelli)

} # end scen_num for loop

# 

time_df
round(time_df)