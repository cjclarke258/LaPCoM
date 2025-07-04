########################################################################################################################
# MEAN ARI WITH SD FOR LaPCoM AND THE FOUR COMPETING MODELS ACROSS ALL REPETITIONS
########################################################################################################################

rm(list = ls())

cols = khroma::color("light")(9)

scenarios = c("I", "II", "III", "IV", "V")

ARI_LaPCoM = rep(NA, 10)
ARI_Durante = rep(NA, 10)
ARI_Rebafka = rep(NA, 10)
ARI_Mantziou = rep(NA, 10)
ARI_Signorelli = rep(NA, 10)

ARI_df = data.frame(ARI_mean_CJ = rep(NA, 5),
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
  load(paste0("/Users/cjclarke/Desktop/PhD/PhD Project 1/LaPCoM_Code/Simulation_Study_2/Post_Processing/SS2", scenario,
              "_CJ_PP/pp_res_", scenario, ".Rdata"))
  output_LaPCoM_PP = get(paste0("pp_res_", scenario))
  for (sim in 1:10) {
    load(paste0("/Volumes/My Passport/PhD/PhD Project 1/LaPCoM_Code/Simulation_Study_2/Scripts_Results/SS2_S", scenario, 
                "_CJ/Output_SS2_S", scenario, "_CJ_Simulation_", sim, ".Rdata"))
    output_LaPCoM_sim = get("output")
    true_C = output_LaPCoM_sim$sim_multi$true_params$C_vector
    est_C = output_LaPCoM_PP$posterior_means_modes$C[sim, ]
    ARI_LaPCoM[sim] = mclust::adjustedRandIndex(true_C, est_C)
    rm(output_LaPCoM_sim)
  } # end sim for loop
  ARI_df$ARI_mean_CJ[scen_num] = mean(ARI_LaPCoM)
  ARI_df$ARI_sd_CJ[scen_num] = sd(ARI_LaPCoM)
  rm(output_LaPCoM_PP)
  
  # ====================================================================================================================
  # Durante - PopNet
  # ====================================================================================================================
  # for (sim in 1:10) {
  #   load(paste0("/Users/cjclarke/Desktop/PhD/PhD Project 1/LaPCoM_Code/Simulation_Study_2/Scripts_Results/SS2_S", scenario, 
  #               "_Durante/Output_SS2_S", scenario, "_Durante_Simulation_", sim, ".Rdata"))
  #   output_Durante_sim = get("output")
  #   if (length(output_Durante_sim$res_Durante) == 2) {
  #     res_Durante = apply(output_Durante_sim$res_Durante$Group, 2, function(x) {match(x, unique(x))}) # reorder clustering vectors
  #   } else {
  #     res_Durante = apply(output_Durante_sim$res_Durante, 2, function(x) {match(x, unique(x))}) # reorder clustering vectors
  #   } # end if else statement
  #   psm_Durante = mcclust::comp.psm(t(res_Durante))
  #   cluster_order = unique(mcclust.ext::minVI(psm_Durante, method = "all", cls.draw = t(res_Durante), include.greedy = F)$cl[4, ])
  #   Durante_optimal_clustering_networks = cluster_order[mcclust.ext::minVI(psm_Durante, method = "all", cls.draw = t(res_Durante), include.greedy = F)$cl[1, ]]
  #   ARI_Durante[sim] = mclust::adjustedRandIndex(output_Durante_sim$sim_multi$true_params$C_vector, Durante_optimal_clustering_networks)
  #   rm(output_Durante_sim)
  # } # end sim for loop
  # ARI_df$ARI_mean_Durante[scen_num] = mean(ARI_Durante)
  # ARI_df$ARI_sd_Durante[scen_num] = sd(ARI_Durante)
  
  # ====================================================================================================================
  # Rebafka - graphclust
  # ====================================================================================================================
  # for (sim in 1:10) {
  #   load(paste0("/Users/cjclarke/Desktop/PhD/PhD Project 1/LaPCoM_Code/Simulation_Study_2/Scripts_Results/SS2_S", scenario, 
  #               "_Rebafka/Output_SS2_S", scenario, "_Rebafka_Simulation_", sim, ".Rdata"))    
  #   output_Rebafka_sim = get("output")
  #   ARI_Rebafka[sim] = mclust::adjustedRandIndex(output_Rebafka_sim$sim_multi$true_params$C_vector, output_Rebafka_sim$res_Rebafka$res$graphGroups)
  #   rm(output_Rebafka_sim)
  # } # end sim for loop
  # ARI_df$ARI_mean_Rebafka[scen_num] = mean(ARI_Rebafka)
  # ARI_df$ARI_sd_Rebafka[scen_num] = sd(ARI_Rebafka)
  # rm(output)
  
  # ====================================================================================================================
  # Mantziou
  # ====================================================================================================================
  # for (sim in 1:10) {
  #   load(paste0("/Users/cjclarke/Desktop/PhD/PhD Project 1/LaPCoM_Code/Simulation_Study_2/Scripts_Results/SS2_S", scenario,
  #               "_Mantziou/Output_SS2_S", scenario, "_Mantziou_Simulation_", sim, ".Rdata")) 
  #   output_Mantziou_sim = get("output")
  #   psm_Mantziou = mcclust::comp.psm(output_Mantziou_sim$res_Mantziou$z)
  #   cluster_order = unique(mcclust.ext::minVI(psm_Mantziou, method = "all", cls.draw = output_Mantziou_sim$res_Mantziou$z, include.greedy = F)$cl[4, ])
  #   Mantziou_optimal_clustering_networks = cluster_order[mcclust.ext::minVI(psm_Mantziou, method = "all", cls.draw = output_Mantziou_sim$res_Mantziou$z, include.greedy = F)$cl[1, ]]
  #   ARI_Mantziou[sim] = mclust::adjustedRandIndex(output_Mantziou_sim$sim_multi$true_params$C_vector, Mantziou_optimal_clustering_networks)
  #   rm(output_Mantziou_sim)
  # } # end sim for loop 
  # ARI_df$ARI_mean_Mantziou[scen_num] = mean(ARI_Mantziou)
  # ARI_df$ARI_sd_Mantziou[scen_num] = sd(ARI_Mantziou)

  # ====================================================================================================================
  # Signorelli and Wit 
  # ====================================================================================================================
  # for (sim in 1:10) {
  #   load(paste0("/Users/cjclarke/Desktop/PhD/PhD Project 1/LaPCoM_Code/Simulation_Study_2/Scripts_Results/SS2_S", scenario, 
  #               "_Signorelli/Output_SS2_S", scenario, "_Signorelli_Simulation_", sim, ".Rdata"))
  #   output_SignorelliWit_sim = get("output")
  #   ARI_Signorelli[sim] = mclust::adjustedRandIndex(output_SignorelliWit_sim$sim_multi$true_params$C_vector, output_SignorelliWit_sim$res_Signorelli$clusters_Signorelli)
  #   rm(output_SignorelliWit_sim)
  # } # end sim for loop (Signorelli)
  # ARI_df$ARI_mean_Signorelli[scen_num] = mean(ARI_Signorelli)
  # ARI_df$ARI_sd_Signorelli[scen_num] = sd(ARI_Signorelli)

} # end scen_num for loop

# ----------
round(ARI_df, 2)