########################################################################################################################
# POSTERIOR SAMPLES OF G_PLUS FOR LaPCoM AND THE FOUR COMPETING MODELS ACROSS ALL REPETITIONS
########################################################################################################################

rm(list = ls())

cols = khroma::color("light")(9)

scenarios = c("I", "II", "III", "IV", "V")

G_LaPCoM_chain = rep(NA, 7000)
G_Durante_chain = rep(NA, 7000)
G_Rebafka_chain = rep(NA, 10)
G_Mantziou_chain = rep(NA, 7000)
G_Signorelli_chain = rep(NA, 10)

G_plus_df = data.frame(G_CJ = rep(NA, 5),
                       CI_L_CJ = rep(NA, 5),
                       CI_U_CJ = rep(NA, 5),
                       G_Durante = rep(NA, 5),
                       CI_L_Durante = rep(NA, 5),
                       CI_U_Durante = rep(NA, 5),
                       G_Rebafka = rep(NA, 5),
                       CI_L_Rebafka = rep(NA, 5),
                       CI_U_Rebafka = rep(NA, 5),
                       G_Mantziou = rep(NA, 5),
                       CI_L_Mantziou = rep(NA, 5),
                       CI_U_Mantziou = rep(NA, 5),
                       G_Signorelli = rep(NA, 5),
                       CI_L_Signorelli = rep(NA, 5),
                       CI_U_Signorelli = rep(NA, 5))

for (scen_num in 1:length(scenarios)) {
  
  scenario = scenarios[scen_num]
  
  # ====================================================================================================================
  # LaPCoM  
  # ====================================================================================================================
  load(paste0("/Users/cjclarke/Desktop/PhD/PhD Project 1/LaPCoM_Code/Simulation_Study_2/Post_Processing/SS2", scenario,
              "_CJ_PP/pp_res_", scenario, ".Rdata"))
  output = get(paste0("pp_res_", scenario))
  G_plus_chain_all_reps = NA
  for (sim in 1:10) {
    G_plus_chain_all_reps = c(G_plus_chain_all_reps, output$chains$G_plus[[sim]])
  } # end sim for loop
  G_plus_chain_all_reps = G_plus_chain_all_reps[-1]
  G_plus_df$G_CJ[scen_num] = DescTools::Mode(G_plus_chain_all_reps)[1]
  G_plus_df$CI_L_CJ[scen_num] = bayestestR::ci(G_plus_chain_all_reps, ci = 0.8)$CI_low
  G_plus_df$CI_U_CJ[scen_num] = bayestestR::ci(G_plus_chain_all_reps, ci = 0.8)$CI_high
  
  # ====================================================================================================================
  # Durante - PopNet
  # ====================================================================================================================
  # for (sim in 1:10) {
  #   load(paste0("/Users/cjclarke/Desktop/PhD/PhD Project 1/LaPCoM_Code/Simulation_Study_2/Scripts_Results/SS2_S", scenario, 
  #               "_Durante/Output_SS2_S", scenario, "_Durante_Simulation_", sim, ".Rdata"))
  #   indices = matrix(NA, nrow = 10, ncol = 2)
  #   start = 1
  #   for (i in 1:10) {
  #     end = start + 700 - 1
  #     indices[i, ] = c(start, end)
  #     start = end + 1
  #   } # end i for loop
  #   if (length(output$res_Durante) == 2) {
  #     G_Durante_chain[indices[sim, 1]:indices[sim, 2]] = apply(output$res_Durante$Group, 2, function(x) {length(unique(x))})
  #   } else {
  #     G_Durante_chain[indices[sim, 1]:indices[sim, 2]] = apply(output$res_Durante, 2, function(x) {length(unique(x))})
  #   } # end if else statement
  #   rm(output)
  # } # end sim for loop 
  # G_plus_df$G_Durante[scen_num] = DescTools::Mode(G_Durante_chain)[1]
  # G_plus_df$CI_L_Durante[scen_num] = bayestestR::ci(G_Durante_chain, ci = 0.8)$CI_low
  # G_plus_df$CI_U_Durante[scen_num] = bayestestR::ci(G_Durante_chain, ci = 0.8)$CI_high
  
  # ====================================================================================================================
  # Rebafka - graphclust
  # ====================================================================================================================
  # for (sim in 1:10) {
  #   load(paste0("/Users/cjclarke/Desktop/PhD/PhD Project 1/LaPCoM_Code/Simulation_Study_2/Scripts_Results/SS2_S", scenario, 
  #               "_Rebafka/Output_SS2_S", scenario, "_Rebafka_Simulation_", sim, ".Rdata"))
  #   G_Rebafka_chain[sim] = length(unique(output$res_Rebafka$res$graphGroups))
  #   rm(output)
  # } # end sim for loop 
  # G_plus_df$G_Rebafka[scen_num] = DescTools::Mode(G_Rebafka_chain)[1]
  # G_plus_df$CI_L_Rebafka[scen_num] = bayestestR::ci(G_Rebafka_chain, ci = 0.8)$CI_low
  # G_plus_df$CI_U_Rebafka[scen_num] = bayestestR::ci(G_Rebafka_chain, ci = 0.8)$CI_high

  # ====================================================================================================================
  # Mantziou
  # ====================================================================================================================
  # for (sim in 1:10) {
  #   load(paste0("/Users/cjclarke/Desktop/PhD/PhD Project 1/LaPCoM_Code/Simulation_Study_2/Scripts_Results/SS2_S", scenario,
  #               "_Mantziou/Output_SS2_S", scenario, "_Mantziou_Simulation_", sim, ".Rdata"))
  #   indices = matrix(NA, nrow = 10, ncol = 2)
  #   start = 1
  #   for (i in 1:10) {
  #     end = start + 700 - 1
  #     indices[i, ] = c(start, end)
  #     start = end + 1
  #   } # end i for loop
  #   G_Mantziou_chain[indices[sim, 1]:indices[sim, 2]] = output$res_Mantziou$C_plus
  #   rm(output)
  # } # end sim for loop 
  # G_plus_df$G_Mantziou[scen_num] = DescTools::Mode(G_Mantziou_chain)[1]
  # G_plus_df$CI_L_Mantziou[scen_num] = bayestestR::ci(G_Mantziou_chain, ci = 0.8)$CI_low
  # G_plus_df$CI_U_Mantziou[scen_num] = bayestestR::ci(G_Mantziou_chain, ci = 0.8)$CI_high
  
  # ====================================================================================================================
  # Signorelli and Wit 
  # ====================================================================================================================
  # for (sim in 1:10) {
  #   load(paste0("/Users/cjclarke/Desktop/PhD/PhD Project 1/LaPCoM_Code/Simulation_Study_2/Scripts_Results/SS2_S", scenario, 
  #               "_Signorelli/Output_SS2_S", scenario, "_Signorelli_Simulation_", sim, ".Rdata"))
  #   G_Signorelli_chain[sim] = length(unique(output$res_Signorelli$clusters_Signorelli))
  #   rm(output)
  # } # end sim for loop 
  # G_plus_df$G_Signorelli[scen_num] = DescTools::Mode(G_Signorelli_chain)[1]
  # G_plus_df$CI_L_Signorelli[scen_num] = bayestestR::ci(G_Signorelli_chain, ci = 0.8)$CI_low
  # G_plus_df$CI_U_Signorelli[scen_num] = bayestestR::ci(G_Signorelli_chain, ci = 0.8)$CI_high

} # end scen_num for loop

# ----------
round(G_plus_df)