########################################################################################################################
# MEAN ARI WITH SD FOR LaPCoM AND MANTZIOU'S MODEL ACROSS ALL REPETITIONS (NODE CLUSTERING)
########################################################################################################################

rm(list = ls())

scenarios = c("I", "II", "III", "IV", "V")

ARI_LaPCoM = vector("list", 5)
ARI_Mantziou = vector("list", 5)

for (scen_num in 1:length(scenarios)) {
  
  scenario = scenarios[scen_num]
  
  true_G = c(2, 2, 3, 4, 4)[scen_num]
  ARI_LaPCoM[[scen_num]] = matrix(NA, nrow = true_G, ncol = 10)
  ARI_Mantziou[[scen_num]] = matrix(NA, nrow = true_G, ncol = 10)
  
  # ====================================================================================================================
  # LaPCoM
  # ====================================================================================================================
  load(paste0("/Volumes/My Passport/PhD/PhD Project 1/LaPCoM_Code/Simulation_Study_2/Post_Processing/SS2", scenario,
              "_CJ_PP/pp_res_", scenario, ".Rdata"))
  output_LaPCoM_PP = get(paste0("pp_res_", scenario)); rm(list = paste0("pp_res_", scenario))
  
  for (sim in 1:10) {
    
    load(paste0("/Volumes/My Passport/PhD/PhD Project 1/LaPCoM_Code/Simulation_Study_2/Scripts_Results/SS2_S", scenario, 
                "_CJ/Output_SS2_S", scenario, "_CJ_Simulation_", sim, ".Rdata"))
    output_LaPCoM_sim = get("output"); rm(output)
    
    # ------------------------------------------------------------------------------------------------------------------
    # FIRST, I NEED TO MATCH THE ESTIMATED CLUSTERS / LATENT SPACES TO THE TRUE ONES
    # ------------------------------------------------------------------------------------------------------------------
    true_C = output_LaPCoM_sim$sim_multi$true_params$C_vector
    sim_est_C = output_LaPCoM_PP$posterior_means_modes$C[sim, ]
    
    true_G = length(unique(true_C))
    sim_est_G = length(unique(sim_est_C))
    
    cluster_matching_networks = as.vector(e1071::matchClasses(table(sim_est_C, true_C, dnn = c("Estimated", "True"))))
    
    for (g_true in 1:true_G) {
      g_est_matches = which(cluster_matching_networks == g_true)
      if (length(g_est_matches) == 0) {
        ARI_LaPCoM[[scen_num]][g_true, sim] = NA
      } else {
        ARI_LaPCoM_temp = rep(NA, length(g_est_matches))
        for (g_est_num in 1:length(g_est_matches)) {
          g_est = g_est_matches[g_est_num]
          true_S_g = output_LaPCoM_sim$sim_multi$true_params$S_vec_list[[g_true]]
          sim_est_S_g = output_LaPCoM_PP$posterior_means_modes$S[[sim]][[g_est]]
          if (is.null(sim_est_S_g)) {
            ARI_LaPCoM_temp[g_est_num] = NA
          } else {
            ARI_LaPCoM_temp[g_est_num] = mclust::adjustedRandIndex(true_S_g, sim_est_S_g)
          } # end if else statement
        } # end g_est_num for loop
        ARI_LaPCoM[[scen_num]][g_true, sim] = mean(ARI_LaPCoM_temp)
      } # end if else statement
    } # end g_true for loop
    rm(output_LaPCoM_sim)
    
  } # end sim for loop
  rm(output_LaPCoM_PP)
  
  cat("LaPCoM Scenario", scenario, ": done\n\n")
  
  # ====================================================================================================================
  # Mantziou
  # ====================================================================================================================
  for (sim in 1:10) {
    load(paste0("/Volumes/My Passport/PhD/PhD Project 1/LaPCoM_Code/Simulation_Study_2/Scripts_Results/SS2_S", scenario,
                "_Mantziou/Output_SS2_S", scenario, "_Mantziou_Simulation_", sim, ".Rdata"))
    output_Mantziou_sim = get("output")
    psm_Mantziou = mcclust::comp.psm(output_Mantziou_sim$res_Mantziou$z)
    cluster_order = unique(mcclust.ext::minVI(psm_Mantziou, method = "all", cls.draw = output_Mantziou_sim$res_Mantziou$z, include.greedy = F)$cl[4, ])
    Mantziou_optimal_clustering_networks = cluster_order[mcclust.ext::minVI(psm_Mantziou, method = "all", cls.draw = output_Mantziou_sim$res_Mantziou$z, include.greedy = F)$cl[1, ]]

    cluster_g_matching = as.vector(e1071::matchClasses(table(Mantziou_optimal_clustering_networks, output_Mantziou_sim$sim_multi$true_params$C_vector, dnn = c("Estimated", "True"))))

    for (g_true in 1:true_G) {
      g_est_matches = which(cluster_g_matching == g_true)
      if (length(g_est_matches) == 0) {
        ARI_Mantziou[[scen_num]][g_true, sim] = NA
      } else {
        ARI_Mantziou_temp = rep(NA, length(g_est_matches))
        for (g_est_num in 1:length(g_est_matches)) {
          g_est = g_est_matches[g_est_num]
          true_S = output_Mantziou_sim$sim_multi$true_params$S_vec_list[[g_true]]
          psm_nodes = mcclust::comp.psm(output_Mantziou_sim$res_Mantziou$c[[g_est]])
          cluster_order_nodes = unique(mcclust.ext::minVI(psm_nodes, method = "all", cls.draw = output_Mantziou_sim$res_Mantziou$c[[g_est]], include.greedy = T)$cl[4, ])
          est_S = cluster_order_nodes[mcclust.ext::minVI(psm_nodes, method = "all", cls.draw = output_Mantziou_sim$res_Mantziou$c[[g_est]], include.greedy = T)$cl[1, ]]
          ARI_Mantziou_temp[g_est_num] = ifelse(is.null(est_S), NA, mclust::adjustedRandIndex(true_S, est_S))
        } # end g_est_num
        ARI_Mantziou[[scen_num]][g_true, sim] = mean(ARI_Mantziou_temp)
      } # end if else statement

    } # end g_true for loop

    rm(output_Mantziou_sim)
  } # end sim for loop

  cat("Mantziou Scenario", scenario, ": done\n\n")
  
} # end scen_num for loop

# ----------------------------------------------------------------------------------------------------------------------
# TAKE A LOOK
# ----------------------------------------------------------------------------------------------------------------------
ARI_LaPCoM
ARI_Mantziou

# 

round(min(apply(ARI_LaPCoM_df, 2, median, na.rm = T)), 2)
round(mean(apply(ARI_LaPCoM_df, 2, IQR, na.rm = T)), 2)

round(median(ARI_LaPCoM_df$pp_IV_Z2[c(1:3, 5:10)]), 2)
round(IQR(ARI_LaPCoM_df$pp_IV_Z2[c(1:3, 5:10)]), 2)

# 


round(min(apply(ARI_Mantziou_df, 2, median, na.rm = T)), 2)
round(mean(apply(ARI_Mantziou_df, 2, IQR, na.rm = T)), 2)

round(min(apply(ARI_Mantziou_df, 2, mean, na.rm = T)), 2)
round(min(apply(ARI_Mantziou_df, 2, sd, na.rm = T)), 2)

save.image("/Volumes/My Passport/PhD/PhD Project 1/LaPCoM_Code/Simulation_Study_2/Results_Plots/ARI_Node_Results_Workspace.Rdata")

# ----------------------------------------------------------------------------------------------------------------------
# PLOT THE RESULTS
# ----------------------------------------------------------------------------------------------------------------------
rm(list = ls())

load("/Volumes/My Passport/PhD/PhD Project 1/LaPCoM_Code/Simulation_Study_2/Results_Plots/ARI_Node_Results_Workspace.Rdata")

library(tidyverse)

cols = khroma::color("muted")(9)

ARI_long = ARI_Mantziou_df %>%
  pivot_longer(cols = everything(), names_to = "Condition", values_to = "ARI") %>%
  extract(Condition, into = c("Scenario", "Latent_Space"), regex = "pp_([^_]+)_Z(\\d+)") %>%
  mutate(
    Scenario = factor(Scenario, levels = c("I", "II", "III", "IV", "V")),
    Latent_Space = factor(paste0("Z", Latent_Space), levels = paste0("Z", 1:4))
  )

library(ggplot2)

model_cols = unname(cols[1:4])
pdf("/Volumes/My Passport/PhD/PhD Project 1/LaPCoM_Code/Simulation_Study_2/Results_Plots/SS2_ARI_Nodes_Mantziou.pdf", width = 14, height = 6)
ggplot(ARI_long, aes(x = Latent_Space, y = ARI, fill = Latent_Space)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~ Scenario, nrow = 1, scales = "free_x") +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  scale_x_discrete(
    labels = c(
      Z1 = expression(Z[1]),
      Z2 = expression(Z[2]),
      Z3 = expression(Z[3]),
      Z4 = expression(Z[4])
    )
  ) +
  scale_fill_manual(values = model_cols[1:4]) + 
  labs(x = NULL, y = "ARI", title = "") +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    strip.background = element_rect(fill = "#f0f0f0"),
    strip.text = element_text(face = "bold", size = 18),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 18) 
  ) 
dev.off()

# 

ARI_long <- ARI_LaPCoM_df %>%
  pivot_longer(cols = everything(), names_to = "Condition", values_to = "ARI") %>%
  extract(Condition, into = c("Scenario", "Latent_Space"), regex = "pp_([^_]+)_Z(\\d+)") %>%
  mutate(
    Scenario = factor(Scenario, levels = c("I", "II", "III", "IV", "V")),
    Latent_Space = factor(paste0("Z", Latent_Space), levels = paste0("Z", 1:4))
  )

pdf("/Volumes/My Passport/PhD/PhD Project 1/LaPCoM_Code/Simulation_Study_2/Results_Plots/SS2_ARI_Nodes_LaPCoM.pdf", width = 14, height = 6)
ggplot(ARI_long, aes(x = Latent_Space, y = ARI, fill = Latent_Space)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~ Scenario, nrow = 1, scales = "free_x") +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  scale_x_discrete(
    labels = c(
      Z1 = expression(Z[1]),
      Z2 = expression(Z[2]),
      Z3 = expression(Z[3]),
      Z4 = expression(Z[4])
    )
  ) +
  scale_fill_manual(values = model_cols[1:4]) + 
  labs(x = NULL, y = "ARI", title = "") +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    strip.background = element_rect(fill = "#f0f0f0"),
    strip.text = element_text(face = "bold", size = 18),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 18) 
  ) 
dev.off()