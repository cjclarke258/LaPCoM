########################################################################################################################
# MEAN ARI WITH SD FOR LaPCoM AND MANTZIOU'S MODEL ACROSS ALL REPETITIONS (NODE CLUSTERING)
########################################################################################################################

rm(list = ls())

cols = khroma::color("muted")(9)

scenarios = LETTERS[1:8]

ARI_LaPCoM = vector("list", 8)

for (scen_num in 1:length(scenarios)) {
  
  scenario = scenarios[scen_num]
  
  true_G = 2
  ARI_LaPCoM[[scen_num]] = matrix(NA, nrow = true_G, ncol = 10)

  # ====================================================================================================================
  # LaPCoM
  # ====================================================================================================================
  load(paste0("/Volumes/My Passport/PhD/PhD Project 1/LaPCoM_Code/Simulation_Study_1/Post_Processing/SS1", scenario, 
              "_PP/pp_res_", scenario, ".Rdata"))
  output_LaPCoM_PP = get(paste0("pp_res_", scenario)); rm(list = paste0("pp_res_", scenario))
  
  for (sim in 1:10) {
    
    load(paste0("/Volumes/My Passport/PhD/PhD Project 1/LaPCoM_Code/Simulation_Study_1/Scripts_Results/SS1", scenario, 
                "_Output/Output_SS1", scenario, "_Simulation_", sim, ".Rdata"))
    output_LaPCoM_sim = get("output"); rm(output)
    
    # ------------------------------------------------------------------------------------------------------------------
    # FIRST, I NEED TO MATCH THE ESTIMATED CLUSTERS / LATENT SPACES TO THE TRUE ONES
    # ------------------------------------------------------------------------------------------------------------------
    true_C = output_LaPCoM_sim$sim_multi$true_params$C_vector
    sim_est_C = output_LaPCoM_PP$with_node_clustering$posterior_means_modes$C[sim, ]
    
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
          sim_est_S_g = output_LaPCoM_PP$with_node_clustering$posterior_means_modes$S[[sim]][[g_est]]
          ARI_LaPCoM_temp[g_est_num] = mclust::adjustedRandIndex(true_S_g, sim_est_S_g)
        } # end g_est_num for loop
        ARI_LaPCoM[[scen_num]][g_true, sim] = mean(ARI_LaPCoM_temp)
      } # end if else statement
    } # end g_true for loop
    rm(output_LaPCoM_sim)
    
  } # end sim for loop
  rm(output_LaPCoM_PP)
  
  cat("LaPCoM Scenario", scenario, ": done\n\n")
  
} # end scen_num for loop

save.image("/Users/cjclarke/Desktop/PhD/PhD Project 1/LaPCoM_Code/Simulation_Study_1/Results_Plots/SS1_ARI_Node_Results_Workspace.Rdata")

# ----------------------------------------------------------------------------------------------------------------------
# TAKE A LOOK
# ----------------------------------------------------------------------------------------------------------------------
ARI_LaPCoM

beepr::beep(sound = 10)

min(apply(ARI_LaPCoM, 2, median), na.rm = T)
max(apply(ARI_LaPCoM, 2, IQR, na.rm = T), na.rm = T)

ARI_LaPCoM[, c(1, 2, 4, 5)]

ARI_LaPCoM[, c(7, 8, 10, 11)]
apply(ARI_LaPCoM[, c(7, 8, 10, 11)], 2, median, na.rm = T)
apply(ARI_LaPCoM[, c(7, 8, 10, 11)], 2, IQR, na.rm = T)

apply(ARI_LaPCoM[, c(19, 20, 22, 23)], 2, median, na.rm = T)
apply(ARI_LaPCoM[, c(19, 20, 22, 23)], 2, IQR, na.rm = T)

# ----------------------------------------------------------------------------------------------------------------------
# PLOT THE RESULTS
# ----------------------------------------------------------------------------------------------------------------------
rm(list = ls())

load("/Users/cjclarke/Desktop/PhD/PhD Project 1/LaPCoM_Code/Simulation_Study_1/Results_Plots/SS1_ARI_Node_Results_Workspace.Rdata")

library(tidyr)
library(ggplot2)
library(dplyr)

# Convert ARI_LaPCoM (list of 8 matrices) into long-format ARI_df
ARI_df = do.call(rbind, lapply(seq_along(ARI_LaPCoM), function(i) {
  mat = ARI_LaPCoM[[i]]
  df = data.frame(
    Scenario = paste0("S", i),
    Z1 = mat[1, ],
    Z2 = mat[2, ],
    Simulation = 1:ncol(mat)
  )
  return(df)
}))

ARI_long = ARI_df %>%
  pivot_longer(
    cols = c(Z1, Z2),
    names_to = "Latent_Space",
    values_to = "ARI"
  ) %>%
  mutate(
    Scenario = factor(Scenario, levels = paste0("S", 1:8), labels = LETTERS[1:8]),
    Latent_Space = factor(Latent_Space, levels = c("Z1", "Z2"), labels = c(expression(Z[1]), expression(Z[2])))
  )

model_cols = unname(cols[1:2])

# Create plot function with the same theme
plot_ARI = function(data, file_name) {
  p = ggplot(data, aes(x = Scenario, y = ARI, fill = Latent_Space)) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8)) +
    scale_fill_manual(values = model_cols[1:2], labels = c(expression(Z[1]), expression(Z[2]))) +
    scale_y_continuous(limits = c(0.7, 1)) +
    # scale_x_discrete(
    #   breaks = seq(1.5, 15.5, 2),  # Adjust the number of breaks based on your data
    #   labels = LETTERS[1:8]
    # ) +
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
      plot.title = element_text(hjust = 0.5, size = 18)  # Center the title
    ) +
    geom_vline(xintercept = seq(1.5, 7.5, 1), linetype = "dotted", color = "grey50") +
    labs(x = "Scenario", y = "ARI", title = "")
  
  ggsave(file_name, plot = p, width = 10, height = 5)
}

plot_ARI(
  ARI_long, 
  "/Volumes/My Passport/PhD/PhD Project 1/LaPCoM_Code/Simulation_Study_1/Results_Plots/SS1_ARI_Node_Results.pdf"
)
