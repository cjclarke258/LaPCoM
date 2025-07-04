########################################################################################################################
# INSPECTION OF THE PROCRUSTES CORRELATIONS OF THE POSTERIOR MEAN LATENT SPACES TO THE TRUE LATENT SPACES
########################################################################################################################

rm(list = ls())

cols = khroma::color("muted")(9)

scenarios = LETTERS[1:8]

PC_monoLaPCM = PC_LaPCoM = array(NA, dim = c(10, 2, 8))

for (scen_num in 1:length(scenarios)) {
  
  scenario = scenarios[scen_num]
  
  load(paste0("/Volumes/My Passport/PhD/PhD Project 1/LaPCoM_Code/Simulation_Study_1/Post_Processing/SS1", scenario, 
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
    
    # ------------------------------------------------------------------------------------------------------------------
    # FIRST, I NEED TO MATCH THE ESTIMATED CLUSTERS / LATENT SPACES TO THE TRUE ONES
    # ------------------------------------------------------------------------------------------------------------------
    true_C = output_monoLaPCM_sim$sim_multi$true_params$C_vector
    sim_est_C = output_monoLaPCM_PP$posterior_means_modes$C[sim, ]
    
    true_G = length(unique(true_C))
    sim_est_G = length(unique(sim_est_C))
    
    cluster_matching_networks = as.vector(e1071::matchClasses(table(sim_est_C, true_C, dnn = c("Estimated", "True"))))
    
    for (g_true in 1:true_G) {
      g_est_matches = which(cluster_matching_networks == g_true)
      if (length(g_est_matches) == 0) {
        PC_monoLaPCM[sim, g_true, scen_num] = NA
      } else {
        PC_monoLaPCM_temp = rep(NA, length(g_est_matches))
        for (g_est_num in 1:length(g_est_matches)) {
          g_est = g_est_matches[g_est_num]
          PC_monoLaPCM_temp[g_est_num] = vegan::protest(output_monoLaPCM_sim$sim_multi$true_params$Z_array[, , g_true],
                                                        output_monoLaPCM_PP$posterior_means_modes$Z[[sim]][[g_est]])$scale
        } # end g_est_num for loop
        PC_monoLaPCM[sim, g_true, scen_num] = mean(PC_monoLaPCM_temp)
      } # end if else statement
    } # end g_true for loop
    rm(output_monoLaPCM_sim)
    
  } # end sim for loop
  rm(output_monoLaPCM_PP)
  
  cat("monoLaPCM Scenario", scenario, ": done\n\n")
  
  # ====================================================================================================================
  # LaPCoM
  # ====================================================================================================================
  for (sim in 1:10) {
    
    load(paste0("/Volumes/My Passport/PhD/PhD Project 1/LaPCoM_Code/Simulation_Study_1/Scripts_Results/SS1", scenario,
                "_Output/Output_SS1", scenario, "_Simulation_", sim, ".Rdata"))
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
        PC_LaPCoM[sim, g_true, scen_num] = NA
      } else {
        PC_LaPCoM_temp = rep(NA, length(g_est_matches))
        for (g_est_num in 1:length(g_est_matches)) {
          g_est = g_est_matches[g_est_num]
          PC_LaPCoM_temp[g_est_num] = vegan::protest(output_LaPCoM_sim$sim_multi$true_params$Z_array[, , g_true],
                                                     output_LaPCoM_PP$posterior_means_modes$Z[[sim]][[g_est]])$scale
        } # end g_est_num for loop
        PC_LaPCoM[sim, g_true, scen_num] = mean(PC_LaPCoM_temp)
      } # end if else statement
    } # end g_true for loop
    rm(output_LaPCoM_sim)
    
  } # end sim for loop
  rm(output_LaPCoM_PP)
  
  cat("LaPCoM Scenario", scenario, ": done\n\n")
  
} # end scen_num for loop

PC_monoLaPCM
PC_LaPCoM

save.image("/Volumes/My Passport/PhD/PhD Project 1/LaPCoM_Code/Simulation_Study_1/Results_Plots/Procrustes_Correlations_Workspace.Rdata")

beepr::beep(sound = 10)

# ====================================================================================================================
# TAKE A LOOK
# ====================================================================================================================
rm(list = ls())

load("/Volumes/My Passport/PhD/PhD Project 1/LaPCoM_Code/Simulation_Study_1/Results_Plots/Procrustes_Correlations_Workspace.Rdata")

round(mean(apply(PC_monoLaPCM, 2, median), na.rm = T), 2)
round(mean(apply(PC_monoLaPCM, 2, IQR, na.rm = T), na.rm = T), 2)

round(mean(apply(PC_LaPCoM, 2, median), na.rm = T), 2)
round(mean(apply(PC_LaPCoM, 2, IQR, na.rm = T), na.rm = T), 2)

# 

median_mat_monoLaPCM = matrix(NA, nrow = 8, ncol = 2)
sd_mat_monoLaPCM = matrix(NA, nrow = 8, ncol = 2)

for (s in 1:8) {
  median_mat_monoLaPCM[s, ] = apply(PC_monoLaPCM[ , , s], 2, median, na.rm = T)
  sd_mat_monoLaPCM[s, ] = apply(PC_monoLaPCM[ , , s], 2, sd, na.rm = T)
}
summary_df_monoLaPCM = data.frame(
  Scenario = paste0("Scenario ", LETTERS[1:8]),
  Median1 = round(median_mat_monoLaPCM[, 1], 2),
  SD1 = round(sd_mat_monoLaPCM[, 1], 2),
  Median2 = round(median_mat_monoLaPCM[, 2], 2),
  SD2 = round(sd_mat_monoLaPCM[, 2], 2)
)
knitr::kable(summary_df_monoLaPCM, format = "latex", digits = 2, booktabs = TRUE)

# 

median_mat_LaPCoM = matrix(NA, nrow = 8, ncol = 2)
sd_mat_LaPCoM = matrix(NA, nrow = 8, ncol = 2)

for (s in 1:8) {
  median_mat_LaPCoM[s, ] = apply(PC_LaPCoM[ , , s], 2, median)
  sd_mat_LaPCoM[s, ] = apply(PC_LaPCoM[ , , s], 2, sd)
}
summary_df_LaPCoM = data.frame(
  Scenario = paste0("Scenario ", LETTERS[1:8]),
  Median1 = round(median_mat_LaPCoM[, 1], 2),
  SD1 = round(sd_mat_LaPCoM[, 1], 2),
  Median2 = round(median_mat_LaPCoM[, 2], 2),
  SD2 = round(sd_mat_LaPCoM[, 2], 2)
)
knitr::kable(summary_df_LaPCoM, format = "latex", digits = 2, booktabs = TRUE)

# ====================================================================================================================
# PLOT
# ====================================================================================================================
library(ggplot2)
library(tidyr)
library(dplyr)

# Colors
cols = khroma::color("muted")(9)

# Convert the arrays into long format data frames
get_long_df = function(pc_array, model_name) {
  df = as.data.frame.table(pc_array, responseName = "Procrustes_Correlation")
  colnames(df) = c("Simulation", "Latent_Space", "Scenario", "Procrustes_Correlation")
  df$Simulation = as.numeric(as.character(df$Simulation))
  df$Latent_Space = factor(df$Latent_Space, labels = c(expression(Z[1]), expression(Z[2])))
  df$Scenario = factor(df$Scenario, labels = LETTERS[1:8])
  df$Model = model_name
  return(df)
}

dimnames(PC_monoLaPCM) = list(
  Simulation = 1:10,
  Latent_Space = 1:2,
  Scenario = LETTERS[1:8]
)

dimnames(PC_LaPCoM) = list(
  Simulation = 1:10,
  Latent_Space = 1:2,
  Scenario = LETTERS[1:8]
)

monoLaPCM_long = get_long_df(PC_monoLaPCM, "monoLaPCM")
LaPCoM_long = get_long_df(PC_LaPCoM, "LaPCoM")

# Combine the data for plotting together if needed
# combined_long = rbind(monoLaPCM_long, LaPCoM_long)

model_cols = unname(cols[1:2])

# Create plot function
plot_procrustes = function(data, model_name, file_name) {
  p = ggplot(data, aes(x = Scenario, y = Procrustes_Correlation, fill = Latent_Space)) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8)) +
    scale_fill_manual(values = model_cols, labels = c(expression(Z[1]), expression(Z[2]))) +
    scale_y_continuous(limits = c(0.85, 1)) +
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
    ) +
    geom_vline(xintercept = seq(1.5, 7.5, 1), linetype = "dotted", color = "grey50") +
    labs(x = "Scenario", y = "Procrustes Correlation", title = model_name)
  
  ggsave(file_name, plot = p, width = 10, height = 5)
}

# Save the plots as PDFs
plot_procrustes(
  monoLaPCM_long, 
  "monoLaPCM", 
  "/Volumes/My Passport/PhD/PhD Project 1/LaPCoM_Code/Simulation_Study_1/Results_Plots/SS1_Procrustes_Correlations_monoLaPCM.pdf"
)
plot_procrustes(
  LaPCoM_long, 
  "LaPCoM", 
  "/Volumes/My Passport/PhD/PhD Project 1/LaPCoM_Code/Simulation_Study_1/Results_Plots/SS1_Procrustes_Correlations_LaPCoM.pdf"
)