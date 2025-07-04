########################################################################################################################
# POSTERIOR SAMPLES OF G_PLUS FOR LaPCoM AND THE FOUR COMPETING MODELS ACROSS ALL REPETITIONS
########################################################################################################################

rm(list = ls())

cols = khroma::color("muted")(9)

scenarios = LETTERS[1:8]

G_plus_chain_all_reps_monoLaPCM = G_plus_chain_all_reps_LaPCoM = vector("list", length = 8)

G_plus_df = data.frame(G_monoLaPCM = rep(NA, 8),
                       CI_L_monoLaPCM = rep(NA, 8),
                       CI_U_monoLaPCM = rep(NA, 8),
                       G_LaPCoM = rep(NA, 8),
                       CI_L_LaPCoM = rep(NA, 8),
                       CI_U_LaPCoM = rep(NA, 8))

for (scen_num in 1:length(scenarios)) {
  
  scenario = scenarios[scen_num]
  
  # ====================================================================================================================
  # monoLaPCM  
  # ====================================================================================================================
  load(paste0("/Volumes/My Passport/PhD/PhD Project 1/LaPCoM_Code/Simulation_Study_1/Post_Processing/SS1", scenario, 
              "_PP/pp_res_", scenario, ".Rdata"))
  output = get(paste0("pp_res_", scenario))
  G_plus_chain_all_reps = NA
  for (sim in 1:10) {
    G_plus_chain_all_reps = c(G_plus_chain_all_reps, output$without_node_clustering$chains$G_plus[[sim]])
  } # end sim for loop
  G_plus_chain_all_reps = G_plus_chain_all_reps[-1]
  G_plus_chain_all_reps_monoLaPCM[[scen_num]] = G_plus_chain_all_reps
  G_plus_df$G_monoLaPCM[scen_num] = DescTools::Mode(G_plus_chain_all_reps)[1]
  G_plus_df$CI_L_monoLaPCM[scen_num] = bayestestR::ci(G_plus_chain_all_reps, ci = 0.8)$CI_low
  G_plus_df$CI_U_monoLaPCM[scen_num] = bayestestR::ci(G_plus_chain_all_reps, ci = 0.8)$CI_high
  
  # ====================================================================================================================
  # LaPCoM  
  # ====================================================================================================================
  load(paste0("/Volumes/My Passport/PhD/PhD Project 1/LaPCoM_Code/Simulation_Study_1/Post_Processing/SS1", scenario, 
              "_PP/pp_res_", scenario, ".Rdata"))
  output = get(paste0("pp_res_", scenario))
  G_plus_chain_all_reps = NA
  for (sim in 1:10) {
    G_plus_chain_all_reps = c(G_plus_chain_all_reps, output$with_node_clustering$chains$G_plus[[sim]])
  } # end sim for loop
  G_plus_chain_all_reps = G_plus_chain_all_reps[-1]
  G_plus_chain_all_reps_LaPCoM[[scen_num]] = G_plus_chain_all_reps
  G_plus_df$G_LaPCoM[scen_num] = DescTools::Mode(G_plus_chain_all_reps)[1]
  G_plus_df$CI_L_LaPCoM[scen_num] = bayestestR::ci(G_plus_chain_all_reps, ci = 0.8)$CI_low
  G_plus_df$CI_U_LaPCoM[scen_num] = bayestestR::ci(G_plus_chain_all_reps, ci = 0.8)$CI_high
} # end scen_num for loop

# ----------
round(G_plus_df)

library(ggplot2)
library(dplyr)
library(tidyr)

plot_data = data.frame()

# Loop through scenarios and compute proportions
for (scen_num in seq_along(scenarios)) {
  scenario = scenarios[scen_num]
  
  monoLaPCM_vals = as.vector(G_plus_chain_all_reps_monoLaPCM[[scen_num]])
  LaPCoM_vals = as.vector(G_plus_chain_all_reps_LaPCoM[[scen_num]])
  
  monoLaPCM_props = table(factor(monoLaPCM_vals, levels = 1:8)) / length(monoLaPCM_vals)
  LaPCoM_props = table(factor(LaPCoM_vals, levels = 1:8)) / length(LaPCoM_vals)
  
  df = data.frame(
    G_plus = rep(1:8, times = 2),
    Proportion = c(as.numeric(monoLaPCM_props), as.numeric(LaPCoM_props)),
    Model = rep(c("monoLaPCM", "LaPCoM"), each = 8),
    Scenario = paste("Scenario", scenario)
  )
  
  plot_data = rbind(plot_data, df)
} # end scen_num for loop
rm(scen_num)

plot_data$G_plus = factor(plot_data$G_plus)
plot_data$Model = factor(plot_data$Model, levels = c("monoLaPCM", "LaPCoM"))

model_colors = unname(cols[1:2])

p = ggplot(plot_data, aes(x = G_plus, y = Proportion, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "white") +
  facet_wrap(~ Scenario, nrow = 2) +
  scale_fill_manual(values = model_colors) +
  labs(x = expression(G['+']), y = "Posterior Proportion") +
  theme_minimal(base_size = 16) +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 14),
    panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.8)
  )

# Save the plot as a PDF
ggsave("/Volumes/My Passport/PhD/PhD Project 1/LaPCoM_Code/Simulation_Study_1/Results_Plots/SS1_G_plus_Results.pdf",
       plot = p, width = 16, height = 8)

