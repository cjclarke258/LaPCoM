########################################################################################################################
# PRIMARY SCHOOL (COUNT) DATA APPLICATION (POSTERIOR PREDICTIVE CHECKS)
########################################################################################################################

rm(list = ls())

setwd("/Volumes/My Passport/PhD/PhD Project 1/LaPCoM_Code/Applications/Primary_School_Count/")
source("LaPCoM_Initialisation_Functions.R")

# ======================================================================================================================
# LOAD DATA
# ======================================================================================================================
# ACTUAL DATA
psc_multi = readRDS("primary_school_adj_matrices_count.rds")

# LaPCoM RESULTS
load("pp_res_PSC.Rdata")

# ======================================================================================================================
# EXTRACT THE LAST 100 ALPHAS, TAUS AND LATENT SPACES FOR EACH SIMULATION
# ======================================================================================================================
best_chain_num = 8
num_PPCs = 500
num_samps = dim(primary_school_pp$pp_res_primary_school$chains$Z[[best_chain_num]])[4]

LaPCoM_PPC_alphas = tail(primary_school_pp$pp_res_primary_school$chains$alpha[[best_chain_num]], num_PPCs)
LaPCoM_PPC_Z = primary_school_pp$pp_res_primary_school$chains$Z[[best_chain_num]][, , , tail(1:num_samps, num_PPCs)]
LaPCoM_PPC_C = primary_school_pp$pp_res_primary_school$chains$C[[best_chain_num]][, tail(1:num_samps, num_PPCs)]

# ======================================================================================================================
# GENERATE 100 POSTERIOR PREDICTIVE MULTIPLEXES
# ======================================================================================================================
N = dim(psc_multi)[1]
M = dim(psc_multi)[3]

# multis_PPC = probs = array(0, dim = c(N, N, M, num_PPCs))
# 
# set.seed(783)
# 
# prog_bar = progress::progress_bar$new(format = paste0("Generate Networks: [:bar] :percent"), total = 500, clear = F, width = 60)
# 
# for (ppc_i in 1:num_PPCs) {
# 
#   # ------------------------------------------------------------------------------------------------------------------
#   # GENERATE THE NETWORKS
#   # ------------------------------------------------------------------------------------------------------------------
#   for (m in 1:M) {
#     g = LaPCoM_PPC_C[m, ppc_i] # which cluster network m belongs to
#     D = as.matrix(Rfast::Dist(LaPCoM_PPC_Z[, , g, ppc_i], square = T)) # distance matrix for Z_g
#     eta = LaPCoM_PPC_alphas[ppc_i] - D
#     # probs[, , m, ppc_i] = (exp(eta)) / (1 + exp(eta))
#     probs[, , m, ppc_i][upper.tri(probs[, , m, ppc_i])] = (exp(eta[upper.tri(eta)])) / (1 + exp(eta[upper.tri(eta)]))
#     probs[, , m, ppc_i] = probs[, , m, ppc_i] + t(probs[, , m, ppc_i])
#     diag(probs[, , m, ppc_i]) = 0
#     multis_PPC[, , m, ppc_i][upper.tri(multis_PPC[, , m, ppc_i])] = rpois((N * (N - 1)) / 2, exp(eta))
#     # multis_PPC[, , m, ppc_i][upper.tri(multis_PPC[, , m, ppc_i])] = rbinom(((N * (N - 1)) / 2), 1, probs[, , m, ppc_i][upper.tri(probs[, , m, ppc_i])])
#     multis_PPC[, , m, ppc_i] = multis_PPC[, , m, ppc_i] + t(multis_PPC[, , m, ppc_i])
#     rm(D, g, eta)
#   } # end m for loop
# 
#   prog_bar$tick()
#   rm(m)
# 
# } # end ppc_i for loop
# rm(ppc_i)
# 
# saveRDS(multis_PPC, "PSC_multis_PPCs.RDS")

# ======================================================================================================================
# USE SUMMARY STATISTICS TO COMPARE THE SIMULATED POSTERIOR PREDICTIVE MULTIPLEXES TO THE TRUE MULTIPLEX
# ======================================================================================================================
multis_PPC = readRDS("PSC_multis_PPCs.RDS")
# 
net_dists = matrix(NA, nrow = M, ncol = num_PPCs)
mean_abs_diff_counts = matrix(NA, nrow = M, ncol = num_PPCs)

prog_bar = progress::progress_bar$new(format = paste0("500 PPCs for each network: [:bar] :percent"), total = 16, clear = F, width = 60)

for (m in 1:M) {

  # ----------------------------------------------------------------------------------------------------------------
  # MEAN ABSOLUTE DIFFERENCE IN COUNTS
  # ----------------------------------------------------------------------------------------------------------------
  mean_abs_diff_counts[m, ] = apply(multis_PPC[, , m, ], 3, function(x) { mean(abs(x - psc_multi[, , m])) })

  # ----------------------------------------------------------------------------------------------------------------
  # NETWORK-DISTANCE
  # ----------------------------------------------------------------------------------------------------------------
  gra1 = igraph::graph_from_adjacency_matrix(psc_multi[, , m], mode = "undirected")
  el1 = igraph::as_edgelist(gra1)

  net_dists[m, ] = sapply(seq_len(num_PPCs), function(ppc_i) {
    gra2 = igraph::graph_from_adjacency_matrix(multis_PPC[, , m, ppc_i], mode = "undirected")
    el2 = igraph::as_edgelist(gra2)
    return(nature_dist(el1, el2, gra1, gra2, w1 = 0.45, w2 = 0.45, w3 = 0.1, net_type = "count"))
  })
  rm(gra1, el1)
  
  prog_bar$tick()

} # end m for loop
rm(m)

save.image(file = "LaPCoM_PPCs_PSC.Rdata")

beepr::beep(sound = 10)

# ======================================================================================================================
# PLOTS
# ======================================================================================================================
rm(list = ls())

setwd("/Volumes/My Passport/PhD/PhD Project 1/LaPCoM_Code/Applications/Primary_School_Count/")

load("LaPCoM_PPCs_PSC.Rdata")
psc_multi = readRDS("primary_school_adj_matrices_count.rds")
network_names = dimnames(psc_multi)[[3]]

library(tidyverse)
library(ggplot2)
library(tidyr)
library(dplyr)

# ----------------------------------------------------------------------------------------------------------------------
# MEAN ABSOLUTE DIFFERENCE IN COUNTS
# ----------------------------------------------------------------------------------------------------------------------
round(mean(mean_abs_diff_counts), 2)

df_plot = as.data.frame(mean_abs_diff_counts)
df_plot$Network = factor(network_names, levels = network_names)  # network_names should have length 16
df_long = df_plot %>%
  pivot_longer(cols = -Network, names_to = "Replication", values_to = "MAD")

p = ggplot(df_long, aes(x = 0, y = MAD)) +
  geom_boxplot(fill = "white", color = "black", outlier.shape = 16) +  # Boxplot for each network
  facet_wrap(~ Network, ncol = 4) +  # 4x4 grid
  labs(x = NULL, y = "Mean Absolute Difference") +  # Labels
  theme_minimal(base_size = 16) +  # Minimal theme
  theme(
    strip.text = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.8),
    axis.text.x = element_blank()
  )
ggsave("PPC_Plots_PSC/PSC_PPC_MAD.pdf", p, height = 8, width = 8)

# ----------------------------------------------------------------------------------------------------------------------
# NETWORK DISTANCES
# ----------------------------------------------------------------------------------------------------------------------
round(mean(net_dists), 2)

df_plot = as.data.frame(net_dists)
df_plot$Network = factor(network_names, levels = network_names)  # network_names should have length 16
df_long = df_plot %>%
  pivot_longer(cols = -Network, names_to = "Replication", values_to = "Net_Dist")

p = ggplot(df_long, aes(x = 0, y = Net_Dist)) +
  geom_boxplot(fill = "white", color = "black", outlier.shape = 16) +  # Boxplot for each network
  facet_wrap(~ Network, ncol = 4) +  # 4x4 grid
  labs(x = NULL, y = "Network Distance") +  # Labels
  theme_minimal(base_size = 16) +  # Minimal theme
  theme(
    strip.text = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.8),
    axis.text.x = element_blank()
  )
ggsave("PPC_Plots_PSC/PSC_PPC_Net_Dists.pdf", p, height = 8, width = 8)

# ----------------------------------------------------------------------------------------------------------------------
# LOG FREQUENCY OF COUNTS
# ----------------------------------------------------------------------------------------------------------------------
# df_all = bind_rows(lapply(1:M, function(m) {
#   data_mat = log_freq_counts[[m]][, 1]  # Extract the first column for each network
#   data_frame(Frequency = data_mat, Network = network_names[m])  # Create data frame with Frequency and Network columns
# }))
# 
# # Plot the data with ggplot
# ggplot(df_all, aes(x = 0, y = Frequency)) +
#   geom_boxplot(outlier.shape = NA, fill = "grey80", color = "black") +
#   facet_wrap(~ Network, ncol = 4) +
#   labs(x = "Network", y = "Log Frequency of Zeroes") +
#   theme_minimal(base_size = 16) +
#   theme(
#     strip.text = element_text(size = 18, face = "bold"),
#     axis.text = element_text(size = 14),
#     axis.title = element_text(size = 16),
#     plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
#     panel.grid.major = element_line(color = "gray90"),
#     panel.grid.minor = element_blank(),
#     legend.position = "none",
#     panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.8)
#   )
# 
# # 
# 
# df_all = bind_rows(lapply(1:M, function(m) {
#   data_mat = log_freq_counts[[m]][, 1]  # Extract the first column for each network
#   data_frame(Frequency = data_mat, Network = network_names[m])  # Create data frame with Frequency and Network columns
# }))
# 
# df_truth = bind_rows(lapply(1:M, function(m) {
#   data_frame(Frequency = true_count_freqs_list[[m]], Network = network_names[m])  # Add the truth frequency
# }))
# 
# ggplot(df_all, aes(x = 0, y = Frequency)) +
#   geom_boxplot(outlier.shape = NA, fill = "grey80", color = "black") +  # Boxplot
#   geom_point(data = df_truth, aes(x = 0, y = Frequency), color = "pink", shape = 4, size = 3, linewidth = 2) +  # Pink cross for truth
#   facet_wrap(~ Network, ncol = 4) +
#   labs(x = "Network", y = "Log Frequency of Zeroes") +
#   theme_minimal(base_size = 16) +
#   theme(
#     strip.text = element_text(size = 18, face = "bold"),
#     axis.text = element_text(size = 14),
#     axis.title = element_text(size = 16),
#     plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
#     panel.grid.major = element_line(color = "gray90"),
#     panel.grid.minor = element_blank(),
#     legend.position = "none",
#     panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.8)
#   )

# ----------------------------------------------------------------------------------------------------------------------
# TABLE OF RESULTS (THESE BOXPLOTS ARE SO NARROW, THEY AREN'T REALLY HELPFUL TO USE, VISUALLY)
# ----------------------------------------------------------------------------------------------------------------------
mean_only_stats = function(x) { mean(x) }

mean_mean_abs_diff_counts = apply(mean_abs_diff_counts, 1, mean_only_stats)
mean_net_dists = apply(net_dists, 1, mean_only_stats)

summary_means = data.frame( Network = 1:16,
                            MeanAbsDiffCounts = mean_mean_abs_diff_counts,
                            NetDists = mean_net_dists)

summary_means = round(summary_means, 4)

library(xtable)

xt = xtable(summary_means, 
             caption = "Posterior Predictive Check Means by Network", 
             label = "tab:ppc_means")

addtorow = list()
addtorow$pos = list(0)
addtorow$command = paste0(
  "\\toprule\n",
  "Network & Mean Absolute Diff Counts & Network Distances \\\\\n",
  "\\midrule\n"
)

print(xt, 
      include.rownames = FALSE, 
      booktabs = TRUE,
      add.to.row = addtorow,
      sanitize.colnames.function = identity)

# ----------------------------------------------------------------------------------------------------------------------
# COUNT DISTRIBUTION OF THE ZEROES
# ----------------------------------------------------------------------------------------------------------------------
zero_counts = list()

for (m in 1:M) {
  net_name = network_names[m]
  
  # True network zero count
  true_net = psc_multi[, , m]
  true_zero_count = sum(true_net == 0)
  
  zero_counts[[length(zero_counts) + 1]] = data.frame(
    network = net_name,
    type = "True",
    zero_count = true_zero_count,
    log_zero_count = log1p(true_zero_count),
    ppc_id = NA
  )
  
  # PPC networks
  for (ppc_i in 1:num_PPCs) {
    pred_net = multis_PPC[, , m, ppc_i]
    pred_zero_count = sum(pred_net == 0)
    
    zero_counts[[length(zero_counts) + 1]] = data.frame(
      network = net_name,
      type = "PPC",
      zero_count = pred_zero_count,
      log_zero_count = log1p(pred_zero_count),
      ppc_id = ppc_i
    )
  } # end ppc_i for loop
  rm(ppc_i)
} # end m for loop
rm(m)

zero_df = bind_rows(zero_counts)

library(ggplot2)
library(dplyr)

pdf("PPC_Plots_PSC/PSC_PPC_Zeroes.pdf", width = 15, height = 6)
ggplot(zero_df, aes(x = network, y = log_zero_count)) +
  geom_boxplot(data = filter(zero_df, type == "PPC"),
               fill = "gray85", color = "black", outlier.shape = NA) +
  geom_point(data = filter(zero_df, type == "True"),
             color = "magenta", size = 2) +
  coord_cartesian(ylim = c(10.825, 10.95)) +
  labs(x = "Network", y = "log(1 + # of zero edges)",
       title = "") +
  theme_minimal(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14)
  ) 
dev.off()

# ----------------------------------------------------------------------------------------------------------------------
# CDF PLOTS [COUNTS > 1]
# ----------------------------------------------------------------------------------------------------------------------
# par(mfrow = c(4, 4))
# for (m in 1:M) {
#   
#   true_net = psc_multi[, , m]
#   true_vec = as.vector(true_net)
#   # true_ecdf = ecdf(log(true_net + 1)) # when x = 1, Fn(x) = 0.95, so the plot is very bunch upwards, so consider 0 alone
#   true_ecdf = ecdf(log(true_net[true_net > 0]))
#   
#   plot(true_ecdf, main = network_names[m], verticals = T, pch = 16, cex = 0.5, col = "magenta", lwd = 2)
#   legend("bottomright", legend = "True ECDF", col = "magenta", lty = 1, lwd = 2, bty = "n")
#   
#   for (ppc_i in 1:num_PPCs) {
#     
#     pred_net = multis_PPC[, , m, ppc_i]
#     pred_vec = as.vector(pred_net)
#     pred_ecdf = ecdf(log(pred_net[pred_net > 0]))
#     
#     plot(pred_ecdf, verticals = T, pch = 16, cex = 0.5, col = adjustcolor("black", 0.4), add = T)
#     
#   } # end ppc_i for loop
#   
#   plot(true_ecdf, verticals = T, pch = 16, cex = 0.5, col = "magenta", lwd = 2, add = T) # on top
#   
#   prog_bar$tick()
#   
# } # end m for loop
# rm(m)

library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)

all_ecdf_data = list()

for (m in 1:M) {
  net_name = network_names[m]
  
  # True network
  true_vals = log(as.vector(psc_multi[,,m])[psc_multi[,,m] > 0])
  true_ecdf = ecdf(true_vals)
  true_df = data.frame(
    x = sort(unique(true_vals)),
    y = true_ecdf(sort(unique(true_vals))),
    network = net_name,
    type = "True",
    ppc_id = NA
  )
  all_ecdf_data[[length(all_ecdf_data) + 1]] = true_df
  
  # PPC networks
  for (ppc_i in 1:num_PPCs) {
    pred_vals = log(as.vector(multis_PPC[,,m,ppc_i])[multis_PPC[,,m,ppc_i] > 0])
    if (length(pred_vals) == 0) next  # Skip if no positive entries
    pred_ecdf = ecdf(pred_vals)
    x_vals = sort(unique(pred_vals))
    pred_df = data.frame(
      x = x_vals,
      y = pred_ecdf(x_vals),
      network = net_name,
      type = "PPC",
      ppc_id = ppc_i
    )
    all_ecdf_data[[length(all_ecdf_data) + 1]] = pred_df
  } # end ppc_i network
  rm(ppc_i)
} # end m for loop
rm(m)

ecdf_df = bind_rows(all_ecdf_data)

pdf("PPC_Plots_PSC/PSC_PPC_ECDF.pdf", height = 8, width = 8)
ggplot(ecdf_df, aes(x = x, y = y, group = interaction(type, ppc_id))) +
  geom_step(data = filter(ecdf_df, type == "PPC"),
            color = "black", alpha = 0.2) +
  geom_step(data = filter(ecdf_df, type == "True"),
            color = "magenta", linewidth = 0.75) +
  facet_wrap(~network, ncol = 4) +
  theme_minimal() +
  labs(x = "", y = "ECDF", title = "") +
  theme_minimal(base_size = 16) +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    strip.background = element_rect(fill = "#f0f0f0"),
    strip.text = element_text(face = "bold", size = 12),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 18) 
  ) 
dev.off()
