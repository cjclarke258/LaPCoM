########################################################################################################################
# PRIMARY SCHOOL (COUNT) DATA APPLICATION (POSTERIOR PREDICTIVE CHECKS) [USING LATENTNET AS A COMPARISON]
########################################################################################################################

rm(list = ls())

setwd("/Volumes/My Passport/PhD/PhD Project 1/LaPCoM_Code/Applications/Primary_School_Count/")
source("LaPCoM_Initialisation_Functions.R")

# ======================================================================================================================
# LOAD DATA
# ======================================================================================================================
psc_multi = readRDS("primary_school_adj_matrices_count.rds")

# ======================================================================================================================
# USE LATENTNET ON EACH NETWORK WITH G=1,...,15 USING BIC TO SELECT THE BEST ONE
# ======================================================================================================================
set.seed(7836867)
N = dim(psc_multi)[1]
M = dim(psc_multi)[3]
# G_max = 9
# bic_mat = matrix(NA, nrow = M, ncol = G_max)
# for (m in 1:M) {
#   cat("\nNetwork", m, "\n---------\n")
#   for (g in 1:G_max) {
#     res_lpcm = latentnet::ergmm(network::as.network(psc_multi[, , m]) ~ euclidean(d = 2, G = g),
#                                 control = latentnet::ergmm.control(burnin = 5000, sample.size = 5000, interval = 50))
#     sum_lpcm = summary(res_lpcm)
#     bic_mat[m, g] = sum_lpcm$bic$overall
#     cat("g =", g, ": Done\n")
#   } # end g for loop
# } # end m for loop
# 
# save(bic_mat, file = "psc_PPC_latentnet_res.Rdata")

# ======================================================================================================================
# FIND BEST FITTING LATENTNET MODEL FOR EACH NETWORK
# ======================================================================================================================
# K_opt = apply(bic_mat, 1, which.min)
# K_opt

# ======================================================================================================================
# REFIT EACH BEST ONE AND EXTRACT LAST 100 ALPHAS AND ZS
# ======================================================================================================================
# best_fits = vector("list", length = M)
# for (m in 1:M) {
#   best_fits[[m]] = latentnet::ergmm(network::as.network(psc_multi[, , m]) ~ euclidean(d = 2, G = K_opt[m]),
#                                     control = latentnet::ergmm.control(burnin = 5000, sample.size = 5000, interval = 50))
#   cat("\nNetwork", m, ": Done")
# } # end m for loop

# FOR THIS APPLICATION, I HAD TO RUN THE ABOVE CODE SEPARATELY ON SONIC

best_fits = vector("list", length = M)
K_opt = rep(NA, M)

for (m in 1:M) {
  
  base_name = paste0("PSC_PPC_latentnet_Output/Updated_PSC_PPC_latentnet_", m)
  rds_file = paste0(base_name, ".RDS")
  rdata_file = paste0(base_name, ".Rdata")
  
  if (file.exists(rds_file)) {
    res = readRDS(rds_file)
    message("Read ", rds_file)
  } else if (file.exists(rdata_file)) {
    load(rdata_file)
    res = output
    message("Loaded ", rdata_file)
  } else {
    warning("No file found for ", base_name)
  } # end if else statement
  
  best_fits[[m]] = res$best_fit
  K_opt[m] = res$K_opt
  
} # end m for loop

latentnet_res = list(best_fits = best_fits, K_opt = K_opt)
save(latentnet_res, file = "psc_PPC_latentnet_res.Rdata")

beepr::beep(sound = 10)

# ======================================================================================================================
# PLOT ALL THE LATENT SPACES
# ======================================================================================================================
rm(list = ls())

setwd("/Volumes/My Passport/PhD/PhD Project 1/LaPCoM_Code/Applications/Primary_School_Count/")

psc_multi = readRDS("primary_school_adj_matrices_count.rds")
M = dim(psc_multi)[3]
N = dim(psc_multi)[1]

load("psc_PPC_latentnet_res.Rdata")

pdf("PPC_Plots_PSC/psc_latentnet_Latent_Spaces.pdf", height = 8, width = 8)
par(mfrow = c(4, 4))
par(mar = c(1.75, 1.75, 1.5, 0.5) + 0.4) # bltr
for (m in 1:M) {
  plot(latentnet_res$best_fits[[m]], main = dimnames(psc_multi)[[3]][m], las = 1, print.formula = F, cex.main = 1.8, vertex.cex = 1.75,
       cluster.col = khroma::color("light")(9))
  legend("topright", legend = paste0("K = ", latentnet_res$K_opt[m]), bty = "n")
} # end m for loop
dev.off()

# ======================================================================================================================
# GENERATE 500 NEW MULTIPLEXES FROM THE POSTERIOR PREDICTIVE DISTRIBUTION
# ======================================================================================================================
rm(list = ls())

setwd("/Volumes/My Passport/PhD/PhD Project 1/LaPCoM_Code/Applications/Primary_School_Count/")

psc_multi = readRDS("primary_school_adj_matrices_count.rds")
M = dim(psc_multi)[3]
N = dim(psc_multi)[1]

load("psc_PPC_latentnet_res.Rdata")

num_PPCs = 500

latentnet_alphas = matrix(NA, nrow = num_PPCs, ncol = M)
latentnet_latent_spaces = array(NA, dim = c(N, 2, num_PPCs, M))

for (m in 1:M) {
  fit_m = latentnet_res$best_fits[[m]]
  latentnet_alphas[, m] = tail(fit_m$sample$beta, num_PPCs)
  latentnet_latent_spaces[, , , m] = fit_m$sample$Z[tail(1:dim(fit_m$sample$Z)[1], num_PPCs), , ]
} # end m for loop

set.seed(651465)

probs = array(0, dim = c(N, N, M, num_PPCs))
multiplexes = array(0, dim = c(N, N, M, num_PPCs))

prog_bar = progress::progress_bar$new(format = paste0("Generate Networks: [:bar] :percent"), total = 500, clear = F, width = 60)

for (ppc_i in 1:num_PPCs) {
  for (m in 1:M) {
    D = as.matrix(Rfast::Dist(latentnet_latent_spaces[, , ppc_i, m], square = T)) # distance matrix for Z_g
    eta = latentnet_alphas[ppc_i, m] - D
    probs[, , m, ppc_i][upper.tri(probs[, , m, ppc_i])] = (exp(eta[upper.tri(eta)])) / (1 + exp(eta[upper.tri(eta)]))
    probs[, , m, ppc_i] = probs[, , m, ppc_i] + t(probs[, , m, ppc_i])
    diag(probs[, , m, ppc_i]) = 0
    multiplexes[, , m, ppc_i][upper.tri(multiplexes[, , m, ppc_i])] = rpois((N * (N - 1)) / 2, exp(eta))
    # multiplexes[, , m, ppc_i][upper.tri(multiplexes[, , m, ppc_i])] = rbinom(((N * (N - 1)) / 2), 1, probs[, , m, ppc_i][upper.tri(probs[, , m, ppc_i])])
    multiplexes[, , m, ppc_i] = multiplexes[, , m, ppc_i] + t(multiplexes[, , m, ppc_i])
  } # end m for loop
  prog_bar$tick()
  rm(m)
} # end i for loop 
rm(ppc_i)

latentnet_PPD_data = list(multiplexes = multiplexes, probs = probs)
save(latentnet_PPD_data, file = "psc_latentnet_PPD_data.Rdata")

beepr::beep(sound = 10)

# ======================================================================================================================
# RUN PPCS
# ======================================================================================================================
rm(list = ls())

setwd("/Volumes/My Passport/PhD/PhD Project 1/LaPCoM_Code/Applications/Primary_School_Count/")

psc_multi = readRDS("primary_school_adj_matrices_count.rds")
M = dim(psc_multi)[3]
N = dim(psc_multi)[1]

load("psc_PPC_latentnet_res.Rdata")
load("psc_latentnet_PPD_data.Rdata")

num_PPCs = 500

source("LaPCoM_Initialisation_Functions.R")

multiplexes = latentnet_PPD_data$multiplexes

log_freq_counts = vector("list", length = M)
net_dists = matrix(NA, nrow = M, ncol = num_PPCs)
mean_abs_diff_counts = matrix(NA, nrow = M, ncol = num_PPCs)

prog_bar = progress::progress_bar$new(format = paste0("500 PPCs for each network: [:bar] :percent"), total = 16, clear = F, width = 60)

for (m in 1:M) {
  
  # ----------------------------------------------------------------------------------------------------------------
  # FREQUENCY OF COUNTS
  # ----------------------------------------------------------------------------------------------------------------
  count_tab = t(log(apply(multiplexes[, , m, ] + 1, 3, function(x) { tabulate(x, max(psc_multi) + 1) })))
  colnames(count_tab) = 0:max(psc_multi)
  log_freq_counts[[m]] = count_tab
  
  # ----------------------------------------------------------------------------------------------------------------
  # MEAN ABSOLUTE DIFFERENCE IN COUNTS
  # ----------------------------------------------------------------------------------------------------------------
  mean_abs_diff_counts[m, ] = apply(multiplexes[, , m, ], 3, function(x) { mean(abs(x - psc_multi[, , m])) })
  
  # ----------------------------------------------------------------------------------------------------------------
  # NETWORK-DISTANCE
  # ----------------------------------------------------------------------------------------------------------------
  gra1 = igraph::graph_from_adjacency_matrix(psc_multi[, , m], mode = "undirected")
  el1 = igraph::as_edgelist(gra1)
  
  net_dists[m, ] = sapply(seq_len(num_PPCs), function(ppc_i) {
    tryCatch({
      gra2 = igraph::graph_from_adjacency_matrix(multiplexes[, , m, ppc_i], mode = "undirected")
      el2 = igraph::as_edgelist(gra2)
      nature_dist(el1, el2, gra1, gra2, w1 = 0.45, w2 = 0.45, w3 = 0.1, net_type = "count")
    }, error = function(e) {
      return(NA)
    })
  })
  
  rm(gra1, el1, count_tab)
  
  prog_bar$tick()
  
} # end m for loop

beepr::beep(sound = 10)

save.image(file = "latentnet_PPCs_PSC")

# ======================================================================================================================
# PLOTS TO COMPARE TO LaPCoM RESULTS
# ======================================================================================================================
rm(list = ls())

setwd("/Volumes/My Passport/PhD/PhD Project 1/LaPCoM_Code/Applications/Primary_School_Count/")

load("latentnet_PPCs_PSC")
psc_multi = readRDS("primary_school_adj_matrices_count.rds")
network_names = dimnames(psc_multi)[[3]]

library(tidyverse)
library(ggplot2)
library(tidyr)
library(dplyr)

# ----------------------------------------------------------------------------------------------------------------------
# MEAN ABSOLUTE DIFFERENCE IN COUNTS
# ----------------------------------------------------------------------------------------------------------------------
df_plot = as.data.frame(mean_abs_diff_counts) # no NAs
df_plot$Network = factor(network_names, levels = network_names)  # network_names should have length 16
df_long = df_plot %>%
  pivot_longer(cols = -Network, names_to = "Replication", values_to = "MAD")

ggplot(df_long, aes(x = 0, y = MAD)) +
  geom_boxplot(fill = "white", color = "black", outlier.shape = 16) +  # Boxplot for each network
  facet_wrap(~ Network, ncol = 4) +  # 4x4 grid
  labs(x = NULL, y = "Mean Absolute Difference") +  # Labels
  theme_minimal(base_size = 14) +  # Minimal theme
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.8),
    axis.text.x = element_blank()
  )

# ----------------------------------------------------------------------------------------------------------------------
# NETWORK DISTANCES
# ----------------------------------------------------------------------------------------------------------------------
df_plot = as.data.frame(net_dists)
df_plot$Network = factor(network_names, levels = network_names)  # network_names should have length 16
df_long = df_plot %>%
  pivot_longer(cols = -Network, names_to = "Replication", values_to = "Net_Dist")

ggplot(df_long, aes(x = 0, y = Net_Dist)) +
  geom_boxplot(fill = "white", color = "black", outlier.shape = 16) +  # Boxplot for each network
  facet_wrap(~ Network, ncol = 4) +  # 4x4 grid
  labs(x = NULL, y = "Network Distance") +  # Labels
  theme_minimal(base_size = 14) +  # Minimal theme
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.8),
    axis.text.x = element_blank()
  )

# ----------------------------------------------------------------------------------------------------------------------
# LOG FREQUENCY OF COUNTS
# ----------------------------------------------------------------------------------------------------------------------
show_dim = sum(rowSums(sapply(log_freq_counts, colSums, na.rm = T)) > 0)
show_dim = 8

true_count_freqs_list = lapply(1:M, function(m) {
  true_count_freqs = log(tabulate(psc_multi[, , m] + 1, max(psc_multi) + 1))[1:show_dim]
  true_count_freqs[is.infinite(true_count_freqs)] = NA
  return(true_count_freqs)
})

df_all = data.frame()
df_truth = data.frame()

for (m in 1:M) {
  
  data_mat = log_freq_counts[[m]][, 1:show_dim]
  data_mat[data_mat == -Inf] = NA  # Handle -Inf
  
  df_temp = as.data.frame(data_mat)
  df_long = df_temp %>%
    pivot_longer(cols = everything(), names_to = "Count", values_to = "Frequency") %>%
    mutate(
      Count = as.integer(gsub("V", "", Count)),  # Now Count = 1 to 8
      Network = network_names[m]
    )
  
  df_all = bind_rows(df_all, df_long)
  
  df_truth = bind_rows(df_truth, data.frame(
    Count = 0:(show_dim - 1),
    Frequency = true_count_freqs_list[[m]],
    Network = network_names[m]
  ))
}

df_all$Network = factor(df_all$Network, levels = network_names)
df_truth$Network = factor(df_truth$Network, levels = network_names)

ggplot(df_all, aes(x = factor(Count), y = Frequency)) +
  geom_boxplot(outlier.shape = NA, fill = "grey80", color = "black") +
  geom_point(data = df_truth, aes(x = factor(Count), y = Frequency, color = "True Value"), size = 2) +
  facet_wrap(~ Network, ncol = 4) +
  scale_y_continuous(limits = c(0, 12), breaks = seq(0, 12, 2)) +
  scale_x_discrete(labels = 1:show_dim) +
  scale_color_manual(name = "", values = c("True Value" = "red")) +
  labs(x = "Count", y = "Log Frequency") +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.8)
  )
