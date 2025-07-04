# ======================================================================================================================
# DATA
# ======================================================================================================================
rm(list = ls())

setwd("/Volumes/My Passport/PhD/PhD Project 1/LaPCoM_Code/Applications/Krackhardt/")

load("LaPCoM_Krack_PPCs.Rdata")
krack_PPC_LaPCoM = list(aucs_nets = aucs_nets,
                         aucs_probs = aucs_probs,
                         square_diff_densities = square_diff_densities,
                         net_dists = net_dists,
                         ham_dists = ham_dists,
                         F1s = F1s)

load("latentnet_PPC_res.Rdata")
krack_PPC_latentnet = list(aucs_nets = aucs_nets,
                         aucs_probs = aucs_probs,
                         square_diff_densities = square_diff_densities,
                         net_dists = net_dists,
                         ham_dists = ham_dists,
                         F1s = F1s)

PPC_info = list(krack_PPC_LaPCoM = krack_PPC_LaPCoM, krack_PPC_latentnet = krack_PPC_latentnet)
save(PPC_info, file = "Krackhardt_PPC_LaPCoM_latentnet_Comparison.Rdata")

# ======================================================================================================================
# PLOTS
# ======================================================================================================================
rm(list = ls())
load("Krackhardt_PPC_LaPCoM_latentnet_Comparison.Rdata")
load("Krackhardt_Multiplex_Directed.Rdata")

M = N = 21
cols = khroma::color("muted")(9)

# png("Krack_PPC_Plots/Krack_Combined_PPC_Density.png", height = 420 * 2, width = 480 * 3)
# par(mfrow = c(3, 7))
# par(mar = c(7, 5, 6, 1) + 0.1) # bltr
# for (m in 1:M) { 
#   temp_df = data.frame(latentnet = as.vector(PPC_info$krack_PPC_latentnet$square_diff_densities[m, ]),
#                        LaPCoM = as.vector(PPC_info$krack_PPC_LaPCoM$square_diff_densities[m, ]))
#   boxplot(temp_df, xlab = "", ylab = "", ylim = c(0, 0.35), col = cols, main = LETTERS[m], las = 1, cex.axis = 1.5, cex.main = 2) 
#   # text(labels = paste("Median:\n", round(median(temp_df$LaPCoM), 3)), x = 1, y = 0.375, cex = 1.75, col = cols[1])
#   # text(labels = paste("Median:\n", round(median(temp_df$latentnet), 3)), x = 2, y = 0.375, cex = 1.75, col = cols[2])
# } # end m for loop
# dev.off()
# 
# png("Krack_PPC_Plots/Krack_Combined_PPC_AUC_probs.png", height = 420 * 2, width = 480 * 3)
# par(mfrow = c(3, 7))
# par(mar = c(7, 5, 6, 1) + 0.1) # bltr
# for (m in 1:M) { 
#   temp_df = data.frame(latentnet = as.vector(PPC_info$krack_PPC_latentnet$aucs_probs[m, ]),
#                        LaPCoM = as.vector(PPC_info$krack_PPC_LaPCoM$aucs_probs[m, ]))
#   boxplot(temp_df, xlab = "", ylab = "", ylim = c(0, 1), col = cols,main = LETTERS[m], las = 1, cex.axis = 1.5, cex.main = 2) 
#   abline(h = network::network.density(network::as.network(krack_multi[, , m])), col = "red")
#   # text(labels = paste(round(median(temp_df$LaPCoM), 3)), x = 1, y = max(temp_df$LaPCoM) + 0.1, cex = 1.5, col = cols[1])
#   # text(labels = paste(round(median(temp_df$latentnet), 3)), x = 2, y = max(temp_df$latentnet) + 0.1, cex = 1.5, col = cols[2])
#   } # end m for loop
# dev.off()
# 
# png("Krack_PPC_Plots/Krack_Combined_PPC_Net_Distances.png", height = 420 * 2, width = 480 * 3)
# par(mfrow = c(3, 7))
# par(mar = c(7, 5, 6, 1) + 0.1) # bltr
# for (m in 1:M) { 
#   temp_df = data.frame(latentnet = as.vector(PPC_info$krack_PPC_latentnet$net_dists[m, ]),
#                        LaPCoM = as.vector(PPC_info$krack_PPC_LaPCoM$net_dists[m, ]))
#   boxplot(temp_df, xlab = "", ylab = "", ylim = c(0, 0.6), col = cols, main = LETTERS[m], las = 1, cex.axis = 1.5, cex.main = 2) 
#   # text(labels = paste("Median:\n", round(median(temp_df$LaPCoM), 3)), x = 1, y = 0.375, cex = 1.75, col = cols[1])
#   # text(labels = paste("Median:\n", round(median(temp_df$latentnet), 3)), x = 2, y = 0.375, cex = 1.75, col = cols[2])
# } # end m for loop
# dev.off()
# 
# png("Krack_PPC_Plots/Krack_Combined_PPC_Ham_Distances.png", height = 420 * 2, width = 480 * 3)
# par(mfrow = c(3, 7))
# par(mar = c(7, 5, 6, 1) + 0.1) # bltr
# for (m in 1:M) { 
#   temp_df = data.frame(latentnet = as.vector(PPC_info$krack_PPC_latentnet$ham_dists[m, ]),
#                        LaPCoM = as.vector(PPC_info$krack_PPC_LaPCoM$ham_dists[m, ]))
#   boxplot(temp_df, xlab = "", ylab = "", ylim = c(0, 0.2), col = cols, main = LETTERS[m], las = 1, cex.axis = 1.5, cex.main = 2) 
#   # text(labels = paste("Median:\n", round(median(temp_df$LaPCoM), 3)), x = 1, y = 0.375, cex = 1.75, col = cols[1])
#   # text(labels = paste("Median:\n", round(median(temp_df$latentnet), 3)), x = 2, y = 0.375, cex = 1.75, col = cols[2])
# } # end m for loop
# dev.off()
# 
# png("Krack_PPC_Plots/Krack_Combined_PPC_F1s.png", height = 420 * 2, width = 480 * 3)
# par(mfrow = c(3, 7))
# par(mar = c(7, 5, 6, 1) + 0.1) # bltr
# for (m in 1:M) { 
#   temp_df = data.frame(latentnet = as.vector(PPC_info$krack_PPC_latentnet$F1s[m, ]),
#                        LaPCoM = as.vector(PPC_info$krack_PPC_LaPCoM$F1s[m, ]))
#   boxplot(temp_df, xlab = "", ylab = "", ylim = c(0, 1), col = cols, main = LETTERS[m], las = 1, cex.axis = 1.5, cex.main = 2) 
#   # text(labels = paste("Median:\n", round(median(temp_df$LaPCoM), 3)), x = 1, y = 0.375, cex = 1.75, col = cols[1])
#   # text(labels = paste("Median:\n", round(median(temp_df$latentnet), 3)), x = 2, y = 0.375, cex = 1.75, col = cols[2])
# } # end m for loop
# dev.off()

# 

min(apply(PPC_info$krack_PPC_LaPCoM$F1s, 1, median))
which.min(apply(PPC_info$krack_PPC_LaPCoM$F1s, 1, median))

mean(apply(PPC_info$krack_PPC_LaPCoM$F1s[-which.min(apply(PPC_info$krack_PPC_LaPCoM$F1s, 1, median)), ], 1, median))

# 

max(apply(PPC_info$krack_PPC_LaPCoM$square_diff_densities, 1, median))
which.max(apply(PPC_info$krack_PPC_LaPCoM$square_diff_densities, 1, median))

mean(apply(PPC_info$krack_PPC_LaPCoM$square_diff_densities[-which.max(apply(PPC_info$krack_PPC_LaPCoM$square_diff_densities, 1, median)), ], 1, median))

# 

mean(apply(PPC_info$krack_PPC_LaPCoM$net_dists, 1, median))

mean(apply(PPC_info$krack_PPC_latentnet$net_dists, 1, IQR, na.rm = T))
mean(apply(PPC_info$krack_PPC_LaPCoM$net_dists, 1, IQR))

# 

mean(apply(PPC_info$krack_PPC_LaPCoM$ham_dists, 1, median))

# 

library(ggplot2)
library(dplyr)
library(tidyr)

model_cols = unname(cols)

plot_ppc_metric = function(ppc_latentnet, ppc_LaPCoM, ylab, title, ylim_max) {
  
  latentnet_df = as.data.frame(ppc_latentnet)
  LaPCoM_df   = as.data.frame(ppc_LaPCoM)
  
  latentnet_df$Network = LETTERS[1:21]
  latentnet_df$Model = "latentnet"
  
  LaPCoM_df$Network = LETTERS[1:21]
  LaPCoM_df$Model = "LaPCoM"
  
  plot_df = bind_rows(latentnet_df, LaPCoM_df) %>%
    pivot_longer(cols = -c(Network, Model), names_to = "Rep", values_to = "Value")
  plot_df$Model = factor(plot_df$Model, levels = c("latentnet", "LaPCoM"))
  
  ggplot(plot_df, aes(x = Model, y = Value, fill = Model)) +
    geom_boxplot(width = 0.6, outlier.size = 0.8, color = "black") +
    facet_wrap(~ Network, nrow = 3) +
    scale_fill_manual(values = model_cols) +
    coord_cartesian(ylim = c(0, ylim_max)) +
    labs(x = NULL, y = ylab) +
    theme_minimal(base_size = 16) +
    theme(
      strip.text = element_text(size = 18, face = "bold"),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.8),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

} # end plot_ppc_metric function

pdf("Krack_PPC_Plots/Krack_Combined_PPC_Net_Distances.pdf", width = 16, height = 8)
plot_ppc_metric(ppc_latentnet = PPC_info$krack_PPC_latentnet$net_dists, 
                ppc_LaPCoM = PPC_info$krack_PPC_LaPCoM$net_dists, 
                ylab = "Distance from Observed Network", 
                title = "Posterior Predictive Check: Network Distances",
                ylim_max = 0.6)
dev.off()

pdf("Krack_PPC_Plots/Krack_Combined_PPC_Ham_Distances.pdf", width = 16, height = 8)
plot_ppc_metric(ppc_latentnet = PPC_info$krack_PPC_latentnet$ham_dists, 
                ppc_LaPCoM = PPC_info$krack_PPC_LaPCoM$ham_dists, 
                ylab = "Hamming Distance", 
                title = "Posterior Predictive Check: Hamming Distance",
                ylim_max = 0.15)
dev.off()

pdf("Krack_PPC_Plots/Krack_Combined_PPC_F1s.pdf", width = 16, height = 8)
plot_ppc_metric(ppc_latentnet = PPC_info$krack_PPC_latentnet$F1s, 
                ppc_LaPCoM = PPC_info$krack_PPC_LaPCoM$F1s, 
                ylab = expression(F[1]), 
                title = paste0("Posterior Predictive Check: ", expression(F[1])),
                ylim_max = 1)
dev.off()

pdf("Krack_PPC_Plots/Krack_Combined_PPC_Density.pdf", width = 16, height = 8)
plot_ppc_metric(ppc_latentnet = PPC_info$krack_PPC_latentnet$square_diff_densities, 
                ppc_LaPCoM = PPC_info$krack_PPC_LaPCoM$square_diff_densities, 
                ylab = "Squared Difference in Density", 
                title = "Posterior Predictive Check: Squared Difference in Density",
                ylim_max = 0.4)
dev.off()

# 

pdf("Krack_PPC_Plots/Krack_Combined_PPC_AUC_probs.pdf", width = 16, height = 8)

ppc_latentnet = PPC_info$krack_PPC_latentnet$aucs_probs
ppc_LaPCoM = PPC_info$krack_PPC_LaPCoM$aucs_probs
ylab = "Area Under the Curve (AUC)"
title = "Posterior Predictive Check: AUC"

latentnet_df = as.data.frame(ppc_latentnet)
LaPCoM_df   = as.data.frame(ppc_LaPCoM)

latentnet_df$Network = LETTERS[1:21]
latentnet_df$Model = "latentnet"

LaPCoM_df$Network = LETTERS[1:21]
LaPCoM_df$Model = "LaPCoM"

plot_df = bind_rows(latentnet_df, LaPCoM_df) %>%
  pivot_longer(cols = -c(Network, Model), names_to = "Rep", values_to = "Value")
plot_df$Model = factor(plot_df$Model, levels = c("latentnet", "LaPCoM"))

true_df = data.frame(
  Network = LaPCoM_df$Network,
  TrueDensity = apply(krack_multi, 3, function(x) { network::network.density(network::as.network(x)) })
)

ggplot(plot_df, aes(x = Model, y = Value, fill = Model)) +
  geom_boxplot(width = 0.6, outlier.size = 0.8, color = "black") +
  facet_wrap(~ Network, nrow = 3) +
  geom_hline(data = true_df, aes(yintercept = TrueDensity), color = "red", linewidth = 1) +
  scale_fill_manual(values = model_cols) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = NULL, y = ylab) +
  theme_minimal(base_size = 16) +
  theme(
    strip.text = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.8),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

dev.off()
