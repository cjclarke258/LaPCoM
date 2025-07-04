# ======================================================================================================================
# DATA
# ======================================================================================================================
rm(list = ls())

setwd("/Volumes/My Passport/PhD/PhD Project 1/LaPCoM_Code/Applications/Aarhus/")

load("LaPCoM_Aarhus_PPCs.Rdata")
Aarhus_PPC_LaPCoM = list(aucs_nets = aucs_nets,
                         aucs_probs = aucs_probs,
                         square_diff_densities = square_diff_densities,
                         net_dists = net_dists,
                         ham_dists = ham_dists,
                         F1s = F1s)

load("Aarhus_latentnet_PPC_res.Rdata")
Aarhus_PPC_latentnet = list(aucs_nets = aucs_nets,
                         aucs_probs = aucs_probs,
                         square_diff_densities = square_diff_densities,
                         net_dists = net_dists,
                         ham_dists = ham_dists,
                         F1s = F1s)

PPC_info = list(Aarhus_PPC_LaPCoM = Aarhus_PPC_LaPCoM, Aarhus_PPC_latentnet = Aarhus_PPC_latentnet)
save(PPC_info, file = "Aarhus_PPC_LaPCoM_latentnet_Comparison.Rdata")

# ======================================================================================================================
# PLOTS
# ======================================================================================================================
rm(list = ls())
load("Aarhus_PPC_LaPCoM_latentnet_Comparison.Rdata")
aarhus_multi = readRDS("aarhus_adjacency_matrices_updated.rds")
M = dim(aarhus_multi)[3]
N = dim(aarhus_multi)[1]
cols = khroma::color("muted")(9)

# pdf("Aarhus_PPC_Plots/Updated_Aarhus_Combined_PPC_Density.pdf", height = 2.75, width = 10)
# par(mfrow = c(1, 5))
# par(mar = c(4, 4, 2, 0.5) + 0.4) # bltr
# for (m in 1:M) { 
#   temp_df = data.frame(latentnet = as.vector(PPC_info$Aarhus_PPC_latentnet$square_diff_densities[m, ]),
#                        LaPCoM = as.vector(PPC_info$Aarhus_PPC_LaPCoM$square_diff_densities[m, ]))
#   boxplot(temp_df, xlab = "", ylab = "", ylim = c(0, 0.35), col = cols, main = dimnames(aarhus_multi)[[3]][m], las = 1, cex.axis = 1, cex.main = 1.5) 
# } # end m for loop
# dev.off()
# 
# pdf("Aarhus_PPC_Plots/Updated_Aarhus_Combined_PPC_AUC_probs.pdf", height = 2.75, width = 10)
# par(mfrow = c(1, 5))
# par(mar = c(4, 4, 2, 0.5) + 0.4) # bltr
# for (m in 1:M) { 
#   temp_df = data.frame(latentnet = as.vector(PPC_info$Aarhus_PPC_latentnet$aucs_probs[m, ]),
#                        LaPCoM = as.vector(PPC_info$Aarhus_PPC_LaPCoM$aucs_probs[m, ]))
#   boxplot(temp_df, xlab = "", ylab = "", ylim = c(0, 1), col = cols,main = dimnames(aarhus_multi)[[3]][m], las = 1, cex.axis = 1, cex.main = 1.5) 
#   abline(h = network::network.density(network::as.network(aarhus_multi[, , m])), col = "red")
# } # end m for loop
# dev.off()
# 
# pdf("Aarhus_PPC_Plots/Updated_Aarhus_Combined_PPC_Net_Distances.pdf", height = 2.75, width = 10)
# par(mfrow = c(1, 5))
# par(mar = c(4, 4, 2, 0.5) + 0.4) # bltr
# for (m in 1:M) { 
#   temp_df = data.frame(latentnet = as.vector(PPC_info$Aarhus_PPC_latentnet$net_dists[m, ]),
#                        LaPCoM = as.vector(PPC_info$Aarhus_PPC_LaPCoM$net_dists[m, ]))
#   boxplot(temp_df, xlab = "", ylab = "", ylim = c(0, 0.6), col = cols, main = dimnames(aarhus_multi)[[3]][m], las = 1, cex.axis = 1, cex.main = 1.5) 
# } # end m for loop
# dev.off()
# 
# pdf("Aarhus_PPC_Plots/Updated_Aarhus_Combined_PPC_Ham_Distances.pdf", height = 2.75, width = 10)
# par(mfrow = c(1, 5))
# par(mar = c(4, 4, 2, 0.5) + 0.4) # bltr
# for (m in 1:M) { 
#   temp_df = data.frame(latentnet = as.vector(PPC_info$Aarhus_PPC_latentnet$ham_dists[m, ]),
#                        LaPCoM = as.vector(PPC_info$Aarhus_PPC_LaPCoM$ham_dists[m, ]))
#   boxplot(temp_df, xlab = "", ylab = "", ylim = c(0, 0.04), col = cols, main = dimnames(aarhus_multi)[[3]][m], las = 1, cex.axis = 1, cex.main = 1.5) 
# } # end m for loop
# dev.off()
# 
# pdf("Aarhus_PPC_Plots/Updated_Aarhus_Combined_PPC_F1s.pdf", height = 2.75, width = 10)
# par(mfrow = c(1, 5))
# par(mar = c(4, 4, 2, 0.5) + 0.4) # bltr
# for (m in 1:M) { 
#   temp_df = data.frame(latentnet = as.vector(PPC_info$Aarhus_PPC_latentnet$F1s[m, ]),
#                        LaPCoM = as.vector(PPC_info$Aarhus_PPC_LaPCoM$F1s[m, ]))
#   boxplot(temp_df, xlab = "", ylab = "", ylim = c(0, 1), col = cols, main = dimnames(aarhus_multi)[[3]][m], las = 1, cex.axis = 1, cex.main = 1.5) 
# } # end m for loop
# dev.off()

# 

mean(apply(PPC_info$Aarhus_PPC_LaPCoM$F1s[1:3, ], 1, median))
mean(apply(PPC_info$Aarhus_PPC_latentnet$F1s[1:3, ], 1, median))

mean(apply(PPC_info$Aarhus_PPC_LaPCoM$F1s[1:3, ], 1, IQR))
mean(apply(PPC_info$Aarhus_PPC_latentnet$F1s[1:3, ], 1, IQR))

apply(PPC_info$Aarhus_PPC_LaPCoM$F1s[4:5, ], 1, median)
apply(PPC_info$Aarhus_PPC_latentnet$F1s[4:5, ], 1, median)

apply(PPC_info$Aarhus_PPC_LaPCoM$F1s[4:5, ], 1, IQR)
apply(PPC_info$Aarhus_PPC_latentnet$F1s[4:5, ], 1, IQR)

# 

mean(apply(PPC_info$Aarhus_PPC_LaPCoM$square_diff_densities[c(1:3, 5), ], 1, median))
mean(apply(PPC_info$Aarhus_PPC_latentnet$square_diff_densities[c(1:3, 5), ], 1, median))

mean(apply(PPC_info$Aarhus_PPC_LaPCoM$square_diff_densities[c(1:3, 5), ], 1, IQR))
mean(apply(PPC_info$Aarhus_PPC_latentnet$square_diff_densities[c(1:3, 5), ], 1, IQR))

median(PPC_info$Aarhus_PPC_LaPCoM$square_diff_densities[4, ])
median(PPC_info$Aarhus_PPC_latentnet$square_diff_densities[4, ])

IQR(PPC_info$Aarhus_PPC_LaPCoM$square_diff_densities[4, ])
IQR(PPC_info$Aarhus_PPC_latentnet$square_diff_densities[4, ])

# 

mean(apply(PPC_info$Aarhus_PPC_LaPCoM$net_dists, 1, median))
mean(apply(PPC_info$Aarhus_PPC_latentnet$net_dists, 1, median), na.rm = T)

mean(apply(PPC_info$Aarhus_PPC_LaPCoM$net_dists, 1, IQR))
mean(apply(PPC_info$Aarhus_PPC_latentnet$net_dists, 1, IQR, na.rm = T))

# 

mean(apply(PPC_info$Aarhus_PPC_LaPCoM$ham_dists, 1, median))
mean(apply(PPC_info$Aarhus_PPC_latentnet$ham_dists, 1, median))

mean(apply(PPC_info$Aarhus_PPC_LaPCoM$ham_dists, 1, IQR))
mean(apply(PPC_info$Aarhus_PPC_latentnet$ham_dists, 1, IQR))

# 

library(ggplot2)
library(dplyr)
library(tidyr)

model_cols = unname(cols)

plot_ppc_metric = function(ppc_latentnet, ppc_LaPCoM, ylab, title, ylim_max) {
  
  latentnet_df = as.data.frame(ppc_latentnet)
  LaPCoM_df   = as.data.frame(ppc_LaPCoM)
  
  latentnet_df$Network = LETTERS[1:5]
  latentnet_df$Model = "latentnet"
  
  LaPCoM_df$Network = LETTERS[1:5]
  LaPCoM_df$Model = "LaPCoM"
  
  plot_df = bind_rows(latentnet_df, LaPCoM_df) %>%
    pivot_longer(cols = -c(Network, Model), names_to = "Rep", values_to = "Value")
  plot_df$Model = factor(plot_df$Model, levels = c("latentnet", "LaPCoM"))
  
  ggplot(plot_df, aes(x = Model, y = Value, fill = Model)) +
    geom_boxplot(width = 0.6, outlier.size = 0.8, color = "black") +
    facet_wrap(~ Network, nrow = 1) +
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

pdf("Aarhus_PPC_Plots/Aarhus_Combined_PPC_Net_Distances.pdf", width = 15, height = 5)
plot_ppc_metric(ppc_latentnet = PPC_info$Aarhus_PPC_latentnet$net_dists, 
                ppc_LaPCoM = PPC_info$Aarhus_PPC_LaPCoM$net_dists, 
                ylab = "Distance from Observed Network", 
                title = "Posterior Predictive Check: Network Distances",
                ylim_max = 0.6)
dev.off()

pdf("Aarhus_PPC_Plots/Aarhus_Combined_PPC_Ham_Distances.pdf", width = 15, height = 5)
plot_ppc_metric(ppc_latentnet = PPC_info$Aarhus_PPC_latentnet$ham_dists, 
                ppc_LaPCoM = PPC_info$Aarhus_PPC_LaPCoM$ham_dists, 
                ylab = "Hamming Distance", 
                title = "Posterior Predictive Check: Hamming Distance",
                ylim_max = 0.035)
dev.off()

pdf("Aarhus_PPC_Plots/Aarhus_Combined_PPC_F1s.pdf", width = 15, height = 5)
plot_ppc_metric(ppc_latentnet = PPC_info$Aarhus_PPC_latentnet$F1s, 
                ppc_LaPCoM = PPC_info$Aarhus_PPC_LaPCoM$F1s, 
                ylab = expression(F[1]), 
                title = paste0("Posterior Predictive Check: ", expression(F[1])),
                ylim_max = 1)
dev.off()

pdf("Aarhus_PPC_Plots/Aarhus_Combined_PPC_Density.pdf", width = 15, height = 5)
plot_ppc_metric(ppc_latentnet = PPC_info$Aarhus_PPC_latentnet$square_diff_densities, 
                ppc_LaPCoM = PPC_info$Aarhus_PPC_LaPCoM$square_diff_densities, 
                ylab = "Squared Difference in Density", 
                title = "Posterior Predictive Check: Squared Difference in Density",
                ylim_max = 0.5)
dev.off()

# 

pdf("Aarhus_PPC_Plots/Aarhus_Combined_PPC_AUC_probs.pdf", width = 15, height = 5)

ppc_latentnet = PPC_info$Aarhus_PPC_latentnet$aucs_probs
ppc_LaPCoM = PPC_info$Aarhus_PPC_LaPCoM$aucs_probs
ylab = "Area Under the Curve (AUC)"
title = "Posterior Predictive Check: AUC"

latentnet_df = as.data.frame(ppc_latentnet)
LaPCoM_df   = as.data.frame(ppc_LaPCoM)

latentnet_df$Network = LETTERS[1:5]
latentnet_df$Model = "latentnet"

LaPCoM_df$Network = LETTERS[1:5]
LaPCoM_df$Model = "LaPCoM"

plot_df = bind_rows(latentnet_df, LaPCoM_df) %>%
  pivot_longer(cols = -c(Network, Model), names_to = "Rep", values_to = "Value")
plot_df$Model = factor(plot_df$Model, levels = c("latentnet", "LaPCoM"))

true_df = data.frame(
  Network = LaPCoM_df$Network,
  TrueDensity = apply(aarhus_multi, 3, function(x) { network::network.density(network::as.network(x)) })
)

ggplot(plot_df, aes(x = Model, y = Value, fill = Model)) +
  geom_boxplot(width = 0.6, outlier.size = 0.8, color = "black") +
  facet_wrap(~ Network, nrow = 1) +
  geom_hline(data = true_df, aes(yintercept = TrueDensity), color = "red", linewidth = 1) +
  scale_fill_manual(values = model_cols) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = NULL, y = ylab) +
  theme_minimal(base_size = 16) +
  theme(
    strip.text = element_text(size = 19, face = "bold"),
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
