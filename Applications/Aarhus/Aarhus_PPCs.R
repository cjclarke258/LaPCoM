########################################################################################################################
# AARHUS DATA APPLICATION (POSTERIOR PREDICTIVE CHECKS)
########################################################################################################################

rm(list = ls())

setwd("/Volumes/My Passport/PhD/PhD Project 1/LaPCoM_Code/Applications/Aarhus/")
source("LaPCoM_Initialisation_Functions.R")

# ======================================================================================================================
# LOAD DATA
# ======================================================================================================================
# ACTUAL DATA
aarhus_multi = readRDS("aarhus_adjacency_matrices_updated.rds")

# LaPCoM RESULTS
# load("pp_res_aarhus.Rdata")
load("pp_res_aarhus_PPMethodC.Rdata")

# ======================================================================================================================
# EXTRACT THE LAST 100 ALPHAS, TAUS AND LATENT SPACES FOR EACH SIMULATION
# ======================================================================================================================
best_chain_num = 31 #CHANGE THIS ACCORDINGLY
num_PPCs = 500
num_samps = dim(aarhus_pp$pp_res_aarhus$chains$Z[[best_chain_num]])[4]

LaPCoM_PPC_alphas = tail(aarhus_pp$pp_res_aarhus$chains$alpha[[best_chain_num]], num_PPCs)
LaPCoM_PPC_Z = aarhus_pp$pp_res_aarhus$chains$Z[[best_chain_num]][, , , tail(1:num_samps, num_PPCs)]
LaPCoM_PPC_C = aarhus_pp$pp_res_aarhus$chains$C[[best_chain_num]][, tail(1:num_samps, num_PPCs)]

# ======================================================================================================================
# GENERATE 100 POSTERIOR PREDICTIVE MULTIPLEXES
# ======================================================================================================================
N = dim(aarhus_multi)[1] 
M = dim(aarhus_multi)[3]

multis_PPC = probs = array(0, dim = c(N, N, M, num_PPCs))

set.seed(117487)

for (ppc_i in 1:num_PPCs) {
  
  # ------------------------------------------------------------------------------------------------------------------
  # GENERATE THE NETWORKS
  # ------------------------------------------------------------------------------------------------------------------
  for (m in 1:M) {
    g = LaPCoM_PPC_C[m, ppc_i] # which cluster network m belongs to
    D = as.matrix(Rfast::Dist(LaPCoM_PPC_Z[, , g, ppc_i], square = T)) # distance matrix for Z_g
    eta = LaPCoM_PPC_alphas[ppc_i] - D
    # probs[, , m, ppc_i] = (exp(eta)) / (1 + exp(eta))
    probs[, , m, ppc_i][upper.tri(probs[, , m, ppc_i])] = (exp(eta[upper.tri(eta)])) / (1 + exp(eta[upper.tri(eta)]))
    probs[, , m, ppc_i] = probs[, , m, ppc_i] + t(probs[, , m, ppc_i])
    diag(probs[, , m, ppc_i]) = 0
    # multis_PPC[, , m, ppc_i] = matrix(rbinom(N^2, 1, probs[, , m, ppc_i]), ncol = N, nrow = N)
    multis_PPC[, , m, ppc_i][upper.tri(multis_PPC[, , m, ppc_i])] = rbinom(((N * (N - 1)) / 2), 1, probs[, , m, ppc_i][upper.tri(probs[, , m, ppc_i])])
    multis_PPC[, , m, ppc_i] = multis_PPC[, , m, ppc_i] + t(multis_PPC[, , m, ppc_i])
    rm(D, g, eta)
  } # end m for loop
  rm(m)

} # end ppc_i for loop
rm(ppc_i)

# ======================================================================================================================
# USE SUMMARY STATISTICS TO COMPARE THE SIMULATED POSTERIOR PREDICTIVE MULTIPLEXES TO THE TRUE MULTIPLEX
# ======================================================================================================================
square_diff_densities = matrix(NA, nrow = M, ncol = num_PPCs)
aucs_probs = matrix(NA, nrow = M, ncol = num_PPCs)
aucs_nets = matrix(NA, nrow = M, ncol = num_PPCs)
net_dists = matrix(NA, nrow = M, ncol = num_PPCs)
ham_dists = matrix(NA, nrow = M, ncol = num_PPCs)
F1s = matrix(NA, nrow = M, ncol = num_PPCs)

prog_bar = progress::progress_bar$new(format = paste0("500 PPCs: [:bar] :percent"), total = 500, clear = F, width = 60)

for (ppc_i in 1:num_PPCs) {
  
  for (m in 1:M) {
    
    # ----------------------------------------------------------------------------------------------------------------
    # DENSITY
    # ----------------------------------------------------------------------------------------------------------------
    true_dens = network::network.density(network::as.network(aarhus_multi[, , m]))
    synth_dens = network::network.density(network::as.network(multis_PPC[, , m, ppc_i]))
    square_diff_densities[m, ppc_i] = (true_dens - synth_dens)^2
    
    # ----------------------------------------------------------------------------------------------------------------
    # PRECISION-RECALL AUC WITH PROBABILITIES
    # ----------------------------------------------------------------------------------------------------------------
    pred_obj = ROCR::prediction(as.vector(probs[, , m, ppc_i]), as.vector(aarhus_multi[, , m]))
    pr_obj = ROCR::performance(pred_obj, "prec", "rec")
    # plot(pr_obj)
    aucs_probs[m, ppc_i] = ROCR::performance(pred_obj, "aucpr")@y.values[[1]][1]
    rm(pred_obj, pr_obj)
    
    # ----------------------------------------------------------------------------------------------------------------
    # PRECISION-RECALL AUC WITH NETWORKS
    # ----------------------------------------------------------------------------------------------------------------
    pred_obj_nets = ROCR::prediction(as.vector(multis_PPC[, , m, ppc_i]), as.vector(aarhus_multi[, , m]))
    pr_obj_nets = ROCR::performance(pred_obj_nets, "prec", "rec")
    # plot(pr_obj)
    aucs_nets[m, ppc_i] = ROCR::performance(pred_obj_nets, "aucpr")@y.values[[1]][1]
    
    # ----------------------------------------------------------------------------------------------------------------
    # NETWORK-DISTANCE
    # ----------------------------------------------------------------------------------------------------------------
    gra1 = igraph::graph_from_adjacency_matrix(aarhus_multi[, , m], "undirected")
    gra2 = igraph::graph_from_adjacency_matrix(multis_PPC[, , m, ppc_i], "undirected")
    el1 = cbind(igraph::as_edgelist(gra1))
    el2 = cbind(igraph::as_edgelist(gra2))
    net_dists[m, ppc_i] = nature_dist(el1, el2, gra1, gra2, w1 = 0.45, w2 = 0.45, w3 = 0.1, net_type = "binary")
    
    # ----------------------------------------------------------------------------------------------------------------
    # HAMMING-DISTANCE
    # ----------------------------------------------------------------------------------------------------------------
    net_array = array(NA, dim = c(N, N, 2))
    net_array[, , 1] = aarhus_multi[, , m]
    net_array[, , 2] = multis_PPC[, , m, ppc_i]
    
    ham_dists[m, ppc_i] = sna::hdist(net_array, g1 = 1, g2 = 2, normalize = T)
    
    # ----------------------------------------------------------------------------------------------------------------
    # F1-SCORE
    # ----------------------------------------------------------------------------------------------------------------
    F1s[m, ppc_i] = MLmetrics::F1_Score(aarhus_multi[, , m], multis_PPC[, , m, ppc_i])
    
  } # end m for loop
  rm(m)
  
  prog_bar$tick()
  
} # end ppc_i for loop
rm(ppc_i)

rm(el1, el2, gra1, gra2, pr_obj_nets, pred_obj_nets, prog_bar, true_dens, synth_dens, alpha_fun, entropy, nature_dist, nnd, node_distance)

save.image(file = "LaPCoM_Aarhus_PPCs.Rdata")
beepr::beep(sound = 10)

# ======================================================================================================================
# PLOTS
# ======================================================================================================================
rm(list = ls())

load("LaPCoM_Aarhus_PPCs.Rdata")

pdf("Aarhus_PPC_Plots/Aarhus_PPC_Density.pdf", height = 2.75, width = 10)
par(mfrow = c(1, 5))
par(mar = c(4, 4, 2, 0.5) + 0.4) # bltr
for (m in 1:M) { 
  boxplot(as.vector(square_diff_densities[m, ]), 
          xlab = c(paste("Median =", round(median(square_diff_densities[m, ]), 3))), ylab = "", ylim = c(0, 0.3),
          main = dimnames(aarhus_multi)[[3]][m], las = 1, cex.axis = 1, cex.lab = 1.2, cex.main = 1.5) 
  # legend("topright", legend = paste("Median =", round(median(square_diff_densities[m, ]), 3)), bty = "n", cex = 1)
} # end m for loop
dev.off()

pdf("Aarhus_PPC_Plots/Aarhus_PPC_AUC_probs.pdf", height = 2.75, width = 10)
par(mfrow = c(1, 5))
par(mar = c(4, 4, 2, 0.5) + 0.4) # bltr
for (m in 1:M) {
  boxplot(as.vector(aucs_probs[m, ]), 
          xlab = paste("Median =", round(median(aucs_probs[m, ]), 2)), ylab = "", ylim = c(0, 1),
          main = dimnames(aarhus_multi)[[3]][m], las = 1, cex.axis = 1, cex.lab = 1.2, cex.main = 1.5) 
  # legend("topright", legend = paste("Median =", round(median(aucs_probs[m, ]), 2)), bty = "n", cex = 1)
  abline(h = network::network.density(network::as.network(aarhus_multi[, , m])), col = "red")
} # end m for loop
dev.off()

pdf("Aarhus_PPC_Plots/Aarhus_PPC_AUC_nets.pdf", height = 2.75, width = 10)
par(mfrow = c(1, 5))
par(mar = c(4, 4, 2, 0.5) + 0.4) # bltr
for (m in 1:M) {
  boxplot(as.vector(aucs_nets[m, ]), 
          xlab = paste("Median =", round(median(aucs_nets[m, ]), 2)), ylab = "", ylim = c(0, 1),
          main = dimnames(aarhus_multi)[[3]][m], las = 1, cex.axis = 1, cex.lab = 1.2, cex.main = 1.5) 
  # legend("topright", legend = paste("Median =", round(median(aucs_nets[m, ]), 2)), bty = "n", cex = 1)
  abline(h = network::network.density(network::as.network(aarhus_multi[, , m])), col = "red")
} # end m for loop
dev.off()

pdf("Aarhus_PPC_Plots/Aarhus_PPC_Net_Distances.pdf", height = 2.75, width = 10)
par(mfrow = c(1, 5))
par(mar = c(4, 4, 2, 0.5) + 0.4) # bltr
for (m in 1:M) {
  boxplot(as.vector(net_dists[m, ]), 
          xlab = paste("Median =", round(median(net_dists[m, ]), 2)), ylab = "", ylim = c(0, 0.5),
          main = dimnames(aarhus_multi)[[3]][m], las = 1, cex.axis = 1, cex.lab = 1.2, cex.main = 1.5) 
  # legend("bottomright", legend = paste("Median =", round(median(net_dists[m, ]), 2)), bty = "n", cex = 1)
} # end m for loop
dev.off()

pdf("Aarhus_PPC_Plots/Aarhus_PPC_Ham_Distances.pdf", height = 2.75, width = 10)
par(mfrow = c(1, 5))
par(mar = c(4, 4, 2, 0.5) + 0.4) # bltr
for (m in 1:M) {
  boxplot(as.vector(ham_dists[m, ]), 
          xlab = paste("Median =", round(median(ham_dists[m, ]), 2)), ylab = "", ylim = c(0, 0.2),
          main = dimnames(aarhus_multi)[[3]][m], las = 1, cex.axis = 1, cex.lab = 1.2, cex.main = 1.5) 
  # legend("topright", legend = paste("Median =", round(median(ham_dists[m, ]), 2)), bty = "n", cex = 1)
} # end m for loop
dev.off()

pdf("Aarhus_PPC_Plots/Aarhus_PPC_F1s.pdf", height = 2.75, width = 10)
par(mfrow = c(1, 5))
par(mar = c(4, 4, 2, 0.5) + 0.4) # bltr
for (m in 1:M) {
  boxplot(as.vector(F1s[m, ]), 
          xlab = paste("Median =", round(median(F1s[m, ]), 2)), ylab = "", ylim = c(0, 1),
          main = dimnames(aarhus_multi)[[3]][m], las = 1, cex.axis = 1, cex.lab = 1.2, cex.main = 1.5) 
  # legend("bottomright", legend = paste("Median =", round(median(F1s[m, ]), 2)), bty = "n", cex = 1)
} # end m for loop
dev.off()
