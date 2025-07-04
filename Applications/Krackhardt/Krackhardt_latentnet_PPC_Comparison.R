########################################################################################################################
# KRACKHARDT DATA APPLICATION (POSTERIOR PREDICTIVE CHECKS) [USING LATENTNET AS A COMPARISON]
########################################################################################################################

rm(list = ls())

setwd("/Volumes/My Passport/PhD/PhD Project 1/LaPCoM_Code/Applications/Krackhardt/")
source("LaPCoM_Initialisation_Functions.R")

# ======================================================================================================================
# LOAD DATA
# ======================================================================================================================
load("Krackhardt_Multiplex_Directed.Rdata")

# ======================================================================================================================
# USE LATENTNET ON EACH NETWORK WITH G=1,...,7 USING BIC TO SELECT THE BEST ONE
# ======================================================================================================================
# set.seed(168956)
# M = N = 21
# bic_mat = matrix(NA, nrow = M, ncol = 3)
# for (m in 1:M) {
#   cat("\nNetwork", m, "\n---------\n")
#   for (g in 1:3) {
#     res_lpcm = latentnet::ergmm(network::as.network(krack_multi[, , m]) ~ euclidean(d = 2, G = g),
#                                 control = latentnet::ergmm.control(burnin = 1500000, sample.size = 5000, interval = 50))
#     sum_lpcm = summary(res_lpcm)
#     bic_mat[m, g] = sum_lpcm$bic$overall
#     cat("g =", g, ": Done\n")
#   } # end g for loop
# } # end m for loop
# 
# save(bic_mat, file = "Krackhardt_PPC_latentnet_res.Rdata")

# ======================================================================================================================
# FIND BEST FITTING LATENTNET MODEL FOR EACH NETWORK
# ======================================================================================================================
load("Krackhardt_PPC_latentnet_res.Rdata")
K_opt = apply(latentnet_res$bic_mat, 1, which.min)

# ======================================================================================================================
# REFIT EACH BEST ONE AND EXTRACT LAST 100 ALPHAS AND ZS
# ======================================================================================================================
# best_fits = vector("list", length = M)
# for (m in 1:M) {
#   best_fits[[m]] = latentnet::ergmm(network::as.network(krack_multi[, , m]) ~ euclidean(d = 2, G = K_opt[m]),
#                                     control = latentnet::ergmm.control(burnin = 1500000, sample.size = 5000, interval = 50))
#   cat("\nNetwork", m, ": Done")
# } # end m for loop
# 
# latentnet_res = list(bic_mat = bic_mat, best_fits = best_fits, K_opt = K_opt)
# save(latentnet_res, file = "Krackhardt_PPC_latentnet_res.Rdata")

# ======================================================================================================================
# PLOT ALL THE LATENT SPACES
# ======================================================================================================================
load("Krackhardt_PPC_latentnet_res.Rdata")

M = N = length(latentnet_res$best_fits)

pdf("Krack_PPC_Plots/Krack_latentnet_Latent_Spaces.pdf", height = 8, width = 16)
par(mfrow = c(3, 7))
par(mar = c(1.75, 1.75, 1.5, 0.5) + 0.4) # bltr
for (m in 1:M) {
  plot(latentnet_res$best_fits[[m]], main = LETTERS[m], las = 1, print.formula = F, cex.main = 1.8, vertex.cex = 2,
       cluster.col = khroma::color("muted")(9))
  legend("bottomleft", legend = paste0("K = ", K_opt[m]), bty = "n")
} # end m for loop
dev.off()

# ======================================================================================================================
# GENERATE 500 NEW MULTIPLEXES FROM THE POSTERIOR PREDICTIVE DISTRIBUTION
# ======================================================================================================================
rm(list = ls())

setwd("/Users/cjclarke/Desktop/PhD/PhD Project 1/LaPCoM_Code/Applications/Krackhardt/")
load("Krackhardt_PPC_latentnet_res.Rdata")

M = N = length(latentnet_res$best_fits)

num_PPCs = 500

latentnet_alphas = matrix(NA, nrow = num_PPCs, ncol = M)
latentnet_latent_spaces = array(NA, dim = c(N, 2, num_PPCs, M))

for (m in 1:M) {
  fit_m = latentnet_res$best_fits[[m]]
  latentnet_alphas[, m] = tail(fit_m$sample$beta, num_PPCs)
  latentnet_latent_spaces[, , , m] = fit_m$sample$Z[tail(1:dim(fit_m$sample$Z)[1], num_PPCs), , ]
} # end m for loop

set.seed(45741)

probs = array(NA, dim = c(N, N, M, num_PPCs))
multiplexes = array(NA, dim = c(N, N, M, num_PPCs))

for (ppc_i in 1:num_PPCs) {
  for (m in 1:M) {
    D = as.matrix(Rfast::Dist(latentnet_latent_spaces[, , ppc_i, m], square = T)) # distance matrix for Z_g
    eta = latentnet_alphas[ppc_i, m] - D
    probs[, , m, ppc_i] = (exp(eta)) / (1 + exp(eta))
    diag(probs[, , m, ppc_i]) = 0
    multiplexes[, , m, ppc_i] = matrix(rbinom(N^2, 1, probs[, , m, ppc_i]), ncol = N, nrow = N)
  } # end m for loop
} # end i for loop 

latentnet_PPD_data = list(multiplexes = multiplexes, probs = probs)
save(latentnet_PPD_data, file = "latentnet_PPD_data.Rdata")

# ======================================================================================================================
# RUN PPCS
# ======================================================================================================================
source("LaPCoM_Initialisation_Functions.R")
load("Krackhardt_Multiplex.Rdata")

square_diff_densities = matrix(NA, nrow = M, ncol = num_PPCs)
aucs_probs = matrix(NA, nrow = M, ncol = num_PPCs)
aucs_nets = matrix(NA, nrow = M, ncol = num_PPCs)
net_dists = matrix(NA, nrow = M, ncol = num_PPCs)
ham_dists = matrix(NA, nrow = M, ncol = num_PPCs)
F1s = matrix(NA, nrow = M, ncol = num_PPCs)
  
for (ppc_i in 1:num_PPCs) {
  
  prog_bar = progress::progress_bar$new(format = paste0("Multiplex ", ppc_i, ": [:bar] :percent"), total = num_PPCs, clear = F, width = 60)
  
  for (m in 1:M) {
    
    # ----------------------------------------------------------------------------------------------------------------
    # DENSITY
    # ----------------------------------------------------------------------------------------------------------------
    true_dens = network::network.density(network::as.network(krack_multi[, , m]))
    synth_dens = network::network.density(network::as.network(multiplexes[, , m, ppc_i]))
    square_diff_densities[m, ppc_i] = (true_dens - synth_dens)^2
    
    # ----------------------------------------------------------------------------------------------------------------
    # PRECISION-RECALL AUC WITH PROBABILITIES
    # ----------------------------------------------------------------------------------------------------------------
    pred_obj = ROCR::prediction(as.vector(probs[, , m, ppc_i]), as.vector(krack_multi[, , m]))
    pr_obj = ROCR::performance(pred_obj, "prec", "rec")
    # plot(pr_obj)
    aucs_probs[m, ppc_i] = ROCR::performance(pred_obj, "aucpr")@y.values[[1]][1]
    rm(pred_obj, pr_obj)
    
    # ----------------------------------------------------------------------------------------------------------------
    # PRECISION-RECALL AUC WITH NETWORKS
    # ----------------------------------------------------------------------------------------------------------------
    pred_obj_nets = ROCR::prediction(as.vector(multiplexes[, , m, ppc_i]), as.vector(krack_multi[, , m]))
    pr_obj_nets = ROCR::performance(pred_obj_nets, "prec", "rec")
    # plot(pr_obj)
    aucs_nets[m, ppc_i] = ROCR::performance(pred_obj_nets, "aucpr")@y.values[[1]][1]
    
    # ----------------------------------------------------------------------------------------------------------------
    # NETWORK-DISTANCE
    # ----------------------------------------------------------------------------------------------------------------
    if (sum(multiplexes[, , m, ppc_i]) == 0) {
      net_dists[m, ppc_i] = NA
    } else {
      gra1 = igraph::graph_from_adjacency_matrix(krack_multi[, , m], "directed")
      gra2 = igraph::graph_from_adjacency_matrix(multiplexes[, , m, ppc_i], "directed")
      el1 = cbind(igraph::as_edgelist(gra1))
      el2 = cbind(igraph::as_edgelist(gra2))
      net_dists[m, ppc_i] = nature_dist(el1, el2, gra1, gra2, w1 = 0.45, w2 = 0.45, w3 = 0.1, net_type = "binary")
    } # end if else statement
    
    # ----------------------------------------------------------------------------------------------------------------
    # HAMMING-DISTANCE
    # ----------------------------------------------------------------------------------------------------------------
    net_array = array(NA, dim = c(N, N, 2))
    net_array[, , 1] = krack_multi[, , m]
    net_array[, , 2] = multiplexes[, , m, ppc_i]
    
    ham_dists[m, ppc_i] = sna::hdist(net_array, g1 = 1, g2 = 2, normalize = T)
    
    # ----------------------------------------------------------------------------------------------------------------
    # F1-SCORE
    # ----------------------------------------------------------------------------------------------------------------
    F1s[m, ppc_i] = MLmetrics::F1_Score(krack_multi[, , m], multiplexes[, , m, ppc_i])
    
  } # end m for loop
  
  prog_bar$tick()
  
} # end ppc_i for loop

beepr::beep(sound = 10)

save.image(file = "latentnet_PPC_res.Rdata")

# ======================================================================================================================
# PLOTS TO COMPARE TO LaPCoM RESULTS
# ======================================================================================================================
png("Krack_PPC_Plots/latentnet_Density.png", height = 420 * 2, width = 480 * 3)
par(mfrow = c(3, 7))
par(mar = c(5, 5, 4, 1) + 0.1) # bltr
for (m in 1:M) { 
  boxplot(as.vector(square_diff_densities[m, ]), xlab = "", ylab = "", ylim = c(0, 0.4),
          main = paste("Network", LETTERS[m]), las = 1, cex.axis = 2, cex.lab = 1.5, cex.main = 2) 
  legend("top", legend = paste("Median =", round(median(square_diff_densities[m, ]), 3)), bty = "n", cex = 1.5)
} # end m for loop
dev.off()

png("Krack_PPC_Plots/latentnet_AUC_probs.png", height = 420 * 2, width = 480 * 3)
par(mfrow = c(3, 7))
par(mar = c(5, 5, 4, 1) + 0.1) # bltr
for (m in 1:M) {
  boxplot(as.vector(aucs_probs[m, ]), xlab = "", ylab = "", ylim = c(0, 1),
          main = paste("Network", LETTERS[m]), las = 1, cex.axis = 2, cex.lab = 1.5, cex.main = 2) 
  abline(h = network::network.density(network::as.network(krack_multi[, , m])), col = "red")
  legend("bottom", legend = paste("Median =", round(median(aucs_probs[m, ]), 2)), bty = "n", cex = 1.5)
} # end m for loop
dev.off()

png("Krack_PPC_Plots/latentnet_AUC_nets.png", height = 420 * 2, width = 480 * 3)
par(mfrow = c(3, 7))
par(mar = c(5, 5, 4, 1) + 0.1) # bltr
for (m in 1:M) {
  boxplot(as.vector(aucs_nets[m, ]), xlab = "", ylab = "", ylim = c(0, 1),
          main = paste("Network", LETTERS[m]), las = 1, cex.axis = 2, cex.lab = 1.5, cex.main = 2) 
  abline(h = network::network.density(network::as.network(krack_multi[, , m])), col = "red")
  legend("bottom", legend = paste("Median =", round(median(aucs_nets[m, ]), 2)), bty = "n", cex = 1.5)
} # end m for loop
dev.off()

png("Krack_PPC_Plots/latentnet_Net_Distances.png", height = 420 * 2, width = 480 * 3)
par(mfrow = c(3, 7))
par(mar = c(5, 5, 4, 1) + 0.1) # bltr
for (m in 1:M) {
  boxplot(as.vector(net_dists[m, ]), xlab = "", ylab = "", ylim = c(0, 0.7),
          main = paste("Network", LETTERS[m]), las = 1, cex.axis = 2, cex.lab = 1.5, cex.main = 2) 
  legend("top", legend = paste("Median =", round(median(net_dists[m, ]), 2)), bty = "n", cex = 1.5)
} # end m for loop
dev.off()

png("Krack_PPC_Plots/latentnet_Ham_Distances.png", height = 420 * 2, width = 480 * 3)
par(mfrow = c(3, 7))
par(mar = c(5, 5, 4, 1) + 0.1) # bltr
for (m in 1:M) {
  boxplot(as.vector(ham_dists[m, ]), xlab = "", ylab = "", ylim = c(0, 0.2),
          main = paste("Network", LETTERS[m]), las = 1, cex.axis = 2, cex.lab = 1.5, cex.main = 2) 
  legend("top", legend = paste("Median =", round(median(ham_dists[m, ]), 2)), bty = "n", cex = 1.5)
} # end m for loop
dev.off()

png("Krack_PPC_Plots/latentnet_F1s.png", height = 420 * 2, width = 480 * 3)
par(mfrow = c(3, 7))
par(mar = c(5, 5, 4, 1) + 0.1) # bltr
for (m in 1:M) {
  boxplot(as.vector(F1s[m, ]), xlab = "", ylab = "", ylim = c(0, 1),
          main = paste("Network", LETTERS[m]), las = 1, cex.axis = 2, cex.lab = 1.5, cex.main = 2) 
  legend("bottom", legend = paste("Median =", round(median(F1s[m, ]), 2)), bty = "n", cex = 1.5)
} # end m for loop
dev.off()
