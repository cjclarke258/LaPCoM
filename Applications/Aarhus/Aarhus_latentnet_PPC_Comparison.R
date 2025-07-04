########################################################################################################################
# AARHUS DATA APPLICATION (POSTERIOR PREDICTIVE CHECKS) [USING LATENTNET AS A COMPARISON]
########################################################################################################################

rm(list = ls())

setwd("/Volumes/My Passport/PhD/PhD Project 1/LaPCoM_Code/Applications/Aarhus/")
source("LaPCoM_Initialisation_Functions.R")

# ======================================================================================================================
# LOAD DATA
# ======================================================================================================================
aarhus_multi = readRDS("aarhus_adjacency_matrices_updated.rds")

# ======================================================================================================================
# USE LATENTNET ON EACH NETWORK WITH G=1,...,7 USING BIC TO SELECT THE BEST ONE
# ======================================================================================================================
set.seed(3668142)
N = dim(aarhus_multi)[1]
M = dim(aarhus_multi)[3]
G_max = 7
bic_mat = matrix(NA, nrow = M, ncol = G_max)
for (m in 1:M) {
  cat("\nNetwork", m, "\n---------\n")
  for (g in 1:G_max) {
    res_lpcm = latentnet::ergmm(network::as.network(aarhus_multi[, , m]) ~ euclidean(d = 2, G = g),
                                control = latentnet::ergmm.control(burnin = 1500000, sample.size = 5000, interval = 50))
    sum_lpcm = summary(res_lpcm)
    bic_mat[m, g] = sum_lpcm$bic$overall
    cat("g =", g, ": Done\n")
  } # end g for loop
} # end m for loop

save(bic_mat, file = "Aarhus_PPC_latentnet_res.Rdata")

# ======================================================================================================================
# FIND BEST FITTING LATENTNET MODEL FOR EACH NETWORK
# ======================================================================================================================
K_opt = apply(bic_mat, 1, which.min)

# ======================================================================================================================
# REFIT EACH BEST ONE AND EXTRACT LAST 100 ALPHAS AND ZS
# ======================================================================================================================
best_fits = vector("list", length = M)
for (m in 1:M) {
  best_fits[[m]] = latentnet::ergmm(network::as.network(aarhus_multi[, , m]) ~ euclidean(d = 2, G = K_opt[m]),
                                    control = latentnet::ergmm.control(burnin = 1500000, sample.size = 5000, interval = 50))
  cat("\nNetwork", m, ": Done")
} # end m for loop

latentnet_res = list(bic_mat = bic_mat, best_fits = best_fits, K_opt = K_opt)
save(latentnet_res, file = "Aarhus_PPC_latentnet_res.Rdata")

beepr::beep(sound = 10)

# ======================================================================================================================
# PLOT ALL THE LATENT SPACES
# ======================================================================================================================
rm(list = ls())

setwd("/Volumes/My Passport/PhD/PhD Project 1/LaPCoM_Code/Applications/Aarhus/")

aarhus_multi = readRDS("aarhus_adjacency_matrices_updated.rds")
M = dim(aarhus_multi)[3]
N = dim(aarhus_multi)[1]

load("Aarhus_PPC_latentnet_res.Rdata")

pdf("Aarhus_PPC_Plots/Aarhus_latentnet_Latent_Spaces.pdf", height = 2, width = 10)
par(mfrow = c(1, 5))
par(mar = c(1.75, 1.75, 1.5, 0.5) + 0.4) # bltr
for (m in 1:M) {
  plot(latentnet_res$best_fits[[m]], main = dimnames(aarhus_multi)[[3]][m], las = 1, print.formula = F, cex.main = 1.8, vertex.cex = 1.75,
       cluster.col = khroma::color("muted")(9))
  legend("bottomright", legend = paste0("K = ", latentnet_res$K_opt[m]), bty = "n")
} # end m for loop
dev.off()

# ======================================================================================================================
# GENERATE 500 NEW MULTIPLEXES FROM THE POSTERIOR PREDICTIVE DISTRIBUTION
# ======================================================================================================================
rm(list = ls())

setwd("/Volumes/My Passport/PhD/PhD Project 1/LaPCoM_Code/Applications/Aarhus/")

aarhus_multi = readRDS("aarhus_adjacency_matrices_updated.rds")
M = dim(aarhus_multi)[3]
N = dim(aarhus_multi)[1]

load("Aarhus_PPC_latentnet_res.Rdata")

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

for (ppc_i in 1:num_PPCs) {
  for (m in 1:M) {
    D = as.matrix(Rfast::Dist(latentnet_latent_spaces[, , ppc_i, m], square = T)) # distance matrix for Z_g
    eta = latentnet_alphas[ppc_i, m] - D
    # probs[, , m, ppc_i] = (exp(eta)) / (1 + exp(eta))
    probs[, , m, ppc_i][upper.tri(probs[, , m, ppc_i])] = (exp(eta[upper.tri(eta)])) / (1 + exp(eta[upper.tri(eta)]))
    probs[, , m, ppc_i] = probs[, , m, ppc_i] + t(probs[, , m, ppc_i])
    diag(probs[, , m, ppc_i]) = 0
    # multiplexes[, , m, ppc_i] = matrix(rbinom(N^2, 1, probs[, , m, ppc_i]), ncol = N, nrow = N)
    multiplexes[, , m, ppc_i][upper.tri(multiplexes[, , m, ppc_i])] = rbinom(((N * (N - 1)) / 2), 1, probs[, , m, ppc_i][upper.tri(probs[, , m, ppc_i])])
    multiplexes[, , m, ppc_i] = multiplexes[, , m, ppc_i] + t(multiplexes[, , m, ppc_i])
  } # end m for loop
} # end i for loop 

latentnet_PPD_data = list(multiplexes = multiplexes, probs = probs)
save(latentnet_PPD_data, file = "Aarhus_latentnet_PPD_data.Rdata")

beepr::beep(sound = 10)

# ======================================================================================================================
# RUN PPCS
# ======================================================================================================================
source("LaPCoM_Initialisation_Functions.R")

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
    true_dens = network::network.density(network::as.network(aarhus_multi[, , m]))
    synth_dens = network::network.density(network::as.network(multiplexes[, , m, ppc_i]))
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
    pred_obj_nets = ROCR::prediction(as.vector(multiplexes[, , m, ppc_i]), as.vector(aarhus_multi[, , m]))
    pr_obj_nets = ROCR::performance(pred_obj_nets, "prec", "rec")
    # plot(pr_obj)
    aucs_nets[m, ppc_i] = ROCR::performance(pred_obj_nets, "aucpr")@y.values[[1]][1]
    
    # ----------------------------------------------------------------------------------------------------------------
    # NETWORK-DISTANCE
    # ----------------------------------------------------------------------------------------------------------------
    if (sum(multiplexes[, , m, ppc_i]) == 0) {
      net_dists[m, ppc_i] = NA
    } else {
      gra1 = igraph::graph_from_adjacency_matrix(aarhus_multi[, , m], "directed")
      gra2 = igraph::graph_from_adjacency_matrix(multiplexes[, , m, ppc_i], "directed")
      el1 = cbind(igraph::as_edgelist(gra1))
      el2 = cbind(igraph::as_edgelist(gra2))
      net_dists[m, ppc_i] = nature_dist(el1, el2, gra1, gra2, w1 = 0.45, w2 = 0.45, w3 = 0.1, net_type = "binary")
    } # end if else statement
    
    # ----------------------------------------------------------------------------------------------------------------
    # HAMMING-DISTANCE
    # ----------------------------------------------------------------------------------------------------------------
    net_array = array(NA, dim = c(N, N, 2))
    net_array[, , 1] = aarhus_multi[, , m]
    net_array[, , 2] = multiplexes[, , m, ppc_i]
    
    ham_dists[m, ppc_i] = sna::hdist(net_array, g1 = 1, g2 = 2, normalize = T)
    
    # ----------------------------------------------------------------------------------------------------------------
    # F1-SCORE
    # ----------------------------------------------------------------------------------------------------------------
    F1s[m, ppc_i] = MLmetrics::F1_Score(aarhus_multi[, , m], multiplexes[, , m, ppc_i])
    
  } # end m for loop
  
  prog_bar$tick()
  
} # end ppc_i for loop

beepr::beep(sound = 10)

save.image(file = "Aarhus_latentnet_PPC_res.Rdata")

# ======================================================================================================================
# PLOTS TO COMPARE TO LaPCoM RESULTS
# ======================================================================================================================
pdf("Aarhus_PPC_Plots/latentnet_Density.pdf", height = 2.75, width = 10)
par(mfrow = c(1, 5))
par(mar = c(4, 4, 2, 0.5) + 0.4) # bltr
for (m in 1:M) { 
  boxplot(as.vector(square_diff_densities[m, ]), 
          xlab = paste("Median =", round(median(square_diff_densities[m, ]), 3)), ylab = "", ylim = c(0, 0.4),
          main = dimnames(aarhus_multi)[[3]][m], las = 1, cex.axis = 1, cex.lab = 1.2, cex.main = 1.5) 
  # legend("top", legend = paste("Median =", round(median(square_diff_densities[m, ]), 3)), bty = "n", cex = 1.5)
} # end m for loop
dev.off()

pdf("Aarhus_PPC_Plots/latentnet_AUC_probs.pdf", height = 2.75, width = 10)
par(mfrow = c(1, 5))
par(mar = c(4, 4, 2, 0.5) + 0.4) # bltr
for (m in 1:M) {
  boxplot(as.vector(aucs_probs[m, ]), 
          xlab = paste("Median =", round(median(aucs_probs[m, ]), 2)), ylab = "", ylim = c(0, 1),
          main = dimnames(aarhus_multi)[[3]][m], las = 1, cex.axis = 1, cex.lab = 1.2, cex.main = 1.5) 
  abline(h = network::network.density(network::as.network(aarhus_multi[, , m])), col = "red")
  # legend("bottom", legend = paste("Median =", round(median(aucs_probs[m, ]), 2)), bty = "n", cex = 1.5)
} # end m for loop
dev.off()

pdf("Aarhus_PPC_Plots/latentnet_AUC_nets.pdf", height = 2.75, width = 10)
par(mfrow = c(1, 5))
par(mar = c(4, 4, 2, 0.5) + 0.4) # bltr
for (m in 1:M) {
  boxplot(as.vector(aucs_nets[m, ]), 
          xlab = paste("Median =", round(median(aucs_nets[m, ]), 2)), ylab = "", ylim = c(0, 1),
          main = dimnames(aarhus_multi)[[3]][m], las = 1, cex.axis = 1, cex.lab = 1.2, cex.main = 1.5) 
  abline(h = network::network.density(network::as.network(aarhus_multi[, , m])), col = "red")
  # legend("bottom", legend = paste("Median =", round(median(aucs_nets[m, ]), 2)), bty = "n", cex = 1.5)
} # end m for loop
dev.off()

pdf("Aarhus_PPC_Plots/latentnet_Net_Distances.pdf", height = 2.75, width = 10)
par(mfrow = c(1, 5))
par(mar = c(4, 4, 2, 0.5) + 0.4) # bltr
for (m in 1:M) {
  boxplot(as.vector(net_dists[m, ]), 
          xlab = paste("Median =", round(median(net_dists[m, ]), 2)), ylab = "", ylim = c(0, 0.7),
          main = dimnames(aarhus_multi)[[3]][m], las = 1, cex.axis = 1, cex.lab = 1.2, cex.main = 1.5) 
  # legend("top", legend = paste("Median =", round(median(net_dists[m, ]), 2)), bty = "n", cex = 1.5)
} # end m for loop
dev.off()

pdf("Aarhus_PPC_Plots/latentnet_Ham_Distances.pdf", height = 2.75, width = 10)
par(mfrow = c(1, 5))
par(mar = c(4, 4, 2, 0.5) + 0.4) # bltr
for (m in 1:M) {
  boxplot(as.vector(ham_dists[m, ]), 
          xlab = paste("Median =", round(median(ham_dists[m, ]), 2)), ylab = "", ylim = c(0, 0.2),
          main = dimnames(aarhus_multi)[[3]][m], las = 1, cex.axis = 1, cex.lab = 1.2, cex.main = 1.5) 
  # legend("top", legend = paste("Median =", round(median(ham_dists[m, ]), 2)), bty = "n", cex = 1.5)
} # end m for loop
dev.off()

pdf("Aarhus_PPC_Plots/latentnet_F1s.pdf", height = 2.75, width = 10)
par(mfrow = c(1, 5))
par(mar = c(4, 4, 2, 0.5) + 0.4) # bltr
for (m in 1:M) {
  boxplot(as.vector(F1s[m, ]), 
          xlab = paste("Median =", round(median(F1s[m, ]), 2)), ylab = "", ylim = c(0, 1),
          main = dimnames(aarhus_multi)[[3]][m], las = 1, cex.axis = 1, cex.lab = 1.2, cex.main = 1.5) 
  # legend("bottom", legend = paste("Median =", round(median(F1s[m, ]), 2)), bty = "n", cex = 1.5)
} # end m for loop
dev.off()
