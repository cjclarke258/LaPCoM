####################################################################################################
# KRACKHARDT DATA APPLICATION (POST-PROCESSING)
####################################################################################################

rm(list = ls())
setwd("/Volumes/My Passport/PhD/PhD Project 1/LaPCoM_Code/Applications/Krackhardt/")

# ==================================================================================================
# LOAD IN ALL NECESSARY FILES AND FUNCTIONS
# ==================================================================================================
source("LaPCoM_PP_MethodB.R")

# ==================================================================================================
# LOAD IN THE DATA
# ==================================================================================================
num_chains = 40
output_combined = vector("list", length = num_chains)
for (chain in 1:num_chains) {
  load(paste0("/Volumes/My Passport/PhD/PhD Project 1/LaPCoM_Code/Applications/Krackhardt/Krack_Output/krack_output_sim_", chain, ".Rdata"))
  output_combined[[chain]] = output
  output_combined[[chain]]$sim_multi$multi = output_combined[[chain]]$krack_multi
  output_combined[[chain]]$res_LaPCoM = output_combined[[chain]]$res_krack
  rm(output)
} # end chain for loop
rm(chain)

# ------------------------------------------------------------------------------------------------
# CALCULATE AVERAGE ACCEPTANCE RATES ACROSS THE REPLICATIONS
# ------------------------------------------------------------------------------------------------
avg_ARs = vector("list", length = 0)

avg_AR_alpha = avg_AR_e = rep(NA, num_chains)
avg_AR_Z = avg_AR_w = matrix(NA, nrow = num_chains, ncol = 12)

for (chain in 1:num_chains) {
  avg_AR_alpha[chain] = mean(output_combined[[chain]]$res_LaPCoM$store_ARs$acc_rate_alpha)
  temp_Z = rowMeans(output_combined[[chain]]$res_LaPCoM$store_ARs$acc_rate_Z)
  avg_AR_Z[chain, ] = c(temp_Z, rep(NA, 12 - length(temp_Z)))
  avg_AR_e[chain] = mean(output_combined[[chain]]$res_LaPCoM$store_ARs$acc_rate_e)
  temp_w = rowMeans(output_combined[[chain]]$res_LaPCoM$store_ARs$acc_rate_w)
  avg_AR_w[chain, ] = c(temp_w, rep(NA, 12 - length(temp_w)))
} # end sim for loop

avg_ARs$alpha = mean(avg_AR_alpha)
avg_ARs$Z = colMeans(avg_AR_Z, na.rm = T)
avg_ARs$e = mean(avg_AR_e)
avg_ARs$w = colMeans(avg_AR_w, na.rm = T)

avg_ARs

# ------------------------------------------------------------------------------------------------
# POST-PROCESSING
# ------------------------------------------------------------------------------------------------
pp_res_krack = post_process_LaPCoM(output = output_combined, 
                               multi = output_combined[[1]]$krack_multi,
                               plots_path = "Krack_PP_Plots_PPMethodB/")

posterior_means_modes = pp_res_krack$posterior_means_modes

# --------------------------------------------------------------------------------------------------------------------
# ONLY CONSIDER CHAINS SUCH THAT G_+^ = G_+^OVERALL
# --------------------------------------------------------------------------------------------------------------------
barplot(table(posterior_means_modes$G_plus_hat_FS), las = 1, ylim = c(0, 35), ylab = "Freq.", xlab = "G_+_hat")
barplot(table(posterior_means_modes$G_plus_hat_WG), las = 1, ylim = c(0, 35), ylab = "Freq.", xlab = "G_+_hat")
G_plus_hat_all_chains = as.numeric(DescTools::Mode(posterior_means_modes$G_plus_hat_WG))
G_plus_hat_all_chains

chains_keep = which(posterior_means_modes$G_plus_hat_WG == G_plus_hat_all_chains)
num_chains_keep = length(chains_keep)

# --------------------------------------------------------------------------------------------------------------------
# SELECT THE "BEST" CHAIN, I.E. THE ONE THAT MAXIMISES THE LOG POSTERIOR
# --------------------------------------------------------------------------------------------------------------------
Rcpp::sourceCpp("LaPCoM_log_LPM_fast.cpp")
source("LaPCoM_TS_Full_Conditionals.R")

multi = output_combined[[1]]$krack_multi
net_type = "binary"
net_mode = "undirected"
M = dim(multi)[3]
N = dim(multi)[1]

G_0 = 2
G_max = 5

m_alpha = 0; s_alpha = 1
a_G = 8; b_G = 18; c_G = 10
l_G = 6; r_G = 3
s_e = 4

K_0 = 2
K_max = 6
a_K = 8; b_K = 18; c_K = 10
l_K = 6; r_K = 3
s_w = 3

mu_mean = rep(0, 2); mu_cov = diag(2)
u_sigma = 11; v_sigma = 2

store_log_post = rep(NA, num_chains_keep)

for (chain_num in 1:num_chains_keep) {
  
  chain_index = chains_keep[chain_num]
  
  store_log_post[chain_num] = log_post(multi, alpha = posterior_means_modes$alpha[[chain_index]], 
                                       Z = posterior_means_modes$Z[[chain_index]], 
                                       C_vector = posterior_means_modes$C[chain_index, ], 
                                       net_type = net_type, 
                                       net_mode = net_mode,
                                       G_plus = G_plus_hat_all_chains, 
                                       tau = posterior_means_modes$tau[[chain_index]], 
                                       e = posterior_means_modes$e[[chain_index]], 
                                       G = posterior_means_modes$G_hat[[chain_index]], 
                                       a_G = a_G, b_G = b_G, c_G = c_G, l_G = l_G, r_G = r_G, 
                                       K_plus = posterior_means_modes$Kg_plus_hat_WG[[chain_index]], 
                                       K = posterior_means_modes$Kg_hat[[chain_index]], 
                                       S_vector = posterior_means_modes$S[[chain_index]], 
                                       pi = posterior_means_modes$pi[[chain_index]], 
                                       mu = posterior_means_modes$mu[[chain_index]], 
                                       Sigma = posterior_means_modes$Sigma[[chain_index]], 
                                       mu_mean = mu_mean, mu_cov = mu_cov, u_sigma = u_sigma, v_sigma = v_sigma, 
                                       w = posterior_means_modes$w[[chain_index]], 
                                       a_K = a_K, b_K = b_K, c_K = c_K, l_K = l_K, r_K = r_K)
  
} # end chain_num for loop
rm(chain_num)

best_chain_ind = which.max(store_log_post)
max(store_log_post)
best_chain_num = chains_keep[best_chain_ind]
best_chain_num

beepr::beep(sound = 10)

# --------------------------------------------------------------------------------------------------------------------
# OBTAIN CLUSTERING FROM THE "BEST" CHAIN
# --------------------------------------------------------------------------------------------------------------------
best_chain_clustering = posterior_means_modes$C[best_chain_num, ]
best_chain_clustering

# --------------------------------------------------------------------------------------------------------------------
# COMPARE IT TO THE OTHER CLUSTERINGS OF THE KEPT CHAINS TO CHECK STABILITY
# --------------------------------------------------------------------------------------------------------------------
apply(posterior_means_modes$C[chains_keep[-best_chain_ind], ], 1, 
      function(col) { mclust::adjustedRandIndex(best_chain_clustering, col) })

# --------------------------------------------------------------------------------------------------------------------
# OBTAIN POSTERIOR ESTIMATES FOR THIS "BEST" CHAIN (WITHOUT AGGREGATING ACROSS CHAINS)
# --------------------------------------------------------------------------------------------------------------------
krack_alpha = posterior_means_modes$alpha[best_chain_num]
# posterior_means_modes$e
# posterior_means_modes$tau
krack_Z = posterior_means_modes$Z[[best_chain_num]]
krack_mu = posterior_means_modes$mu[[best_chain_num]]
krack_Sigma = posterior_means_modes$Sigma[[best_chain_num]]
krack_S = posterior_means_modes$S[[best_chain_num]]
krack_Kg_plus = posterior_means_modes$Kg_plus_hat_WG[[best_chain_num]]
# posterior_means_modes$pi
# posterior_means_modes$w
krack_C = best_chain_clustering
krack_G = G_plus_hat_all_chains

krack_pp = list(krack_alpha = krack_alpha,
                krack_Z = krack_Z,
                krack_mu = krack_mu,
                krack_Sigma = krack_Sigma,
                krack_S = krack_S,
                krack_Kg_plus = krack_Kg_plus,
                krack_C = krack_C,
                best_chain_num = best_chain_num, 
                krack_G = krack_G,
                pp_res_krack = pp_res_krack, 
                avg_ARs = avg_ARs, 
                multi = multi)

save(krack_pp, file = "pp_res_krack_PPMethodB.Rdata")  