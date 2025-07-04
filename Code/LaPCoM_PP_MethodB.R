########################################################################################################################
# FUNCTION TO PERFORM POST-PROCESSING ON THE RESULTING MCMC CHAIN FROM LaPCoM
########################################################################################################################
post_process_LaPCoM = function(output, # a list containing the output of the LaPCoM algorithm for each simulation
                           multi, # the multiplex in array form with dimensions (N, N, M)
                           produce_plots = NULL,
                           plots_path = NULL, # name of folder to store plots
                           cols = NULL, 
                           seed = NULL) {
  
  # ====================================================================================================================
  # NECESSARY CHECKS
  # ====================================================================================================================
  if (missing(output) || !is.list(output) || length(output) < 1) {
    stop("output must be a list of length at least 1 and must be provided.")
  } # end if statement
  
  if (missing(multi) || !is.array(multi)) {
    stop("multi must be an array and must be provided.")
  } # end if statement
  
  if (isTRUE(produce_plots)) {
    
    if (is.null(plots_path)) {
      stop("Please provide a valid 'plots_path' when 'produce_plots' is TRUE.")
    } else {
      if (!dir.exists(plots_path)) { # Create the directory if it doesn't exist
        dir.create(plots_path, recursive = TRUE)
      } # end if statement
    } # end if else statement
    
  } else if (is.null(produce_plots) && !is.null(plots_path)) {
    
    produce_plots = TRUE
    if (!dir.exists(plots_path)) { # Create the directory if it doesn't exist
      dir.create(plots_path, recursive = TRUE)
    } # end if statement
    
  } else if (is.null(produce_plots) && is.null(plots_path)) {
    produce_plots = FALSE
    
  } # end if else statement
  
  if (is.null(cols)) { # check for cols; if NULL, set to default
    cols = rep(khroma::color("light")(9), 4)
    message("cols not specified; using default colour palette.")
  } # end if statement
  
  # ====================================================================================================================
  # SET THE SEED
  # ====================================================================================================================
  set.seed(seed)
  
  # ====================================================================================================================
  # POST PROCESSING PROCEDURE FOR LaPCoM
  # ====================================================================================================================
  M = dim(multi)[3]
  N = dim(multi)[1]
  
  G_max = dim(output[[1]]$res_LaPCoM$store_Zs)[3]
  K_max = dim(output[[1]]$res_LaPCoM$store_mus)[1]
  
  params = list(M = M, N = N, G_max = G_max, K_max = K_max)
  params$Kg_plus_hat = vector("list", length = length(output))
  params$T0_iters_net = vector("list", length = length(output)) # iterations where G+ = G+hat
  params$perm_indices_net = vector("list", length = length(output)) # classifications/iterations that were permutations
  params$T0_rho_net = rep(NA, length(output)) # number of class. seq.s / iterations that were not permutations
  params$perm_indices_node = vector("list", length = length(output)) # classifications/iterations that were permutations
  params$T0_rho_node = vector("list", length = length(output)) # number of class. seq.s / iterations that were not permutations
  
  chains = list(alpha = vector("list", length = length(output)), 
                C_vector = vector("list", length = length(output)),
                G = vector("list", length = length(output)),
                G_plus = vector("list", length = length(output)),
                e = vector("list", length = length(output)),
                tau = vector("list", length = length(output)),
                Z = vector("list", length = length(output)), 
                proc_matched_Z = vector("list", length = length(output)),
                mu = vector("list", length = length(output)),
                Sigma = vector("list", length = length(output)),
                S_vector = vector("list", length = length(output)),
                Kg = vector("list", length = length(output)),
                Kg_plus = vector("list", length = length(output)),
                pi = vector("list", length = length(output)),
                w = vector("list", length = length(output)))
  
  acc_rates = list(alpha = rep(NA, length(output)),
                   e = rep(NA, length(output)),
                   Z = vector("list", length = length(output)),
                   w = vector("list", length = length(output)))
  
  posterior_means_modes = list(alpha = rep(NA, length(output)),
                               C = matrix(NA, nrow = length(output), ncol = M),
                               G_hat = vector("list", length = length(output)),
                               G_plus_hat_WG = rep(NA, length(output)),
                               G_plus_hat_FS = rep(NA, length(output)),
                               e = rep(NA, length(output)),
                               tau = vector("list", length = length(output)),
                               Z = vector("list", length = length(output)),
                               mu = vector("list", length = length(output)),
                               Sigma = vector("list", length = length(output)),
                               S = vector("list", length = length(output)),
                               Kg_hat = vector("list", length = length(output)),
                               Kg_plus_hat_WG = vector("list", length = length(output)),
                               Kg_plus_hat_FS = vector("list", length = length(output)),
                               pi = vector("list", length = length(output)),
                               w = vector("list", length = length(output)))
  
  for (sim in 1:length(output)) {
    
    # ------------------------------------------------------------------------------------------------------------------
    # NETWORK-LEVEL PARAMETERS
    # ------------------------------------------------------------------------------------------------------------------
    
    T_iters_net = dim(output[[sim]]$res_LaPCoM$store_C_vector)[2] # number of iterations to be post-processed
    
    # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    # 1. USE THE WADE AND GHARAMANI METHOD TO OBTAIN THE OPTIMAL NETWORK-LEVEL CLUSTERING (USING ALL POSTERIOR SAMPLES)
    # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    store_C_vector = output[[sim]]$res_LaPCoM$store_C_vector
    
    psm_net = mcclust::comp.psm(t(store_C_vector))
    wade_ghar_net = mcclust.ext::minVI(psm_net, method = "all", cls.draw = t(store_C_vector), include.greedy = F)
    
    posterior_means_modes$C[sim, ] = wade_ghar_net$cl[1, ]
    
    # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    # 2. DETERMINE THE NUMBER OF NON-EMPTY COMPONENTS
    # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    G_plus_hat = posterior_means_modes$G_plus_hat_WG[sim] = length(unique(posterior_means_modes$C[sim, ]))
    
    # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    # 3. RESULTS IN A CLASSIFICATION INDICATING TO WHICH CLUSTERS THE COMPONENT-SPECIFIC PARAMS OF EACH DRAW BELONG
    # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    
    # ..............................................................................................................
    # (A) ISOLATE THE NON-EMPTY COMPONENTS
    # ..............................................................................................................
    non_empty_net_clusts = vector("list", length = T_iters_net)
    for (t in 1:T_iters_net) { non_empty_net_clusts[[t]] = sort(unique(store_C_vector[, t])) } # end t for loop
    
    store_tau_active = store_Zs_active = store_mus_active = store_Sigmas_active = store_S_vector_active = 
      store_Kg_active = store_Kg_plus_active = store_pi_active = store_w_active = acc_rate_Z_active = 
      acc_rate_w_active = vector("list", length = T_iters_net)
    
    for (t in 1:T_iters_net) { 
      
      store_tau_active[[t]] = output[[sim]]$res_LaPCoM$store_tau[non_empty_net_clusts[[t]], t]
      store_Zs_active[[t]] = output[[sim]]$res_LaPCoM$store_Zs[, , non_empty_net_clusts[[t]], t]
      store_mus_active[[t]] = output[[sim]]$res_LaPCoM$store_mus[, , non_empty_net_clusts[[t]], t]
      store_Sigmas_active[[t]] = output[[sim]]$res_LaPCoM$store_Sigmas[, , , non_empty_net_clusts[[t]], t]
      store_S_vector_active[[t]] = output[[sim]]$res_LaPCoM$store_S_vector[, non_empty_net_clusts[[t]], t]
      store_Kg_active[[t]] = output[[sim]]$res_LaPCoM$store_K[non_empty_net_clusts[[t]], t]
      store_Kg_plus_active[[t]] = output[[sim]]$res_LaPCoM$store_K_plus[non_empty_net_clusts[[t]], t]
      store_pi_active[[t]] = output[[sim]]$res_LaPCoM$store_pi[, non_empty_net_clusts[[t]], t]
      store_w_active[[t]] = output[[sim]]$res_LaPCoM$store_w[non_empty_net_clusts[[t]], t]
      acc_rate_Z_active[[t]] = output[[sim]]$res_LaPCoM$store_ARs$acc_rate_Z[non_empty_net_clusts[[t]], t]
      acc_rate_w_active[[t]] = output[[sim]]$res_LaPCoM$store_ARs$acc_rate_w[non_empty_net_clusts[[t]], t]
      
    } # end t for loop
    
    # ..............................................................................................................
    # (B) USING ONE (ANY) LS AS A REFERENCE, PROCRUSTES TRANSFORM ALL OTHERS TO ALIGN SO THAT THEY ARE COMPARABLE 
    # ..............................................................................................................
    if (G_plus_hat == 1) {
      
      proc_matched_Z_temp = vector("list", length = length(store_Zs_active))
      ref_space = store_Zs_active[[1]]
      
      for (t in 1:T_iters_net) {
        proc_matched_Z_temp[[t]] = vegan::procrustes(ref_space, store_Zs_active[[t]])$Yrot
      } # end t for loop
      
    } else { # G_plus_hat > 1
      
      proc_matched_Z_temp = vector("list", length = length(store_Zs_active))
      ref_space = store_Zs_active[[1]][, , 1]
      
      for (t in 1:T_iters_net) {
        
        proc_matched_Z_temp[[t]] = array(NA, dim = dim(store_Zs_active[[t]]))
        
        for (g in 1:dim(store_Zs_active[[t]])[3]) {
          proc_matched_Z_temp[[t]][, , g] = vegan::procrustes(ref_space, store_Zs_active[[t]][, , g])$Yrot
        } # end g for loop
        
      } # end t for loop
      
    } # end if else statement
    
    # ..............................................................................................................
    # (C) ARRANGE THE LATENT SPACES INTO A MATRIX WITH G_PLUS_HAT X T ROWS AND d*N COLUMNS
    # ..............................................................................................................
    if (G_plus_hat > 1) { # only needs to be done if G_plus_hat > 1, otherwise, already in that format
      
      ls_vector = matrix(NA, nrow = 0, ncol = 2*N)
      for (t in 1:T_iters_net) {
        
        # row_means = rbind(row_means, t(apply(store_Zs_active[[t]], 3, rowMeans))) # original code used row means
        ls_vector = rbind(ls_vector, t(apply(proc_matched_Z_temp[[t]], 3, function(mat) { apply(mat, 2, as.vector) })))
        
      } # end t for loop 
      
    } # end if else statement
    
    # ..............................................................................................................
    # (D) CLUSTER THE MATRIX INTO G_PLUS_HAT CLUSTERS USING K-MEANS CLUSTERING
    # ..............................................................................................................
    if (G_plus_hat > 1) { # only needs to be done if G_plus_hat > 1, otherwise, not necessary
      
      kmeans_res_net = kmeans(scale(ls_vector), centers = G_plus_hat, nstart = 100, iter.max = 100)
      
      png(paste0(plots_path, "/Sim", sim, "_Net-Level_PPR.png"))
      plot(ls_vector, pch = 19, las = 1, xlab = "Dim. 1", ylab = "Dim. 2", col = kmeans_res_net$cluster,
           main = "Point Process Representation at the Network-Level")
      dev.off()
      
    } # end if else statement
    
    # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    # 4.
    # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    
    # ..............................................................................................................
    # (A) FOR EACH ITERATION t = 1,...,T0, CONSTRUCT A CLASSIFICATION SEQUENCE RHO_t OF SIZE G_PLUS_HAT 
    # ..............................................................................................................
    if (G_plus_hat == 1) {
      
      class_seqs_net = matrix(1, nrow = T_iters_net, ncol = 1)
      
      posterior_means_modes$C[sim, ] = rep(1, M)
      
    } else { # G_plus_hat > 1
      
      kmeans_res_cluster_temp = kmeans_res_net$cluster
      class_seqs_net = vector("list", length = T_iters_net)
      
      for (t in 1:T_iters_net) {
        
        class_seqs_net[[t]] = kmeans_res_cluster_temp[1:dim(store_Zs_active[[t]])[3]]
        kmeans_res_cluster_temp = kmeans_res_cluster_temp[-(1:dim(store_Zs_active[[t]])[3])]
        
      } # end t for loop
      
    } # end if else statement
    
    # # ..............................................................................................................
    # # (B) FOR EACH ITERATION t = 1,...,T0, CARRY OUT A PERMUTATION TEST, CHECKING IF RHO_t IS A PERMUTATION OF 
    # #     (1, ..., G_PLUS_HAT); IF NOT, REMOVE DRAW t FROM CONSIDERATION IN THE REST OF THE ALGORITHM
    # # ..............................................................................................................
    # perm_indices_net = c()
    # 
    # for (t in 1:T_iters_net) {
    #   
    #   if (length(class_seqs_net[[t]]) == G_plus_hat) {
    #     
    #     if (all(sort(class_seqs_net[[t]]) == 1:G_plus_hat) == T) {
    #       perm_indices_net = c(perm_indices_net, t)
    #     } # end if statement
    #     
    #   } # end if else statement
    #   
    # } # end t for loop
    
    # ..............................................................................................................
    # (B) ALTERNATIVE VERSION (APRIL 4TH, 2025)
    # ..............................................................................................................
    if (G_plus_hat == 1) {
      # N/A
    } else {
      
      perm_indices_net = c()
      
      for (t in 1:T_iters_net) {
        
        if (length(sort(unique(class_seqs_net[[t]]))) == G_plus_hat) {
          
          if (all(sort(unique(class_seqs_net[[t]])) == 1:G_plus_hat)) {
            perm_indices_net = c(perm_indices_net, t)
          } else {
            # doesn't pass
          } # end if else statement
          # doesn't pass
        } else {
          # doesn't pass
        } # end if else statement
        
      } # end t for loop
      
    } # end if else statement
    
    # ..............................................................................................................
    # (C) WE KEEP ONLY THOSE DRAWS/ITERATIONS t THAT PASSED THE PERMUTATION TEST; IF THERE ARE NO DRAWS t THAT 
    #     PASSED THE PERMUTATION TEST, THEN WE ARE UNABLE TO MOVE FORWARD WITH THE POST-PRCOESSING ALGORITHM
    #     NOTE: THE NUMBER OF CLASSIFICATION SEQUENCES OF T0 NOT BEING A PERMUTATION IS DENOTED T0_RHO
    # ..............................................................................................................
    if (length(perm_indices_net) == 0) { # check if there are no permutations (then we cannot continue)
      
      mes = "There are no class. seqs that are perm.s of 1:%s. Unable to proceed with Sim. %s."
      message(sprintf(mes, G_plus_hat, sim))
      next
      
    } else if (length(perm_indices_net) == 1) {
      
      mes = "There is only one class. seq that is a perm. of 1:%s. Unable to proceed with Sim. %s."
      message(sprintf(mes, G_plus_hat, sim))
      next
      
    } # end if else statement
    
    params$perm_indices_net[[sim]] = perm_indices_net
    
    class_seqs_keep_net = class_seqs_net[perm_indices_net]
    
    T0_rho_net = T_iters_net - length(perm_indices_net)
    params$T0_rho_net[sim] = T0_rho_net
    
    # ..............................................................................................................
    # - SUBSET THE APPROPRIATE DRAWS/ITERATIONS FOR EACH PARAMETER
    # ..............................................................................................................
    store_alpha = output[[sim]]$res_LaPCoM$store_alpha # no further post-processing steps for the intercept, alpha
    store_alpha = store_alpha[perm_indices_net]
    chains$alpha[[sim]] = store_alpha
    posterior_means_modes$alpha[sim] = mean(chains$alpha[[sim]])
    
    store_C_vector = store_C_vector[, perm_indices_net]
    chains$C_vector[[sim]] = store_C_vector
    
    store_G_plus = output[[sim]]$res_LaPCoM$store_G_plus
    chains$G_plus[[sim]] = store_G_plus
    
    store_G = output[[sim]]$res_LaPCoM$store_G
    store_G = store_G[perm_indices_net]
    chains$G[[sim]] = store_G
    posterior_means_modes$G_hat[[sim]] = DescTools::Mode(chains$G[[sim]])
    
    store_e = output[[sim]]$res_LaPCoM$store_e
    store_e = store_e[perm_indices_net]
    chains$e[[sim]] = store_e
    posterior_means_modes$e[sim] = mean(chains$e[[sim]])
    
    store_tau_active = store_tau_active[perm_indices_net]
    store_Zs_active = store_Zs_active[perm_indices_net]
    store_mus_active = store_mus_active[perm_indices_net]
    store_Sigmas_active = store_Sigmas_active[perm_indices_net]
    store_S_vector_active = store_S_vector_active[perm_indices_net]
    store_Kg_active = store_Kg_active[perm_indices_net]
    store_Kg_plus_active = store_Kg_plus_active[perm_indices_net]
    store_pi_active = store_pi_active[perm_indices_net]
    store_w_active = store_w_active[perm_indices_net]
    
    # 
    
    acc_rate_alpha = output[[sim]]$res_LaPCoM$store_ARs$acc_rate_alpha
    acc_rate_alpha = acc_rate_alpha[perm_indices_net]
    acc_rates$alpha[sim] = mean(acc_rate_alpha)
    
    acc_rate_Z_active = acc_rate_Z_active[perm_indices_net]
    
    acc_rate_e = output[[sim]]$res_LaPCoM$store_ARs$acc_rate_e
    acc_rate_e = acc_rate_e[perm_indices_net]
    acc_rates$e[sim] = mean(acc_rate_e)
    
    acc_rate_w_active = acc_rate_w_active[perm_indices_net]
    
    # ..............................................................................................................
    # - CREATE RELEVANT PLOTS FOR CERTAIN PARAMETERS
    # ..............................................................................................................
    png(paste0(plots_path, "/Sim", sim, "_Alpha.png"))
    plot(chains$alpha[[sim]], type = "l", las = 1, xlab = "Iteration", ylab = expression(alpha))
    abline(h = posterior_means_modes$alpha[[sim]], col = "magenta", lty = 2, lwd = 2)
    legend("bottomright", legend = "Mean", lty = 2, lwd = 2, col = "magenta", bty = "n")
    dev.off()
    
    png(paste0(plots_path, "/Sim", sim, "_G+.png"))
    barplot(table(chains$G_plus[[sim]]), xlab = latex2exp::TeX("G"))
    dev.off()
    
    png(paste0(plots_path, "/Sim", sim, "_G.png"))
    barplot(table(chains$G[[sim]]), xlab = latex2exp::TeX("G"))
    dev.off()
    
    png(paste0(plots_path, "/Sim", sim, "_e.png"))
    plot(chains$e[[sim]], type = "l", las = 1, xlab = "Iteration", ylab = expression("e"[0]))
    abline(h = posterior_means_modes$e[[sim]], col = "magenta", lty = 2, lwd = 2)
    legend("bottomright", legend = "Mean", lty = 2, lwd = 2, col = "magenta", bty = "n")
    dev.off()
    
    # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    # 5. FOR THE REMAINING T0(1 - T_RHO) DRAWS, A UNIQUE LABELLING IS ACHIEVED BY RE-ORDERING THE DRAWS WE ARE KEEPING 
    #    ACCORDING TO THE CLASSIFICATION SEQUENCES RHO_t (FOR EVERY PARAMETER)
    #    NOTE: THE RE-ORDERED, IDENTIFIED PARAMETERS CAN BE USED FOR FURTHER COMPONENT-SPECIFIC PARAMETER INFERENCE
    # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    store_C_vector_relabelled = matrix(NA, nrow = M, ncol = length(perm_indices_net))
    for (t in 1:length(perm_indices_net)) {
      store_C_vector_relabelled[, t] = class_seqs_keep_net[[t]][store_C_vector[, t]]
    } # end t for loop
    
    store_taus_relabelled = matrix(NA, nrow = G_plus_hat, ncol = length(perm_indices_net))
    for (t in 1:length(perm_indices_net)) {
      for (g in 1:G_plus_hat) {
        store_taus_relabelled[g, t] = store_tau_active[[t]][which(class_seqs_keep_net[[t]] == g)] 
      } # end g for loop
    } # end t for loop
    chains$tau[[sim]] = store_taus_relabelled
    posterior_means_modes$tau[[sim]] = rowMeans(chains$tau[[sim]])
    
    store_Zs_relabelled = array(NA, dim = c(N, 2, G_plus_hat, length(perm_indices_net)))
    for (t in 1:length(perm_indices_net)) {
      if (G_plus_hat == 1) {
        store_Zs_relabelled[, , g, t] = store_Zs_active[[t]]
      } else {
        for (g in 1:G_plus_hat) {
          store_Zs_relabelled[, , g, t] = store_Zs_active[[t]][, , which(class_seqs_keep_net[[t]] == g)] 
        } # end g for loop
      } # end if else statement
    } # end t for loop
    chains$Z[[sim]] = store_Zs_relabelled
    
    store_mus_relabelled = array(NA, dim = c(K_max, 2, G_plus_hat, length(perm_indices_net)))
    for (t in 1:length(perm_indices_net)) {
      if (G_plus_hat == 1) {
        store_mus_relabelled[, , g, t] = store_mus_active[[t]]
      } else {
        for (g in 1:G_plus_hat) {
          store_mus_relabelled[, , g, t] = store_mus_active[[t]][, , which(class_seqs_keep_net[[t]] == g)] 
        } # end g for loop
      } # end if else statement
    } # end t for loop
    
    store_Sigmas_relabelled = array(NA, dim = c(2, 2, K_max, G_plus_hat, length(perm_indices_net)))
    for (t in 1:length(perm_indices_net)) {
      if (G_plus_hat == 1) {
        store_Sigmas_relabelled[, , , g, t] = store_Sigmas_active[[t]]
      } else {
        for (g in 1:G_plus_hat) {
          store_Sigmas_relabelled[, , , g, t] = store_Sigmas_active[[t]][, , , which(class_seqs_keep_net[[t]] == g)] 
        } # end g for loop
      } # end if else statement
    } # end t for loop
    
    store_S_vector_relabelled = array(NA, dim = c(N, G_plus_hat, length(perm_indices_net)))
    for (t in 1:length(perm_indices_net)) {
      if (G_plus_hat == 1) {
        store_S_vector_relabelled[, g, t] = store_S_vector_active[[t]]
      } else {
        for (g in 1:G_plus_hat) {
          store_S_vector_relabelled[, g, t] = store_S_vector_active[[t]][, which(class_seqs_keep_net[[t]] == g)] 
        } # end g for loop
      } # end if else statement
    } # end t for loop
    
    store_Kg_relabelled = matrix(NA, nrow = G_plus_hat, ncol = length(perm_indices_net))
    for (t in 1:length(perm_indices_net)) {
      if (G_plus_hat == 1) {
        store_Kg_relabelled[g, t] = store_Kg_active[[t]]
      } else {
        for (g in 1:G_plus_hat) {
          store_Kg_relabelled[g, t] = store_Kg_active[[t]][which(class_seqs_keep_net[[t]] == g)] 
        } # end g for loop
      } # end if else statement
    } # end t for loop
    
    store_Kg_plus_relabelled = matrix(NA, nrow = G_plus_hat, ncol = length(perm_indices_net))
    for (t in 1:length(perm_indices_net)) {
      if (G_plus_hat == 1) {
        store_Kg_plus_relabelled[g, t] = store_Kg_plus_active[[t]]
      } else {
        for (g in 1:G_plus_hat) {
          store_Kg_plus_relabelled[g, t] = store_Kg_plus_active[[t]][which(class_seqs_keep_net[[t]] == g)] 
        } # end g for loop
      } # end if else statement
    } # end t for loop
    
    store_pi_relabelled = array(NA, dim = c(K_max, G_plus_hat, length(perm_indices_net)))
    for (t in 1:length(perm_indices_net)) {
      if (G_plus_hat == 1) {
        store_pi_relabelled[, g, t] = store_pi_active[[t]]
      } else {
        for (g in 1:G_plus_hat) {
          store_pi_relabelled[, g, t] = store_pi_active[[t]][, which(class_seqs_keep_net[[t]] == g)] 
        } # end g for loop
      } # end if else statement
    } # end t for loop
    
    store_w_relabelled = matrix(NA, nrow = G_plus_hat, ncol = length(perm_indices_net))
    for (t in 1:length(perm_indices_net)) {
      if (G_plus_hat == 1) {
        store_w_relabelled[g, t] = store_w_active[[t]]
      } else {
        for (g in 1:G_plus_hat) {
          store_w_relabelled[g, t] = store_w_active[[t]][which(class_seqs_keep_net[[t]] == g)] 
        } # end g for loop
      } # end if else statement
    } # end t for loop
    chains$w[[sim]] = store_w_relabelled
    posterior_means_modes$w[[sim]] = rowMeans(chains$w[[sim]])
    
    acc_rate_Z_relabelled = matrix(NA, nrow = G_plus_hat, ncol = length(perm_indices_net))
    for (t in 1:length(perm_indices_net)) {
      for (g in 1:G_plus_hat) {
        acc_rate_Z_relabelled[g, t] = acc_rate_Z_active[[t]][which(class_seqs_keep_net[[t]] == g)] 
      } # end g for loop
    } # end t for loop
    acc_rates$Z[[sim]] = rowMeans(acc_rate_Z_relabelled)
    
    acc_rate_w_relabelled = matrix(NA, nrow = G_plus_hat, ncol = length(perm_indices_net))
    for (t in 1:length(perm_indices_net)) {
      for (g in 1:G_plus_hat) {
        acc_rate_w_relabelled[g, t] = acc_rate_w_active[[t]][which(class_seqs_keep_net[[t]] == g)] 
      } # end g for loop
    } # end t for loop
    acc_rates$w[[sim]] = rowMeans(acc_rate_w_relabelled)
    
    # ................................................................................................................
    # - CREATE RELEVANT PLOTS FOR CERTAIN PARAMETERS
    # ................................................................................................................
    for (g in 1:G_plus_hat) {
      
      png(paste0(plots_path, "/Sim", sim, "_w", g, ".png"))
      plot(chains$w[[sim]][g, ], type = "l", las = 1, xlab = "Iteration", ylab = expression(omega[0]))
      abline(h = posterior_means_modes$w[[sim]][g], col = "magenta", lty = 2, lwd = 2)
      legend("bottomright", legend = "Mean", lty = 2, lwd = 2, col = "magenta", bty = "n")
      dev.off()
      
    } # end g for loop
    
    # ................................................................................................................
    # - APPLY A RELABELLING FUNCTION TO THE POSTERIOR OPTIMAL CLUSTERING SOLUTION SO THAT THERE IS A CORRESPONDENCE BW
    #   IT AND THE OTHER PARAMETERS FOR LATER ANALYSIS
    # ................................................................................................................
    WG_order_C = t(apply(store_C_vector_relabelled, 2, 
                         function(col) { as.vector(e1071::matchClasses(table(posterior_means_modes$C[sim, ], col), verbose = F)) }))
    
    stacked_WG_C = apply(WG_order_C, 1, paste, collapse = ",")
    stacked_WG_C_freq = table(stacked_WG_C)
    g_match = as.integer(unlist(stringr::str_split(names(stacked_WG_C_freq)[which.max(stacked_WG_C_freq)], ",")))
    
    posterior_means_modes$C[sim, ] = g_match[posterior_means_modes$C[sim, ]]
    posterior_means_modes$G_plus_hat_FS[sim] = length(unique(posterior_means_modes$C[sim, ]))
    
    # ................................................................................................................
    # - OBTAIN THE POSTERIOR MEAN LATENT SPACES
    # ................................................................................................................
    proc_matched_Z = array(NA, dim = dim(store_Zs_relabelled))
    
    store_mus_relabelled_rotated = array(NA, dim = dim(store_mus_relabelled))
    store_Sigmas_relabelled_rotated = array(NA, dim = dim(store_Sigmas_relabelled))
    
    ref_index = 1
    
    for (g in 1:G_plus_hat) {
      
      proc_matched_Z[, , g, 1] = store_Zs_relabelled[, , g, ref_index]
      store_mus_relabelled_rotated[, , g, 1] = store_mus_relabelled[, , g, 1]
      store_Sigmas_relabelled_rotated[, , , g, 1] = store_Sigmas_relabelled[, , , g, 1]
      
      for (t in 2:dim(proc_matched_Z)[4]) {
        
        proc_Z_fit = vegan::procrustes(proc_matched_Z[, , g, 1], store_Zs_relabelled[, , g, t])
        # proc_matched_Z[, , g, t] = proc_Z_fit$Yrot # pre November 4th 2024
        proc_matched_Z[, , g, t] = sweep(((proc_Z_fit$scale * sweep(store_Zs_relabelled[, , g, t], 2, proc_Z_fit$xmean, "-")) %*% proc_Z_fit$rotation), 2, proc_Z_fit$translation, "+")
        
        # also rotate/translate/scale the mean vectors mu
        mu_scale_rot_trans = sweep((proc_Z_fit$scale * (sweep(store_mus_relabelled[, , g, t], 2, proc_Z_fit$xmean, "-") %*% proc_Z_fit$rotation)), 2, proc_Z_fit$translation, "+")
        store_mus_relabelled_rotated[, , g, t] = mu_scale_rot_trans
        
        for (k in 1:store_Kg_plus_relabelled[g, t]) { # also scale the covariance matrices Sigma
          
          Sigma_scale_rot = (proc_Z_fit$scale * store_Sigmas_relabelled[, , k, g, t])
          store_Sigmas_relabelled_rotated[, , k, g, t] = Sigma_scale_rot
          
        } # end k for loop
        
      } # end t for loop
      
    } # end g for loop
    
    chains$proc_matched_Z[[sim]] = proc_matched_Z
    
    posterior_means_modes$Z[[sim]] = vector("list", length = G_plus_hat)
    
    for (g in 1:G_plus_hat) {
      posterior_means_modes$Z[[sim]][[g]] = apply(proc_matched_Z[, , g, ], 1:2, mean)
    } # end g for loop
    
    # ------------------------------------------------------------------------------------------------------------------
    # NODE-LEVEL PARAMETERS
    # ------------------------------------------------------------------------------------------------------------------
    params$perm_indices_node[[sim]] = vector("list", length = G_plus_hat) # classifications/iterations that were permutations
    params$T0_rho_node[[sim]] = rep(NA, G_plus_hat) # number of class. seq.s / iterations that were not permutations
    
    chains$mu[[sim]] = vector("list", length = G_plus_hat)
    chains$Sigma[[sim]] = vector("list", length = G_plus_hat)
    chains$S_vector[[sim]] = vector("list", length = G_plus_hat)
    chains$Kg[[sim]] = vector("list", length = G_plus_hat)
    chains$Kg_plus[[sim]] = vector("list", length = G_plus_hat)
    chains$pi[[sim]] = vector("list", length = G_plus_hat)
    
    posterior_means_modes$mu[[sim]] = vector("list", length = G_plus_hat)
    posterior_means_modes$Sigma[[sim]] = vector("list", length = G_plus_hat)
    posterior_means_modes$S[[sim]] = vector("list", length = G_plus_hat)
    posterior_means_modes$Kg_hat[[sim]] = rep(NA, G_plus_hat)
    posterior_means_modes$Kg_plus_hat_WG[[sim]] = rep(NA, G_plus_hat)
    posterior_means_modes$Kg_plus_hat_FS[[sim]] = rep(NA, G_plus_hat)
    posterior_means_modes$pi[[sim]] = vector("list", length = G_plus_hat)
    
    for (g in 1:G_plus_hat) {
      
      T_iters_node = dim(store_S_vector_relabelled)[3]
      
      # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
      # 1. USE THE WADE AND GHARAMANI METHOD TO OBTAIN THE OPTIMAL NETWORK-LEVEL CLUSTERING (USING ALL POSTERIOR SAMPLES)
      # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
      psm_node = mcclust::comp.psm(t(store_S_vector_relabelled[, g, ]))
      wade_ghar_node = mcclust.ext::minVI(psm_node, method = "all", cls.draw = t(store_S_vector_relabelled[, g, ]), 
                                          include.greedy = F)
      
      posterior_means_modes$S[[sim]][[g]] = wade_ghar_node$cl[1, ]
      
      # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
      # 2. DETERMINE THE NUMBER OF NON-EMPTY COMPONENTS
      # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
      Kg_plus_hat = params$Kg_plus_hat[[sim]][g] = posterior_means_modes$Kg_plus_hat_WG[[sim]][g] = length(unique(posterior_means_modes$S[[sim]][[g]]))
      
      # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
      # 3. RESULTS IN A CLASSIFICATION INDICATING TO WHICH CLUSTERS THE COMPONENT-SPECIFIC PARAMS OF EACH DRAW BELONG
      # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
      
      # ..............................................................................................................
      # (A) ISOLATE THE NON-EMPTY COMPONENTS
      # ..............................................................................................................
      non_empty_node_clusts = vector("list", length = T_iters_node)
      
      for (t in 1:T_iters_node) {
        non_empty_node_clusts[[t]] = sort(unique(store_S_vector_relabelled[, g, t]))
      } # end t for loop
      
      store_mus_node_level_g_active = vector("list", length = T_iters_node)
      store_Sigmas_node_level_g_active = vector("list", length = T_iters_node)
      store_pi_node_level_g_active = vector("list", length = T_iters_node)
      
      for (t in 1:T_iters_node) { 
        
        store_mus_node_level_g_active[[t]] = store_mus_relabelled_rotated[non_empty_node_clusts[[t]], , g, t]
        store_Sigmas_node_level_g_active[[t]] = store_Sigmas_relabelled_rotated[, , non_empty_node_clusts[[t]], g, t]
        store_pi_node_level_g_active[[t]] = store_pi_relabelled[non_empty_node_clusts[[t]], g, t]
        
      } # end t for loop
      
      # ..............................................................................................................
      # (B) N/A AT THE NODE-LEVEL
      # ..............................................................................................................
      
      # ..............................................................................................................
      # (C) ARRANGE THE CLUSTER MEANS INTO A MATRIX WITH K_g_PLUS_HAT X T ROWS AND p COLUMNS
      # ..............................................................................................................
      mu_means = matrix(NA, nrow = 0, ncol = 2)
      for (t in 1:T_iters_node) { mu_means = rbind(mu_means, store_mus_node_level_g_active[[t]]) } # end t for loop
      
      # ..............................................................................................................
      # (D) CLUSTER THE MATRIX INTO G_PLUS_HAT CLUSTERS USING K-MEANS CLUSTERING
      # ..............................................................................................................
      kmeans_res_node = kmeans(scale(mu_means), centers = Kg_plus_hat, nstart = 100, iter.max = 100)
      
      # ..............................................................................................................
      # - CREATE RELEVANT PLOTS FOR CERTIAN PARAMETERS
      # ..............................................................................................................
      png(paste0(plots_path, "/Sim", sim, "_Node-Level_PPR_LS", g, ".png"))
      plot(mu_means, xlab = "Dim. 1", ylab = "Dim. 2", las = 1, pch = 19, col = kmeans_res_node$cluster,
           main = paste("Point Process Representation at the Node-Level of LS", g))
      dev.off()
      
      # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
      # 4.
      # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
      
      # ..............................................................................................................
      # (A) FOR EACH ITERATION t = 1,...,T0, CONSTRUCT A CLASSIFICATION SEQUENCE RHO_t OF SIZE G_PLUS_HAT
      # ..............................................................................................................
      if (Kg_plus_hat == 1) {
        
        class_seqs_node = matrix(1, nrow = T_iters_node, ncol = 1)
        params$perm_indices_node[[sim]][[g]] = 1:T_iters_node
        T0_rho_node = 0
        params$T0_rho_node[[sim]][g] = T0_rho_node
        
        # 
        
        store_mus_node_level_g_active_temp = matrix(NA, nrow = 0, ncol = 2)
        for (t in 1:T_iters_node) {
          store_mus_node_level_g_active_temp = rbind(store_mus_node_level_g_active_temp, store_mus_node_level_g_active[[t]])
        } # end t for loop
        store_mus_node_level_g_active = store_mus_node_level_g_active_temp; rm(store_mus_node_level_g_active_temp)
        
        store_Sigmas_node_level_g_active_temp = matrix(NA, nrow = 0, ncol = 2)
        for (t in 1:T_iters_node) {
          if (length(dim(store_Sigmas_node_level_g_active[[t]])) == 2) { # only a matrix, not an array
            store_Sigmas_node_level_g_active_temp = rbind(store_Sigmas_node_level_g_active_temp, 
                                                          diag(store_Sigmas_node_level_g_active[[t]]))
          } else {
            store_Sigmas_node_level_g_active_temp = rbind(store_Sigmas_node_level_g_active_temp, 
                                                          t(apply(store_Sigmas_node_level_g_active[[t]], 3, diag)))
          } # end if else statement
        } # end t for loop
        store_Sigmas_node_level_g_active = store_Sigmas_node_level_g_active_temp; rm(store_Sigmas_node_level_g_active_temp)
        
        store_S_vector_node_level_g = store_S_vector_relabelled[, g, ]
        
        posterior_means_modes$S[[sim]][[g]] = rep(1, N)
        
      } else { # Kg_plus_hat > 1
        
        kmeans_res_cluster_temp = kmeans_res_node$cluster
        class_seqs_node = vector("list", length = T_iters_node)
        
        for (t in 1:T_iters_node) {
          
          class_seqs_node[[t]] = kmeans_res_cluster_temp[1:nrow(matrix(store_mus_node_level_g_active[[t]], ncol = 2))]
          kmeans_res_cluster_temp = kmeans_res_cluster_temp[-(1:nrow(matrix(store_mus_node_level_g_active[[t]], ncol = 2)))]
          
        } # end t for loop

      } # end if else statement
      
      # ..............................................................................................................
      # (B) ALTERNATIVE VERSION (APRIL 4TH, 2025)
      #     FOR EACH ITERATION t = 1,...,T0, CARRY OUT A PERMUTATION TEST, CHECKING IF RHO_t IS A PERMUTATION OF
      #     (1, ..., G_PLUS_HAT); IF NOT, REMOVE DRAW t FROM CONSIDERATION IN THE REST OF THE ALGORITHM
      # ..............................................................................................................
      if (Kg_plus_hat == 1) {
        # N/A
      } else {
        
        perm_indices_node = c()
        
        for (t in 1:T_iters_node) {
          
          if (length(sort(unique(class_seqs_node[[t]]))) == Kg_plus_hat) {
            
            if (all(sort(unique(class_seqs_node[[t]])) == 1:Kg_plus_hat)) {
              perm_indices_node = c(perm_indices_node, t)
            } else {
              # doesn't pass
            } # end if else statement
            # doesn't pass
          } else {
            # doesn't pass
          } # end if else statement
          
        } # end t for loop
        
      } # end if else statement
      
      # ..............................................................................................................
      # (C) WE KEEP ONLY THOSE DRAWS/ITERATIONS t THAT PASSED THE PERMUTATION TEST; IF THERE ARE NO DRAWS t THAT
      #     PASSED THE PERMUTATION TEST, THEN WE ARE UNABLE TO MOVE FORWARD WITH THE POST-PRCOESSING ALGORITHM
      #     NOTE: THE NUMBER OF CLASSIFICATION SEQUENCES OF T0 NOT BEING A PERMUTATION IS DENOTED T0_RHO
      # ..............................................................................................................
      if (Kg_plus_hat == 1) {
        # N/A
      } else {
        
        # check if there are no permutations (then we cannot continue with the algorithm)
        if (length(perm_indices_node) == 0) {
          
          mes = "There are no class. seqs that are perm.s of 1:%s in Z_%s. Unable to proceed with Sim. %s."
          message(sprintf(mes, Kg_plus_hat, g, sim))
          next
          
        } else if (length(perm_indices_node) == 1) {
          
          mes = "There is only one class. seq that is perm. of 1:%s in Z_%s. Unable to proceed with Sim. %s."
          message(sprintf(mes, Kg_plus_hat, g, sim))
          next
          
        } # end if else statement
        
        params$perm_indices_node[[sim]][[g]] = perm_indices_node
        
        class_seqs_keep_node = class_seqs_node[perm_indices_node]
        
        T0_rho_node = T_iters_node - length(perm_indices_node)
        params$T0_rho_node[[sim]][g] = T0_rho_node
        
      } # end if else statement
      
      # ..............................................................................................................
      # - SUBSET THE APPROPRIATE DRAWS/ITERATIONS FOR EACH PARAMETER
      # ..............................................................................................................
      if (Kg_plus_hat == 1) {
        # N/A
      } else {
        
        store_mus_node_level_g_active = store_mus_node_level_g_active[perm_indices_node]
        store_Sigmas_node_level_g_active = store_Sigmas_node_level_g_active[perm_indices_node]
        store_pi_node_level_g_active = store_pi_node_level_g_active[perm_indices_node]
        
        store_S_vector_node_level_g = store_S_vector_relabelled[, g, perm_indices_node]
        
      } # end if else statement
      
      # ..............................................................................................................
      # - CREATE RELEVANT PLOTS FOR CERTAIN PARAMETERS
      # ..............................................................................................................
      if (Kg_plus_hat == 1) {
        
        chains$Kg_plus[[sim]] = store_Kg_plus_relabelled[g, ]
        png(paste0(plots_path, "/Sim", sim, "_LS", g, "_Kg+.png"))
        barplot(table(chains$Kg_plus[[sim]]))
        dev.off()
        
        chains$Kg[[sim]] = store_Kg_relabelled[g, ]
        posterior_means_modes$Kg_hat[[sim]][g] = DescTools::Mode(chains$Kg[[sim]])[1] # if more than one, take smaller
        png(paste0(plots_path, "/Sim", sim, "_LS", g, "_Kg.png"))
        barplot(table(chains$Kg[[sim]]))
        dev.off()
        
      } else {
        
        store_Kg_plus_level_g = store_Kg_plus_relabelled[g, perm_indices_node]
        chains$Kg_plus[[sim]] = store_Kg_plus_level_g
        png(paste0(plots_path, "/Sim", sim, "_LS", g, "_Kg+.png"))
        barplot(table(chains$Kg_plus[[sim]]))
        dev.off()
        
        # 
        
        store_Kg_level_g = store_Kg_relabelled[g, perm_indices_node]
        chains$Kg[[sim]] = store_Kg_level_g
        posterior_means_modes$Kg_hat[[sim]][g] = DescTools::Mode(chains$Kg[[sim]])[1]
        png(paste0(plots_path, "/Sim", sim, "_LS", g, "_Kg.png"))
        barplot(table(chains$Kg[[sim]]))
        dev.off()
        
      } # end if else statement
      
      # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
      # 5. FOR THE REMAINING T0(1 - T_RHO) DRAWS, A UNIQUE LABELLING IS ACHIEVED BY RE-ORDERING THE DRAWS WE ARE KEEPING
      #    ACCORDING TO THE CLASSIFICATION SEQUENCES RHO_t (FOR EVERY PARAMETER)
      #    NOTE: THE RE-ORDERED, IDENTIFIED PARAMETERS CAN BE USED FOR FURTHER COMPONENT-SPECIFIC PARAMETER INFERENCE
      # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
      if (Kg_plus_hat == 1) {
        
        store_mus_node_level_g_relabelled = store_mus_node_level_g_active
        chains$mu[[sim]][[g]] = store_mus_node_level_g_relabelled
        
        store_Sigmas_node_level_g_relabelled = store_Sigmas_node_level_g_active
        chains$Sigma[[sim]][[g]] = store_Sigmas_node_level_g_relabelled
        
        chains$pi[[sim]][[g]] = store_pi_node_level_g_active
        posterior_means_modes$pi[[sim]][[g]] = 1
        
        chains$S_vector[[sim]][[g]] = store_S_vector_node_level_g_relabelled = store_S_vector_node_level_g
        
      } else {
        
        store_mus_node_level_g_relabelled = array(NA, dim = c(dim(store_mus_node_level_g_active[[1]])[1],
                                                              dim(store_mus_node_level_g_active[[1]])[2],
                                                              length(store_mus_node_level_g_active)))
        for (t in 1:length(perm_indices_node)) {
          for (k in 1:Kg_plus_hat) {
            inds = which(class_seqs_keep_node[[t]] == k)
            if (length(inds) == 1) { # only one component assigned to this cluster
              store_mus_node_level_g_relabelled[k, , t] = store_mus_node_level_g_active[[t]][inds, ]
            } else { # more than one component assigned to this cluster: merge them
              store_mus_node_level_g_relabelled[k, , t] = colMeans(store_mus_node_level_g_active[[t]][inds, ])
            } # end if else statement
          } # end for k loop
        } # end t for loop
        chains$mu[[sim]][[g]] = store_mus_node_level_g_relabelled
        
        store_Sigmas_node_level_g_relabelled = array(NA, dim = c(dim(store_Sigmas_node_level_g_active[[1]])[1],
                                                                 dim(store_Sigmas_node_level_g_active[[1]])[2],
                                                                 dim(store_Sigmas_node_level_g_active[[1]])[3],
                                                                 length(store_Sigmas_node_level_g_active)))
        for (t in 1:length(perm_indices_node)) {
          for (k in 1:Kg_plus_hat) {
            inds = which(class_seqs_keep_node[[t]] == k)
            if (length(inds) == 1) { # only one component assigned to this cluster
              store_Sigmas_node_level_g_relabelled[, , k, t] = store_Sigmas_node_level_g_active[[t]][, , inds]
            } else { # more than one component assigned to this cluster: merge them
              store_Sigmas_node_level_g_relabelled[, , k, t] = apply(store_Sigmas_node_level_g_active[[t]][, , inds], 1:2, mean)
            } # end if else statement
          } # end for k loop
        } # end t for loop
        chains$Sigma[[sim]][[g]] = store_Sigmas_node_level_g_relabelled

        store_pi_node_level_g_relabelled = matrix(NA, nrow = length(store_pi_node_level_g_active), ncol = Kg_plus_hat)
        for (t in 1:length(perm_indices_node)) {
          for (k in 1:Kg_plus_hat) {
            inds = which(class_seqs_keep_node[[t]] == k)
            if (length(inds) == 1) { # only one component assigned to this cluster
              store_pi_node_level_g_relabelled[t, k] = store_pi_node_level_g_active[[t]][inds]
            } else { # more than one component assigned to this cluster: merge them
              store_pi_node_level_g_relabelled[t, k] = mean(store_pi_node_level_g_active[[t]][inds])
            } # end if else statement
          } # end for k loop
        } # end t for loop
        chains$pi[[sim]][[g]] = store_pi_node_level_g_relabelled
        posterior_means_modes$pi[[sim]][[g]] = colMeans(chains$pi[[sim]][[g]])

        store_S_vector_node_level_g_relabelled = matrix(NA, nrow = nrow(store_S_vector_node_level_g), ncol = ncol(store_S_vector_node_level_g))
        for (t in 1:length(perm_indices_node)) {
          store_S_vector_node_level_g_relabelled[, t] = class_seqs_keep_node[[t]][store_S_vector_node_level_g[, t]]
        } # end t for loop
        chains$S_vector[[sim]][[g]] = store_S_vector_node_level_g_relabelled
        
      } # end if else statement
      
      # ................................................................................................................
      # - APPLY A RELABELLING FUNCTION TO THE POSTERIOR OPTIMAL CLUSTERING SOLUTION SO THAT THERE IS A CORRESPONDENCE BW
      #   IT AND THE OTHER PARAMETERS FOR LATER ANALYSIS
      # ................................................................................................................
      store_S_vector_node_level_g_relabelled = matrix(NA, nrow = nrow(store_S_vector_relabelled[, g, ]), 
                                                      ncol = ncol(store_S_vector_relabelled[, g, ]))
      for (t in 1:dim(store_S_vector_relabelled)[3]) {
        store_S_vector_node_level_g_relabelled[, t] = class_seqs_node[[t]][store_S_vector_relabelled[, g, t]]
      } # end t for loop
      chains$S_vector[[sim]][[g]] = store_S_vector_node_level_g_relabelled
      
      WG_order_S = t(apply(store_S_vector_node_level_g_relabelled, 2, 
                           function(col) { as.vector(e1071::matchClasses(table(posterior_means_modes$S[[sim]][[g]], col), verbose = F)) }))
      
      stacked_WG_S = apply(WG_order_S, 1, paste, collapse = ",")
      stacked_WG_S_freq = table(stacked_WG_S)
      k_match = as.integer(unlist(stringr::str_split(names(stacked_WG_S_freq)[which.max(stacked_WG_S_freq)], ",")))
      
      posterior_means_modes$S[[sim]][[g]] = k_match[posterior_means_modes$S[[sim]][[g]]]
      posterior_means_modes$Kg_plus_hat_FS[[sim]][g] = length(unique(posterior_means_modes$S[[sim]][[g]]))
      
      # ................................................................................................................
      # - OBTAIN THE POSTERIOR MEAN CLUSTER MEANS AND COVARIANCE MATRICES
      # ................................................................................................................
      posterior_means_modes$mu[[sim]][[g]] = matrix(NA, nrow = Kg_plus_hat, ncol = 2)
      posterior_means_modes$Sigma[[sim]][[g]] = matrix(NA, nrow = Kg_plus_hat, ncol = 2)
      
      if (Kg_plus_hat == 1) {
        
        posterior_means_modes$mu[[sim]][[g]] = colMeans(store_mus_node_level_g_relabelled)
        posterior_means_modes$Sigma[[sim]][[g]] = colMeans(store_Sigmas_node_level_g_relabelled)
        
      } else {
        
        for (k in 1:Kg_plus_hat) {
          
          posterior_means_modes$mu[[sim]][[g]][k, ] = rowMeans(store_mus_node_level_g_relabelled[k, , ])
          posterior_means_modes$Sigma[[sim]][[g]][k, ] = rowMeans(apply(store_Sigmas_node_level_g_relabelled[, , k, ], 3, diag))
          
        } # end k for loop
        rm(k)
        
      } # end if else statement
      
    } # end g for loop
    rm(g)
    
    # ------------------------------------------------------------------------------------------------------------------
    # PLOTS
    # ------------------------------------------------------------------------------------------------------------------
    cols_update = rep(cols, 3)
    cols_g = vector("list", length = G_plus_hat)
    
    for (g in 1:G_plus_hat) {
      
      start = ifelse(g == 1, 1, sum(posterior_means_modes$Kg_plus_hat_WG[[sim]][1:(g - 1)]) + 1)
      end = sum(posterior_means_modes$Kg_plus_hat_WG[[sim]][1:g])
      
      cols_g[[g]] = cols_update[start:end]
      
    } # end g for loop
    rm(g)
    
    for (g in 1:G_plus_hat) {
      
      if (is.null(posterior_means_modes$mu[[sim]][[g]])) {
        
        next
        
      } else {
        
        png(paste0(plots_path, "/Sim", sim, "_Z", g, ".png"))
        
        plot(posterior_means_modes$Z[[sim]][[g]], xlab = expression("z"[1]), ylab = expression("z"[2]), pch = 19, las = 1,
             col = adjustcolor(cols_g[[g]][posterior_means_modes$S[[sim]][[g]]], alpha = 0.8),
             # xlim = c(-3, 3), ylim = c(-3, 3), 
             main = paste("Latent Space", g))
        grid()
        
        if (posterior_means_modes$Kg_plus_hat_WG[[sim]][g] == 1) {
          
          mean_point = posterior_means_modes$mu[[sim]][[g]]
          covariance = matrix(c(posterior_means_modes$Sigma[[sim]][[g]][1], 0, 0, posterior_means_modes$Sigma[[sim]][[g]][2]), 2, 2, byrow = T)
          eigen_decomp = eigen(covariance)
          eigenvalues = eigen_decomp$values
          eigenvectors = eigen_decomp$vectors
          axis_length = sqrt(eigenvalues)
          
          arrows(mean_point[1], mean_point[2], mean_point[1] + eigenvectors[1, 1] * axis_length[1],
                 mean_point[2] + eigenvectors[2, 1] * axis_length[1], col = "gray60", length = 0.1, lty = 2, angle = 0)
          
          arrows(mean_point[1], mean_point[2], mean_point[1] - eigenvectors[1, 1] * axis_length[1],
                 mean_point[2] - eigenvectors[2, 1] * axis_length[1], col = "gray60", length = 0.1, lty = 2, angle = 0)
          
          arrows(mean_point[1], mean_point[2], mean_point[1] + eigenvectors[1, 2] * axis_length[2],
                 mean_point[2] + eigenvectors[2, 2] * axis_length[2], col = "gray60", length = 0.1, lty = 2, angle = 0)
          
          arrows(mean_point[1], mean_point[2], mean_point[1] - eigenvectors[1, 2] * axis_length[2],
                 mean_point[2] - eigenvectors[2, 2] * axis_length[2], col = "gray60", length = 0.1, lty = 2, angle = 0)
          
          points(x = posterior_means_modes$mu[[sim]][[g]][1], y = posterior_means_modes$mu[[sim]][[g]][2],
                 pch = 15, col = cols_g[[g]][1], cex = 1.5)
          
          points(x = posterior_means_modes$mu[[sim]][[g]][1], y = posterior_means_modes$mu[[sim]][[g]][2],
                 pch = 7, col = "gray60", bg = cols_g[[g]][1], cex = 1.5)
          
          DescTools::DrawEllipse(x = posterior_means_modes$mu[[sim]][[g]][1], y = posterior_means_modes$mu[[sim]][[g]][2],
                                 radius.x = axis_length[1], radius.y = axis_length[2], border = "gray60", col = NA, lty = 2)
          
        } else {
          
          empirical_means_g = matrix(NA, nrow = posterior_means_modes$Kg_plus_hat_WG[[sim]][g], ncol = 2)
          closest_mean = rep(NA, posterior_means_modes$Kg_plus_hat_WG[[sim]][g])
          
          for (k in 1:posterior_means_modes$Kg_plus_hat_WG[[sim]][g]) {
            
            k_inds = which(posterior_means_modes$S[[sim]][[g]] == k)
            
            if (length(k_inds) == 0) {
              empirical_means_g[k, ] = NA 
            } else if (length(k_inds) == 1) {
              empirical_means_g[k, ] = posterior_means_modes$Z[[sim]][[g]][k_inds, ]
            } else {
              empirical_means_g[k, ] = colMeans(posterior_means_modes$Z[[sim]][[g]][k_inds, ])
            } # end if else statement
            
            dists = rep(NA, posterior_means_modes$Kg_plus_hat_WG[[sim]][g])
            
            for (kk in 1:posterior_means_modes$Kg_plus_hat_WG[[sim]][g]) {
              dists[kk] = sqrt(sum((empirical_means_g[k, ] - posterior_means_modes$mu[[sim]][[g]][kk, ])^2))
            } # end kk for loop
            
            if (all(is.na(dists))) {
              closest_mean[k] = NA
            } else {
              closest_mean[k] = which.min(dists)
            } # end if else statement
            
          } # end k for loop
          
          posterior_means_modes_mu_reordered = matrix(NA, nrow = posterior_means_modes$Kg_plus_hat_WG[[sim]][g], ncol = 2)
          posterior_means_modes_Sigmas_reordered = matrix(NA, nrow = posterior_means_modes$Kg_plus_hat_WG[[sim]][g], ncol = 2)
          
          for (k in 1:posterior_means_modes$Kg_plus_hat_WG[[sim]][g]) {
            
            if (is.na(closest_mean[k])) {
              
              posterior_means_modes_mu_reordered[k, ] = NA
              posterior_means_modes_Sigmas_reordered[k, ] = NA
              
            } else {
              
              posterior_means_modes_mu_reordered[k, ] = posterior_means_modes$mu[[sim]][[g]][closest_mean[k], ]
              posterior_means_modes_Sigmas_reordered[k, ] = posterior_means_modes$Sigma[[sim]][[g]][closest_mean[k], ]
              
            } # end if else statement
            
          } # end k for loop
          
          for (k in 1:posterior_means_modes$Kg_plus_hat_WG[[sim]][g]) {
            
            if (all(is.na(posterior_means_modes_mu_reordered[k, ]))) {
              
            } else {
              
              mean_point = posterior_means_modes_mu_reordered[k, ]
              covariance = matrix(c(posterior_means_modes_Sigmas_reordered[k, 1], 0, 0, posterior_means_modes_Sigmas_reordered[k, 2]), 2, 2, byrow = T)
              eigen_decomp = eigen(covariance)
              eigenvalues = eigen_decomp$values
              eigenvectors = eigen_decomp$vectors
              axis_length = sqrt(eigenvalues)
              
              arrows(mean_point[1], mean_point[2], mean_point[1] + eigenvectors[1, 1] * axis_length[1], 
                     mean_point[2] + eigenvectors[2, 1] * axis_length[1], col = "gray60", length = 0.1, lty = 2, angle = 0)
              
              arrows(mean_point[1], mean_point[2], mean_point[1] - eigenvectors[1, 1] * axis_length[1], 
                     mean_point[2] - eigenvectors[2, 1] * axis_length[1], col = "gray60", length = 0.1, lty = 2, angle = 0)
              
              arrows(mean_point[1], mean_point[2], mean_point[1] + eigenvectors[1, 2] * axis_length[2], 
                     mean_point[2] + eigenvectors[2, 2] * axis_length[2], col = "gray60", length = 0.1, lty = 2, angle = 0)
              
              arrows(mean_point[1], mean_point[2], mean_point[1] - eigenvectors[1, 2] * axis_length[2], 
                     mean_point[2] - eigenvectors[2, 2] * axis_length[2], col = "gray60", length = 0.1, lty = 2, angle = 0)
              
              points(x = posterior_means_modes_mu_reordered[k, 1], y = posterior_means_modes_mu_reordered[k, 2], 
                     pch = 15, col = cols_g[[g]][k], cex = 1.5) # draws WG
              
              points(x = posterior_means_modes_mu_reordered[k, 1], y = posterior_means_modes_mu_reordered[k, 2], 
                     pch = 7, col = "gray60", bg = cols_g[[g]][k], cex = 1.5)
              
              DescTools::DrawEllipse(x = posterior_means_modes_mu_reordered[k, 1], y = posterior_means_modes_mu_reordered[k, 2],
                                     radius.x = axis_length[1], radius.y = axis_length[2], border = "gray60", col = NA, lty = 2)
              
            } # end if else statement
            
          } # end k for loop
          
        } # end if else statement
        
        dev.off()
        
      } # end if else statement
      
    } # end g for loop
    
  } # end sim for loop
  
  return(list(params = params,
              chains = chains,
              acc_rates = acc_rates,
              posterior_means_modes = posterior_means_modes))
  
} # end post_process_LaPCoM function
