########################################################################################################################
# FUNCTION TO PERFORM POST-PROCESSING ON THE RESULTING MCMC CHAIN FROM monoLaPCM 
########################################################################################################################
post_process_monoLaPCM = function(output, # a list of the output of monoLaPCM algorithm for each simulation
                               multi, # the multiplex in array form with dimensions (N, N, M)
                               plots_path, # name of folder to store plots
                               cols = khroma::color("light")(9), seed = 123) {
  
    # ====================================================================================================================
    # SET THE SEED (BC WE USE KMEANS)
    # ====================================================================================================================
    set.seed(seed)
    
    # ====================================================================================================================
    # POST PROCESSING PROCEDURE FOR monoLaPCM (WITHOUT NODE CLUSTERING)
    # ====================================================================================================================
    M = dim(multi)[3]
    N = dim(multi)[1]
    
    G_max = dim(output[[1]]$res_monoLaPCM$store_Zs)[3]

    params = list(M = M, N = N, G_max = G_max)
    params$T0_iters_net = vector("list", length = length(output)) # iterations where G+ = G+hat
    params$perm_indices_net = vector("list", length = length(output)) # classifications/iterations that were permutations
    params$T0_rho_net = rep(NA, length(output)) # number of class. seq.s / iterations that were not permutations

    chains = list(alpha = vector("list", length = length(output)), 
                  C_vector = vector("list", length = length(output)),
                  G = vector("list", length = length(output)),
                  G_plus = vector("list", length = length(output)),
                  e = vector("list", length = length(output)),
                  tau = vector("list", length = length(output)),
                  Z = vector("list", length = length(output)), 
                  proc_matched_Z = vector("list", length = length(output)))
    
    acc_rates = list(alpha = rep(NA, length(output)),
                     e = rep(NA, length(output)),
                     Z = vector("list", length = length(output)))
    
    posterior_means_modes = list(alpha = rep(NA, length(output)),
                                 C = matrix(NA, nrow = length(output), ncol = M),
                                 G_hat = vector("list", length = length(output)),
                                 G_plus_hat = rep(NA, length(output)),
                                 e = rep(NA, length(output)),
                                 tau = vector("list", length = length(output)),
                                 Z = vector("list", length = length(output)))
    
    for (sim in 1:length(output)) {
      
      # ------------------------------------------------------------------------------------------------------------------
      # NETWORK-LEVEL PARAMETERS
      # ------------------------------------------------------------------------------------------------------------------
      
      # ``````````````````````````````````````````````````````````````````````````````````````````````````````````````````
      # AFTER THE MCMC RUN, A SPARSE FINITE MIXTURE IS IDENTIFIED BY POST-PROCESSING THE MCMC DRAWS:
      # ``````````````````````````````````````````````````````````````````````````````````````````````````````````````````
      T_iters_net = dim(output[[sim]]$res_monoLaPCM$store_C_vector)[2] # number of iterations to be post-processed
      
      # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
      # 1. USE THE WADE AND GHARAMANI METHOD TO OBTAIN THE UNIQUE NETWORK-LEVEL CLUSTERING (USING ALL POSTERIOR SAMPLES)
      # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
      store_C_vector = output[[sim]]$res_monoLaPCM$store_C_vector
      
      psm_net = mcclust::comp.psm(t(store_C_vector))
      wade_ghar_net = mcclust.ext::minVI(psm_net, method = "all", cls.draw = t(store_C_vector), include.greedy = F)
      
      posterior_means_modes$C[sim, ] = wade_ghar_net$cl[1, ]
      
      # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
      # 2. DETERMINE THE NUMBER OF NON-EMPTY COMPONENTS
      # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
      G_plus_hat = posterior_means_modes$G_plus_hat[sim] = length(unique(posterior_means_modes$C[sim, ]))
      
      # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
      # 3. A) ISOLATE THE NON-EMPTY COMPONENTS
      # 3. B) USING ONE LS AS A REFERENCE, PROCRUSTES TRANSFORM ALL OTHERS TO ALIGN SO THAT THEY ARE COMPARABLE FOR KMEANS
      # 3. C) ARRANGE THE LATENT SPACES INTO A MATRIX WITH G_PLUS_HAT X T ROWS AND d*N COLUMNS
      # 3. D) CLUSTER THE MATRIX INTO G_PLUS_HAT CLUSTERS USING K-MEANS CLUSTEIRNG
      # NOTE: RESULTS IN A CLASSIFICATION INDICATING TO WHICH CLUSTERS THE COMPONENT-SPECIFIC PARAMS OF EACH DRAW BELONG
      # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
      non_empty_net_clusts = vector("list", length = T_iters_net)
      for (t in 1:T_iters_net) { non_empty_net_clusts[[t]] = sort(unique(store_C_vector[, t])) } # end t for loop
      
      store_tau_active = store_Zs_active = store_mus_active = store_Sigmas_active = store_S_vector_active = 
        store_Kg_active = store_Kg_plus_active = store_pi_active = store_w_active = acc_rate_Z_active = 
        acc_rate_w_active = vector("list", length = T_iters_net)
      
      for (t in 1:T_iters_net) { 
        
        store_tau_active[[t]] = output[[sim]]$res_monoLaPCM$store_tau[non_empty_net_clusts[[t]], t]
        store_Zs_active[[t]] = output[[sim]]$res_monoLaPCM$store_Zs[, , non_empty_net_clusts[[t]], t]
        acc_rate_Z_active[[t]] = output[[sim]]$res_monoLaPCM$store_ARs$acc_rate_Z[non_empty_net_clusts[[t]], t]
        
      } # end t for loop
      
      # 
      
      if (G_plus_hat == 1) {
        
        proc_matched_Z_temp = vector("list", length = length(store_Zs_active))
        ref_space = store_Zs_active[[1]]
        for (t in 1:T_iters_net) {
          proc_matched_Z_temp[[t]] = vegan::procrustes(ref_space, store_Zs_active[[t]])$Yrot
        } # end t for loop
        
      } else {
        
        proc_matched_Z_temp = vector("list", length = length(store_Zs_active))
        ref_space = store_Zs_active[[1]][, , 1]
        for (t in 1:T_iters_net) {
          proc_matched_Z_temp[[t]] = array(NA, dim = dim(store_Zs_active[[t]]))
          for (g in 1:dim(store_Zs_active[[t]])[3]) {
            proc_matched_Z_temp[[t]][, , g] = vegan::procrustes(ref_space, store_Zs_active[[t]][, , g])$Yrot
          } # end g for loop
        } # end t for loop
        
      } # end if else statement
      
      # 
      
      if (G_plus_hat > 1) {
        
        ls_vector = matrix(NA, nrow = 0, ncol = 2*N)
        for (t in 1:T_iters_net) {
          # row_means = rbind(row_means, t(apply(store_Zs_active[[t]], 3, rowMeans))) # original code used row means
          ls_vector = rbind(ls_vector, t(apply(proc_matched_Z_temp[[t]], 3, function(mat) { apply(mat, 2, as.vector) })))
        } # end t for loop 
        kmeans_res_net = kmeans(scale(ls_vector), centers = G_plus_hat, nstart = 100, iter.max = 100)
        png(paste0(plots_path, "/Simulation_", sim, "_LS_kmeans_PPR.png"))
        plot(ls_vector, pch = 19, las = 1, xlab = "Dim. 1", ylab = "Dim. 2", col = kmeans_res_net$cluster)
        dev.off()
        
      } # end if else statement
      
      # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
      # 4. A) FOR EACH ITERATION t = 1,...,T0, CONSTRUCT A CLASSIFICATION SEQUENCE RHO_t OF SIZE G_PLUS_HAT 
      # 4. B) CHECK WHETHER RHO_t IS A PERMUTATION OF (1,...,G_PLUS_HAT); IF NOT, REMOVE THE DRAW t FROM CONSIDERATION
      # 4. C) KEEP ONLY THOSE t THAT ARE PERMUTATIONS FOR ALL PARAMETERS
      # NOTE: THE PROPORTION OF CLASSIFICATION SEQUENCES OF T0 NOT BEING A PERMUTATION IS DENOTED T0_RHO.
      # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
      if (G_plus_hat == 1) {
        
        class_seqs_net = matrix(1, nrow = T_iters_net, ncol = 1)
        
        posterior_means_modes$C[sim, ] = rep(1, M)
        
      } else {
        
        kmeans_res_cluster_temp = kmeans_res_net$cluster
        class_seqs_net = vector("list", length = T_iters_net)
        for (t in 1:T_iters_net) {
          class_seqs_net[[t]] = kmeans_res_cluster_temp[1:dim(store_Zs_active[[t]])[3]]
          kmeans_res_cluster_temp = kmeans_res_cluster_temp[-(1:dim(store_Zs_active[[t]])[3])]
        } # end t for loop
        
      } # end if else statement
      
      # 
      
      perm_indices_net = c()
      for (t in 1:T_iters_net) {
        if (length(class_seqs_net[[t]]) == G_plus_hat) {
          if (all(sort(class_seqs_net[[t]]) == 1:G_plus_hat) == T) {
            perm_indices_net = c(perm_indices_net, t)
          } # end if statement
        } # end if else statement
      } # end t for loop
      
      # check if there are no permutations (then we cannot continue with the algorithm)
      if (length(perm_indices_net) == 0) {
        cat("There are no class. seqs that are permutations of", 1:G_plus_hat, 
            ". Unable to proceed with the algorithm for Simulation", sim, ".\n")
        # move onto the next
        next
      } # end if statement
      params$perm_indices_net[[sim]] = perm_indices_net
      
      class_seqs_keep_net = class_seqs_net[perm_indices_net]
      
      T0_rho_net = T_iters_net - length(perm_indices_net)
      params$T0_rho_net[sim] = T0_rho_net
      
      # THERE ARE NO FURTHER POST-PROCESSING STEPS THAT APPLY TO THE INTERCEPT ALPHA
      
      store_alpha = output[[sim]]$res_monoLaPCM$store_alpha
      store_alpha = store_alpha[perm_indices_net]
      chains$alpha[[sim]] = store_alpha
      posterior_means_modes$alpha[sim] = mean(chains$alpha[[sim]])
      png(paste0(plots_path, "/Simulation_", sim, "_Alpha_Trace_Plot.png"))
      plot(chains$alpha[[sim]], type = "l", las = 1, xlab = "Iteration", ylab = expression(alpha))
      abline(h = posterior_means_modes$alpha[[sim]], col = "magenta", lty = 2, lwd = 2)
      legend("bottomright", legend = "Mean", lty = 2, lwd = 2, col = "magenta", bty = "n")
      dev.off()
      
      # 
      
      store_C_vector = store_C_vector[, perm_indices_net]
      chains$C_vector[[sim]] = store_C_vector
      
      store_G_plus = output[[sim]]$res_monoLaPCM$store_G_plus
      chains$G_plus[[sim]] = store_G_plus
      
      store_G = output[[sim]]$res_monoLaPCM$store_G
      store_G = store_G[perm_indices_net]
      chains$G[[sim]] = store_G
      posterior_means_modes$G_hat[[sim]] = DescTools::Mode(chains$G[[sim]])
      png(paste0(plots_path, "/Simulation_", sim, "_Barplot_G.png"))
      barplot(table(chains$G[[sim]]))
      dev.off()
      
      store_e = output[[sim]]$res_monoLaPCM$store_e
      store_e = store_e[perm_indices_net]
      chains$e[[sim]] = store_e
      posterior_means_modes$e[sim] = mean(chains$e[[sim]])
      png(paste0(plots_path, "/Simulation_", sim, "_e_Trace_Plot.png"))
      plot(chains$e[[sim]], type = "l", las = 1, xlab = "Iteration", ylab = expression("e"[0]))
      abline(h = posterior_means_modes$e[[sim]], col = "magenta", lty = 2, lwd = 2)
      legend("bottomright", legend = "Mean", lty = 2, lwd = 2, col = "magenta", bty = "n")
      dev.off()
      
      store_tau_active = store_tau_active[perm_indices_net]
      store_Zs_active = store_Zs_active[perm_indices_net]
      
      acc_rate_alpha = output[[sim]]$res_monoLaPCM$store_ARs$acc_rate_alpha
      acc_rate_alpha = acc_rate_alpha[perm_indices_net]
      acc_rates$alpha[sim] = mean(acc_rate_alpha)
      
      acc_rate_Z_active = acc_rate_Z_active[perm_indices_net]
      
      # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
      # 5. FOR THE REMAINING T0(1 - T_RHO) DRAWS, A UNIQUE LABELLING IS ACHIEVED BY RE-SORTING THE DRAWS ACCORDING TO THE
      #    CLASSIFICATION SEQUENCES RHO_t (the ones we are keeping)
      # NOTE: THE RE-SORTED, IDENTIFIED DRAWS CAN BE USED FOR FURTHER COMPONENT-SPECIFIC PARAMETER INFERENCE
      # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
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
      
      acc_rate_Z_relabelled = matrix(NA, nrow = G_plus_hat, ncol = length(perm_indices_net))
      for (t in 1:length(perm_indices_net)) {
        for (g in 1:G_plus_hat) {
          acc_rate_Z_relabelled[g, t] = acc_rate_Z_active[[t]][which(class_seqs_keep_net[[t]] == g)] 
        } # end g for loop
      } # end t for loop
      acc_rates$Z[[sim]] = rowMeans(acc_rate_Z_relabelled)
      
      # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
      # APPLY A RELABELLING FUNCTION TO THE POSTERIOR OPTIMAL CLUSTERING SOLUTION SO THAT THERE IS A CORRESPONDENCE BW
      # IT AND THE REMAINING PARAMETERS FOR LATER ANALYSIS
      # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
      WG_order_C = t(apply(store_C_vector_relabelled, 2, function(col) { as.vector(e1071::matchClasses(table(posterior_means_modes$C[sim, ], col), verbose = F)) }))
      
      stacked_WG_C = apply(WG_order_C, 1, paste, collapse = ",")
      stacked_WG_C_freq = table(stacked_WG_C)
      g_match = as.integer(unlist(stringr::str_split(names(stacked_WG_C_freq)[which.max(stacked_WG_C_freq)], ",")))
      
      posterior_means_modes$C[sim, ] = g_match[posterior_means_modes$C[sim, ]]
      
      # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
      # OBTAIN THE POSTERIOR MEAN LATENT POSITION MODEL PARAMETERS (ALPHA AND LATENT SPACES)
      # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
      proc_matched_Z = array(NA, dim = dim(store_Zs_relabelled))
      
      ref_index = 1
      
      for (g in 1:G_plus_hat) {
        
        proc_matched_Z[, , g, 1] = store_Zs_relabelled[, , g, ref_index]
        
        for (t in 2:dim(proc_matched_Z)[4]) {
          proc_Z_fit = vegan::procrustes(proc_matched_Z[, , g, 1], store_Zs_relabelled[, , g, t])
          # proc_matched_Z[, , g, t] = proc_Z_fit$Yrot # pre November 4th 2024
          proc_matched_Z[, , g, t] = sweep(((proc_Z_fit$scale * sweep(store_Zs_relabelled[, , g, t], 2, proc_Z_fit$xmean, "-")) %*% proc_Z_fit$rotation), 2, proc_Z_fit$translation, "+")
        } # end t for loop
      } # end g for loop
      
      chains$proc_matched_Z[[sim]] = proc_matched_Z
      
      posterior_means_modes$Z[[sim]] = vector("list", length = G_plus_hat)
      for (g in 1:G_plus_hat) {
        posterior_means_modes$Z[[sim]][[g]] = apply(proc_matched_Z[, , g, ], 1:2, mean)
      } # end g for loop
      
      # ------------------------------------------------------------------------------------------------------------------
      # PLOTS
      # ------------------------------------------------------------------------------------------------------------------
      cols_update = rep(cols, 3)
      cols_g = vector("list", length = posterior_means_modes$G_plus_hat[sim])
      
      for (g in 1:posterior_means_modes$G_plus_hat[sim]) {
        
        start = ifelse(g == 1, 1, sum(posterior_means_modes$Kg_plus_hat[[sim]][1:(g - 1)]) + 1)
        end = sum(posterior_means_modes$Kg_plus_hat[[sim]][1:g])
        
        cols_g[[g]] = cols_update[start:end]
        
      } # end g for loop
      
      for (g in 1:posterior_means_modes$G_plus_hat[sim]) {
        
        png(paste0(plots_path, "/Simulation_", sim, "_Latent_Space_", g, ".png"))
        
        plot(posterior_means_modes$Z[[sim]][[g]], xlab = expression("z"[1]), ylab = expression("z"[2]), pch = 19, las = 1,
             col = adjustcolor(cols_g[[g]], alpha = 0.8),
             xlim = c(-3, 3), ylim = c(-3, 3))
        grid()
        
        dev.off()
        
      } # end g for loop
      
    } # end sim for loop
    
    return(list(params = params,
                chains = chains,
                acc_rates = acc_rates,
                posterior_means_modes = posterior_means_modes))
    
  } # end post_process_monoLaPCM function