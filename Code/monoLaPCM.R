########################################################################################################################
# FUNCTION TO PERFORM INFERENCE IMPLEMENTING monoLaPCM
########################################################################################################################
monoLaPCM = function(multi, # multiplex in array format with dimensions NxNxM
                  net_mode = NULL, # mode of network: one of "undirected" or "directed"
                  net_type = NULL, # type of network: currently one of "binary", "count" or "pos_real"
                  
                  G_max = NULL, # maximum number of network-level clusters to consider
                  
                  G_0 = NULL, # initial number of network-level clusters
                  a_G = NULL, b_G = NULL, c_G = NULL, # hyperparams for the BNB prior on G-1 (# of net-level comps)
                  
                  l_G = NULL, r_G = NULL, # hyperparameters for the F prior on e (net-level Dirichlet hyperparameter)
                  s_e = NULL, # standard dev. of the proposal distribution of e (net-level Dirichlet hyperparameter)
                  
                  delta_Z = NULL, # scaling factor for the proposal distribution of z (latent positions)
                  
                  m_alpha = NULL, s_alpha = NULL, # hyperparameters for the Normal prior on alpha (LPM intercept)
                  delta_alpha = NULL, # scaling factor for the proposal distribution of alpha (LPM intercept)
                  alpha_prop_sd = NULL, # standard deviation for the proposal distribution of alpha (LPM intercept)
                  
                  beta = NULL, B = NULL, # hyperparameters for the Normal prior on the latent positions
                  
                  thin = NULL, # amount of thinning to incorporate
                  samples = NULL, # number of samples we want at the end of the algorithm (after burn-in and thinning,
                  # but before post-processing, which can reduce the number of samples)
                  
                  cols = NULL, # colours to use for plotting purposes
                  show_plots = NULL, # show plots throughout the algorithm
                  plot_freq = NULL, # how often to show plots throughout the algorithm
                  
                  show_info = NULL, # show the iter number, alpha value and cluster sizes throughout the algorithm
                  log_file = NULL # txt file to fill in with the iteration number, alpha value and cluster sizes 
                  # throughout the algorithm and the acceptance rates at the end of the algorithm
) {
  
  
  # ====================================================================================================================
  # NECESSARY CHECKS
  # ====================================================================================================================
  p = 2 # number of dimensions in the latent space (only considering 2, for now)
  
  # --------------------------------------------------------------------------------------------------------------------
  # ENSURE MULTIPLEX IS IN THE REQUIRED FORMAT
  # --------------------------------------------------------------------------------------------------------------------
  M = dim(multi)[3] # number of networks
  N = dim(multi)[1] # number of nodes (same in each network)
  if (!(is.array(multi) & all(dim(multi) == c(N, N, M)))) {
    stop("mixLPCM requires the multiplex to be in the form of an array with dimensions NxNxM, where N is the number of nodes and M is the number of networks. Please check your formatting and try again.\n")
  } # end if statement
  
  # --------------------------------------------------------------------------------------------------------------------
  # ENSURE NET_MODE IS CORRECTLY IDENTIFIED AS UNDIRECTED/DIRECTED (OR ASSIGN IF NOT SPECIFIED)
  # --------------------------------------------------------------------------------------------------------------------
  if (!exists("net_mode") || is.null(net_mode)) {  # check if net_mode is missing
    
    if (all(apply(multi, 3, isSymmetric))) {
      net_mode = "undirected"
    } else {
      net_mode = "directed"
    } # end if else statement
    
    message(paste0("net_mode not provided. Setting to '", net_mode, "' based on your multiplex input."))
    
  } # end if statement
  
  if (!net_mode %in% c("undirected", "directed")) { # check if net_mode is valid
    stop("Error: 'net_mode' must be one of 'undirected' or 'directed'.")
  } # end if statement
  
  # --------------------------------------------------------------------------------------------------------------------
  # ENSURE NET_TYPE IS CORRECTLY IDENTIFIED AS BINARY/COUNT/POS_REAL (OR ASSIGN IF NOT SPECIFIED)
  # --------------------------------------------------------------------------------------------------------------------
  if (!exists("net_type") || is.null(net_type)) {  # check if net_type is missing
    
    if (all(multi %in% c(0, 1))) {
      net_type = "binary"
    } else if (all(multi >= 0 & multi == floor(multi))) {
      net_type = "count"
    } else if ((all(multi > 0) && any(multi != floor(multi)))) {
      net_type = "pos_real"
    } # end if else statement
    
    message(paste0("net_type not provided. Setting to '", net_type, "' based on your multiplex input."))
    
  } # end if statement
  
  if (!net_type %in% c("binary", "count", "pos_real")) { # check if net_type is valid
    stop("Error: 'net_type' must be one of 'binary', 'count' or 'pos_real'.")
  } # end if statement
  
  # --------------------------------------------------------------------------------------------------------------------
  # CEHCK IF G_MAX IS PROVIDED AND VALID
  # --------------------------------------------------------------------------------------------------------------------
  if (!exists("G_max") || is.null(G_max)) { # check if G_max is missing
    message("G_max not provided. Setting to 5.")
    G_max = 5
  } # end if statement
  
  if (!is.null(G_max) && G_max > M) { # check if G_max is greater than the number of observations to cluster
    stop("Error: 'G_max' must be less than or equal to the number of observations to be clustered, 
         i.e. the number of networks, 'M', in your multiplex.")
  } # end if statement
  
  # --------------------------------------------------------------------------------------------------------------------
  # CHECK IF G_0 IS PROVIDED AND VALID
  # --------------------------------------------------------------------------------------------------------------------
  if (!exists("G_0") || is.null(G_0)) { # check if K_0 is missing
    message("G_0 not provided. Setting to 2.")
    G_0 = 2
  } # end if statement
  
  if (G_0 > M) { # check if G_0 is greater than the number of observations to cluster
    stop("Error: 'G_0' must be less than or equal to the number of observations to be clustered, 
         i.e. the number of networks, 'M', in your multiplex.")
  } # end if statement
  
  # --------------------------------------------------------------------------------------------------------------------
  # CHECK IF A_G, B_G, C_G ARE PROVIDED IN THE MFM MODEL_FRAMEWORK
  # --------------------------------------------------------------------------------------------------------------------
  if (!exists("a_G") || is.null(a_G)) {
    message("a_G not provided. Setting to 8.")
    a_G = 8
  } # end if statement
  if (!exists("b_G") || is.null(b_G)) {
    message("b_G not provided. Setting to 18.")
    b_G = 18
  } # end if statement
  if (!exists("c_G") || is.null(c_G)) {
    message("c_G not provided. Setting to 10.")
    c_G = 10
  } # end if statement
  
  # --------------------------------------------------------------------------------------------------------------------
  # CHECK FOR HYPERPARAMETERS RELATING TO e 
  # --------------------------------------------------------------------------------------------------------------------
  if (!exists("l_G") || is.null(l_G)) {
    message("l_G not provided. Setting to 6.")
    l_G = 6
  } # end if statement
  if (!exists("r_G") || is.null(r_G)) {
    message("r_G not provided. Setting to 3.")
    r_G = 3
  } # end if statement
  
  if (!exists("s_e") || is.null(s_e)) {
    message("s_e not provided. Setting to 4.")
    s_e = 4
  } # end if statement
  
  # --------------------------------------------------------------------------------------------------------------------
  # CHECK FOR THE SCALING FACTOR FOR THE LATENT POSITIONS
  # --------------------------------------------------------------------------------------------------------------------
  if (!exists("delta_Z") || is.null(delta_Z)) { # check if delta_Z is missing
    message("delta_Z not provided. Setting to 0.005.")
    delta_Z = 0.005
  } # end if statement
  
  # --------------------------------------------------------------------------------------------------------------------
  # CHECK FOR THE HYPERPARAMETERS RELATING TO ALPHA (LPM INTERCEPT)
  # --------------------------------------------------------------------------------------------------------------------
  if (!exists("m_alpha") || is.null(m_alpha)) { # check if m_alpha is missing
    message("m_alpha not provided. Setting to 0.")
    m_alpha = 0
  } # end if statement
  
  if (!exists("s_alpha") || is.null(s_alpha)) { # check if s_alpha is missing
    message("s_alpha not provided. Setting to 1.")
    s_alpha = 1
  } # end if statement
  
  if (!exists("delta_alpha") || is.null(delta_alpha)) { # check if delta_alpha is missing
    message("delta_alpha not provided. Setting to 3.")
    delta_alpha = 3
  } # end if statement
  
  if (!exists("alpha_prop_sd") || is.null(alpha_prop_sd)) { # check if alpha_prop_sd is missing
    message("alpha_prop_sd not provided. Setting to 0.01")
    alpha_prop_sd = 0.01
  } # end if statement
  
  # --------------------------------------------------------------------------------------------------------------------
  # CHECK FOR HYPERPARAMETERS RELATING TO THE CLUSTER MEANS, MU
  # --------------------------------------------------------------------------------------------------------------------
  if (!exists("beta") || is.null(beta)) {
    message("beta not provided. Setting to 0-vector.")
    beta = rep(0, p)
  } # end if statement
  
  if (!exists("B") || is.null(B)) {
    message("B not provided. Setting to the identity matrix.")
    B = diag(p)
  } # end if statement
  
  # ====================================================================================================================
  # STORAGE
  # ====================================================================================================================
  iters_keep = samples * thin
  burn_in = 0.3 * iters_keep
  iters = iters_keep + burn_in
  
    # ------------------------------------------------------------------------------------------------------------------
    # NETWORK-LEVEL MIXTURE
    # ------------------------------------------------------------------------------------------------------------------
    store_tau = matrix(NA, nrow = G_max, ncol = iters + 1) # mixing proportions
    
    store_alpha = rep(NA, iters + 1) # LPM intercept
    acc_rate_alpha = rep(NA, iters + 1)
    
    store_Zs = array(NA, c(N, 2, G_max, iters + 1)) # latent positions
    acc_rate_Z = matrix(NA, nrow = G_max, ncol = iters + 1)
    
    store_C_vector = matrix(NA, nrow = M, ncol = iters + 1) # cluster allocations
    store_probs_C = array(NA, dim = c(M, G_max, iters + 1))
    
    store_G = rep(NA, iters + 1) # number of components
    store_G_plus = rep(NA, iters + 1) # number of clusters (active components)
    store_G_indices = matrix(NA, nrow = G_max, ncol = iters + 1) # track which components are active
    
    store_e = rep(NA, iters + 1) # Dirichlet hyperparameter 
    acc_rate_e = rep(NA, iters + 1)
    
  # ====================================================================================================================
  # INITIALISATIONS
  # ====================================================================================================================
  
    # ------------------------------------------------------------------------------------------------------------------
    # NETWORK-LEVEL MIXTURE
    # ------------------------------------------------------------------------------------------------------------------
    store_G[1] = initial_G = G_0 # initial number of components
    store_G_plus[1] = initial_Gplus = G_0 # initial number of clusters (all active at t = 0)
    store_G_indices[, 1] = initial_G_indices = c(1:G_0, rep(0, G_max - G_0))
    
    store_tau[, 1] = initial_tau = c(rep(1/G_0, G_0), rep(0, G_max - G_0)) # mixing proportions evenly spread out
    
    store_e[1] = initial_e = 0.00001
    
    # CALCULATE THE DISTANCES BETWEEN ALL NETWORKS IN THE MULTIPLEX (TO OBTAIN INITIAL CLUSTER ALLOCATIONS)
    dist_mat = matrix(NA, M, M)
    for (m1 in 1:M) {
      for (m2 in 1:M) {
        if (net_type == "count") {
          
          gra1 = igraph::graph_from_adjacency_matrix(multi[, , m1], mode = net_mode, weighted = T)
          gra2 = igraph::graph_from_adjacency_matrix(multi[, , m2], mode = net_mode, weighted = T)
          
          el1 = cbind(igraph::as_edgelist(gra1), igraph::E(gra1)$weight)
          el2 = cbind(igraph::as_edgelist(gra2), igraph::E(gra2)$weight)
          dist_mat[m1, m2] = nature_dist(el1, el2, gra1, gra2, 0.45, 0.45, 0.1, net_type)
          
        } else if (net_type == "binary") {
          
          gra1 = igraph::graph_from_adjacency_matrix(multi[, , m1], mode = net_mode)
          gra2 = igraph::graph_from_adjacency_matrix(multi[, , m2], mode = net_mode)
          
          el1 = igraph::as_edgelist(gra1)
          el2 = igraph::as_edgelist(gra2)
          dist_mat[m1, m2] = nature_dist(el1, el2, gra1, gra2, 0.45, 0.45, 0.1, net_type)
          
        } else if (net_type == "pos_real") {
          
          gra1 = igraph::graph_from_adjacency_matrix(multi[, , m1], mode = net_mode, weighted = T)
          gra2 = igraph::graph_from_adjacency_matrix(multi[, , m2], mode = net_mode, weighted = T)
          
          el1 = igraph::as_edgelist(gra1)
          el2 = igraph::as_edgelist(gra2)
          dist_mat[m1, m2] = nature_dist(el1, el2, gra1, gra2, 0.45, 0.45, 0.1, net_type)
          
        } else {
          cat("\nHelp! I can only handle binary, count or (positive) real networks right now...\n")
          stop()
        } # end net_type check
        rm(gra1, gra2, el1, el2)
      } # end m2 loop
      rm(m2)
    } # end m1 loop
    rm(m1)
    
    mds_init_C = cmdscale(dist_mat, k = 2, eig = T) # 2D representation of dist_mat using MDS
    kmeans_init_C = kmeans(mds_init_C$points, centers = G_0, nstart = 1000) 
    store_C_vector[, 1] = initial_C_vector = kmeans_init_C$cluster
    rm(mds_init_C, kmeans_init_C)
    
    M_comps = rep(0, store_G[1]) # size of network-level clusters
    for (g in 1:store_G[1]) { M_comps[g] = sum(store_C_vector[, 1] == g) }; rm(g)
    store_M_comps = c(M_comps, rep(0, G_max - length(M_comps)))
    
    # OBTAIN THE INITIAL LATENT POSITIONS
    for (g in 1:G_0) {
      g_inds = which(store_C_vector[, 1] == g) # networks in component g
      if (length(g_inds) == 1) { # if there is only one network in cluster g
        gather_geodist_Z = matrix(NA, nrow = N, ncol = N) 
        geo_dist = sna::geodist(multi[, , g_inds[1]], count.paths = F, inf.replace = 5)$gdist # geodesic distances between 
        # all nodes in a network
        gather_geodist_Z = geo_dist
        store_Zs[, , g, 1] = cmdscale(gather_geodist_Z, k = 2, eig = T)$points # 2D representation of the geodesic dists
      } else { # more than one network in cluster g
        gather_geodist_Z = array(NA, c(N, N, length(g_inds))) 
        for (m in 1:length(g_inds)) {
          geo_dist = sna::geodist(multi[, , g_inds[m]], count.paths = F, inf.replace = 5)$gdist # geodesic distances between 
          # all nodes in a network
          gather_geodist_Z[, , m] = geo_dist
        } # end m loop
        avg_geodist = apply(gather_geodist_Z, 1:2, mean) # average geodesic distance of all networks in cluster g
        store_Zs[, , g, 1] = cmdscale(avg_geodist, k = 2, eig = T)$points # 2D representation of the geodesic dists
      } # end if else statement
      rm(gather_geodist_Z, g_inds, geo_dist, avg_geodist)
    } # end g loop
    rm(g)
    initial_latent_spaces_original = store_Zs[, , , 1]
    initial_latent_spaces = initial_latent_spaces_original
    
    for (g in 1:G_0) { # adding some noise so that each chain has a different starting point
      initial_latent_spaces[, , g] = initial_latent_spaces_original[, , g] + 
        rnorm(N * 2, 0, sd(as.vector(initial_latent_spaces_original[, , g])) / 100) # empirical variance / 100
      store_Zs[, , g, 1] = initial_latent_spaces[, , g]
    } # end g for loop
    rm(g)
    
    # CALCULATE THE DISTANCE MATRICES FROM THE INITIAL LATENT SPACES
    D = array(NA, c(N, N, G_max))
    for (g in 1:G_max) {
      D[, , g] = Rfast::Dist(store_Zs[, , g, 1], square = T)
    } # end g for loop
    rm(g)
    
    # OBTAIN THE INITIAL LPM INTERCEPT
    coeffs_init = rep(NA, G_0)
    for (g in 1:G_0) {
      g_inds = which(store_C_vector[, 1] == g)
      coeffs_g = rep(NA, length(g_inds))
      for (m in 1:length(g_inds)) {
        if (net_type == "count") {
          coeffs_g[m] = log(mean(multi[, , g_inds[m]])) + mean(D[, , g])
        } else if (net_type == "binary") {
          coeffs_g[m] = LaplacesDemon::logit(mean(multi[, , g_inds[m]])) + mean(D[, , g])
        } else if (net_type == "pos_real") {
          coeffs_g[m] = -(1 / (mean(multi[, , g_inds[m]]))) + mean(D[, , g])
        } else {
          cat("Uh oh. There was a problem initialising alpha.")
        } # end net_type check
      } # end m loop
      rm(m)
      coeffs_init[g] = mean(coeffs_g)
    } # end g loop
    rm(g, coeffs_g, g_inds)
    
    store_alpha[1] = initial_alpha = mean(coeffs_init)
    rm(coeffs_init)
  
    # --------------------------------------------------------------------------------------------------------------------
    # INITIAL PLOTS
    # --------------------------------------------------------------------------------------------------------------------
    par(mar = c(2, 4, 2, 2) + 0.1)
    
    if (show_plots == T & G_0 <= 4) {
      
      par(mfrow = c(2, 2)) # set up a 2x2 plotting layout
      
      for (g in 1:G_0) {
        
        plot(store_Zs[, , g, 1], col = adjustcolor("grey", alpha = 0.8), pch = 19, las = 1,
             xlab = expression("z"[1]), ylab = expression("z"[2]), 
             # xlim = c(-3, 3), ylim = c(-3, 3), 
             main = paste("Initial Latent Space", g))
        grid()
        
      } # end g for loop
      rm(g)
      
      # fill remaining slots if G_0 < 4 to ensure new plots appear on a fresh screen next time
      for (g in (G_0 + 1):4) {
        plot.new()  # creates an empty plot
      } # end g for loop
      rm(g)
      
    } # end if statement
    
  # --------------------------------------------------------------------------------------------------------------------
  # WRITE TO LOG FILE
  # --------------------------------------------------------------------------------------------------------------------
  write("Initialisation: Done\n", file = log_file, append = T)
  
  ######################################################################################################################
  # MCMC ALGORITHM
  ######################################################################################################################
  
  prog_bar = progress::progress_bar$new(format = "[:bar] :percent (ETR: :eta)", total = iters, clear = F, width = 60)
  
  start_time = Sys.time()
  
  for (t in 1:iters) {
    
    # ==================================================================================================================
    # STEP 1: UPDATE THE NETWORK PARTITION/ALLOCATIONS
    # ==================================================================================================================
    
    # __________________________________________________________________________________________________________________
    # (A) SAMPLE THE ALLOCATION VARIABLE FOR EVERY NETWORK
    # __________________________________________________________________________________________________________________
    store_Cmat = matrix(0, nrow = M, ncol = G_max)
    
    for (m in 1:M) {
      
      # ----------------------------------------------------------------------------------------------------------------
      # (i) Calculate the log normalisation constant
      # ----------------------------------------------------------------------------------------------------------------
      log_norm_const_C_vec = rep(NA, store_G[t])
      for (g in 1:store_G[t]) {
        
        log_norm_const_C_vec[g] = log(store_tau[g, t]) + log_LPM_fast(multi[, , m], store_alpha[t], D[, , g], 
                                                                      net_type, net_mode)
        
        if (is.infinite(log_norm_const_C_vec[g])) {
          log_norm_const_C_vec[g] = log(store_tau[g, t]) + log_LPM(multi[, , m], store_alpha[t], D[, , g], 
                                                                   net_type, net_mode)
        } # end if statement
        
      } # end g loop
      rm(g)
      
      log_norm_const_C = matrixStats::logSumExp(log_norm_const_C_vec)
      rm(log_norm_const_C_vec)
      
      # ----------------------------------------------------------------------------------------------------------------
      # (ii) Calculate the log probabilities
      # ----------------------------------------------------------------------------------------------------------------
      log_probs_C = rep(NA, store_G[t])
      for (g in 1:store_G[t]) {
        
        log_probs_C[g] = log(store_tau[g, t]) + log_LPM_fast(multi[, , m], store_alpha[t], D[, , g], net_type, net_mode) - log_norm_const_C
        
        if (is.infinite(log_probs_C[g])) {
          log_probs_C[g] = log(store_tau[g, t]) + log_LPM(multi[, , m], store_alpha[t], D[, , g], net_type, net_mode) - log_norm_const_C
        } # end if statement
        
      } # end g loop
      rm(g)
      
      # ----------------------------------------------------------------------------------------------------------------
      # (iii) Calculate the probabilities
      # ----------------------------------------------------------------------------------------------------------------
      probs_C = exp(log_probs_C)
      store_probs_C[m, 1:store_G[t], t + 1] = probs_C
      rm(log_probs_C)
      
      if (!all.equal(sum(probs_C), 1)) {
        cat("Uh oh, your network allocation probabilities sum to", sum(probs_C), "when they should sum to 1!\n", 
            "Check on network", m, ".\n")
        stop()
      } # end probs sum if statement 
      
      # ----------------------------------------------------------------------------------------------------------------
      # (iv) Obtain the allocation variable for the current network
      # ----------------------------------------------------------------------------------------------------------------
      C_m = which(rmultinom(1, 1, probs_C) == 1)
      store_C_vector[m, t + 1] = C_m
      store_Cmat[m, C_m] = 1
      
    } # end m loop
    rm(m, C_m)
    
    # __________________________________________________________________________________________________________________
    # (B) DETERMINE RELEVANT QUANTITIES AND RELABEL
    # __________________________________________________________________________________________________________________
    
    # ------------------------------------------------------------------------------------------------------------------
    # (i) Determine M_g (cluster sizes)
    # ------------------------------------------------------------------------------------------------------------------
    M_comps = rep(0, store_G[t])
    for (g in 1:store_G[t]) { M_comps[g] = sum(store_C_vector[, t + 1] == g) }; rm(g)
    store_M_comps = c(M_comps, rep(0, G_max - length(M_comps)))
    
    # ------------------------------------------------------------------------------------------------------------------
    # (ii) Determine G_+ (number of clusters or active components)
    # ------------------------------------------------------------------------------------------------------------------
    store_G_plus[t + 1] = sum(store_M_comps > 0)
    
    # ------------------------------------------------------------------------------------------------------------------
    # (iii) Relabel such that the first G_+ components are non-empty (anything with a g index)
    # ------------------------------------------------------------------------------------------------------------------
    store_G_indices[, t + 1] = c(which(store_M_comps > 0), which(store_M_comps == 0))
    
    store_M_comps = store_M_comps[store_G_indices[, t + 1]]
    
    store_Zs[, , , t] = store_Zs[, , store_G_indices[, t + 1], t]
    store_tau[, t] = store_tau[store_G_indices[, t + 1], t]
    store_Cmat = store_Cmat[, store_G_indices[, t + 1]]
    for (m in 1:M) { store_C_vector[m, t + 1] = which(store_Cmat[m, ] == 1) }; rm(m)
    
    acc_rate_Z[, t] = acc_rate_Z[store_G_indices[, t + 1], t]
    
    # ==================================================================================================================
    # STEP 2: UPDATE THE PARAMETERS OF THE NON-EMPTY COMPONENTS
    # ==================================================================================================================
    
    # __________________________________________________________________________________________________________________
    # (A) SAMPLE THE COMPONENT-SPECIFIC PARAMETERS, I.E. THE LATENT SPACES
    # __________________________________________________________________________________________________________________
    
    for (g in 1:store_G_plus[t + 1]) {
      
      # ----------------------------------------------------------------------------------------------------------------
      # (i) Propose a new Z_g (block update) [from the node-level mixture]
      # ----------------------------------------------------------------------------------------------------------------
      Z_g_prop = store_Zs[, , g, t] + MASS::mvrnorm(N, beta, (delta_Z^2) * B)
      
      if (any(is.na(Z_g_prop))) {
        cat("Your proposed latent space for g =", g, "was not filled in and remains with some NAs.\n")
      } # end NA if statement
      
      # ----------------------------------------------------------------------------------------------------------------
      # # (ii) Calculate the acceptance ratio
      # ----------------------------------------------------------------------------------------------------------------
      log_FC_prop = log_FC_Z_g_monoLaPCM(multi, store_C_vector[, t + 1], g, store_alpha[t], Z_g_prop, D[, , g], beta, B, net_type, net_mode)
      log_FC_current = log_FC_Z_g_monoLaPCM(multi, store_C_vector[, t + 1], g, store_alpha[t], store_Zs[, , g, t], D[, , g], beta, B, net_type, net_mode)
      log_dens_prop = log_dens_Z_g_monoLaPCM(Z_g_prop, store_Zs[, , g, t], (delta_Z ^ 2) * B)
      log_dens_current = log_dens_Z_g_monoLaPCM(store_Zs[, , g, t], Z_g_prop, (delta_Z ^ 2) * B)
      log_acc_ratio = (log_FC_prop - log_FC_current) + (log_dens_current - log_dens_prop)
      rm(log_FC_prop, log_FC_current, log_dens_current, log_dens_prop)
      
      # ----------------------------------------------------------------------------------------------------------------
      # (iii) Accept or reject
      # ----------------------------------------------------------------------------------------------------------------
      
      # ````````````````````````````````````````````````````````````````````````````````````````````````````````````````
      # OPTION 2: USING OFFLINE PROCRUSTES TRANSFORMATIONS
      # ````````````````````````````````````````````````````````````````````````````````````````````````````````````````
      lu = log(runif(1))
      if (lu < log_acc_ratio) { # accept
        store_Zs[, , g, t + 1] = Z_g_prop
        acc_rate_Z[g, t + 1] = 1
      } else { # reject
        store_Zs[, , g, t + 1] = store_Zs[, , g, t]
        acc_rate_Z[g, t + 1] = 0
      } # end if else statement
      rm(lu)
      
    } # end g loop
    rm(g, Z_g_prop)
    
    # __________________________________________________________________________________________________________________
    # CALCULATE THE DISTANCE MATRICES FROM THE UPDATED LATENT SPACES
    # __________________________________________________________________________________________________________________
    D = array(NA, c(N, N, G_max))
    for (g in 1:G_max) {
      D[, , g] = Rfast::Dist(store_Zs[, , g, t + 1], square = T)
    } # end g for loop
    rm(g)
    
    # __________________________________________________________________________________________________________________
    # (B) SAMPLE THE HYPERPARAMETERS AND ANY NON-COMPONENT-SPECIFIC PARAMETERS
    # __________________________________________________________________________________________________________________
    
    # ------------------------------------------------------------------------------------------------------------------
    # Step 5: Sample the non-component-specific parameter, i.e. alpha
    # ------------------------------------------------------------------------------------------------------------------
    
    # ..................................................................................................................
    # (a) Hyperparameters for the informed proposal distribution
    # ..................................................................................................................
    # not doing this anymore
    
    # ..................................................................................................................
    # (b) Propose a new alpha
    # ..................................................................................................................
    alpha_prop = rnorm(1, store_alpha[t], delta_alpha * alpha_prop_sd)
    
    # ..................................................................................................................
    # (c) Calculate the acceptance ratio
    # ..................................................................................................................
    log_FC_prop = log_FC_alpha(alpha_prop, M, store_G_plus[t + 1], store_Cmat, multi, store_Zs[, , , t + 1], 
                               m_alpha, s_alpha, net_type, net_mode)
    log_FC_current = log_FC_alpha(store_alpha[t], M, store_G_plus[t + 1], store_Cmat, multi, store_Zs[, , , t + 1], 
                                  m_alpha, s_alpha, net_type, net_mode)
    log_dens_prop = dnorm(alpha_prop, store_alpha[t], delta_alpha * alpha_prop_sd, log = T)
    log_dens_current = dnorm(store_alpha[t], alpha_prop, delta_alpha * alpha_prop_sd, log = T)
    log_acc_ratio = (log_FC_prop - log_FC_current) + (log_dens_current - log_dens_prop)
    rm(log_FC_prop, log_FC_current, log_dens_current, log_dens_prop)
    
    # ..................................................................................................................
    # (d) Accept or reject
    # ..................................................................................................................
    lu = log(runif(1))
    if (lu < log_acc_ratio) {
      store_alpha[t + 1] = alpha_prop
      acc_rate_alpha[t + 1] = 1
    } else {
      store_alpha[t + 1] = store_alpha[t]
      acc_rate_alpha[t + 1] = 0
    } # end if else statement
    rm(lu, alpha_prop)
    
    # ==================================================================================================================
    # STEP 3: UPDATE G and e
    # ==================================================================================================================
    
    # __________________________________________________________________________________________________________________
    # (A) UPDATE G
    # __________________________________________________________________________________________________________________
    
    # ------------------------------------------------------------------------------------------------------------------
    # (i) Calcualte the log normalisation constant
    # ------------------------------------------------------------------------------------------------------------------
    log_norm_const_G_vec = rep(NA, length(store_G_plus[t + 1]:G_max))
    options_G = store_G_plus[t + 1]:G_max
    for (g in 1:length(options_G)) {
      log_norm_const_G_vec[g] = log_post_comps_G(options_G[g], a_G, b_G, c_G, store_G_plus[t + 1], store_e[t], 
                                                 store_M_comps)
    } # end g for loop
    rm(g)
    
    log_norm_const_G = matrixStats::logSumExp(log_norm_const_G_vec)
    rm(log_norm_const_G_vec)
    
    # ------------------------------------------------------------------------------------------------------------------
    # (ii) Calculate the log probabilities
    # ------------------------------------------------------------------------------------------------------------------
    log_probs_G = rep(NA, length(store_G_plus[t + 1]:G_max))
    for (g in 1:length(options_G)) {
      log_probs_G[g] = log_post_comps_G(options_G[g], a_G, b_G, c_G, store_G_plus[t + 1], store_e[t], store_M_comps) - 
        log_norm_const_G
    } # end g for loop
    rm(g)
    
    # ------------------------------------------------------------------------------------------------------------------
    # (iii) Calculate the probabilities
    # ------------------------------------------------------------------------------------------------------------------
    probs_G = exp(log_probs_G)
    rm(log_probs_G)
    
    if (!all.equal(sum(probs_G), 1)) {
      cat("Uh oh, your G probabilities sum to", sum(probs_G), "when they should sum to 1!\n")
      stop()
    } # end probs sum if statement
    
    # ------------------------------------------------------------------------------------------------------------------
    # (iv) Obtain the new G
    # ------------------------------------------------------------------------------------------------------------------
    ind_G = which(rmultinom(1, 1, probs_G) == 1)
    store_G[t + 1] = options_G[ind_G]
    rm(ind_G)
    
    # __________________________________________________________________________________________________________________
    # (B) UPDATE e
    # __________________________________________________________________________________________________________________
    
    # ------------------------------------------------------------------------------------------------------------------
    # (i) Propose a new e
    # ------------------------------------------------------------------------------------------------------------------
    log_e_prop = rnorm(1, log(store_e[t]), s_e)
    e_prop = exp(log_e_prop)
    
    # ------------------------------------------------------------------------------------------------------------------
    # (ii) Calculate the acceptance ratio
    # ------------------------------------------------------------------------------------------------------------------
    log_FC_prop = log_post_e(e_prop, l_G, r_G, store_G_plus[t + 1], store_M_comps, store_G[t + 1])
    log_FC_current = log_post_e(store_e[t], l_G, r_G, store_G_plus[t + 1], store_M_comps, store_G[t + 1])
    log_dens_prop = dnorm(log_e_prop, log(store_e[t]), s_e, log = T)
    log_dens_current = dnorm(log(store_e[t]), log_e_prop, s_e, log = T)
    log_acc_ratio = (log_FC_prop - log_FC_current) + (log_dens_current - log_dens_prop)
    rm(log_FC_prop, log_FC_current, log_dens_current, log_dens_prop)
    
    # ------------------------------------------------------------------------------------------------------------------
    # (iii) Accept or reject
    # ------------------------------------------------------------------------------------------------------------------
    lu = log(runif(1))
    if (lu < log_acc_ratio) {
      store_e[t + 1] = e_prop
      acc_rate_e[t + 1] = 1
    } else {
      store_e[t + 1] = store_e[t]
      acc_rate_e[t + 1] = 0
    } # end if else statement
    rm(lu, log_e_prop, e_prop)
    
    # ==================================================================================================================
    # STEP 4: ADD G - G_+ EPMTY COMPONENTS AND UPDATE TAU
    # ==================================================================================================================
    
    # __________________________________________________________________________________________________________________
    # (A) ADD G - G_+ EMPTY COMPONENTS
    # __________________________________________________________________________________________________________________
    
    if (store_G[t + 1] > store_G_plus[t + 1]) {
      
      # ----------------------------------------------------------------------------------------------------------------
      # (i) M_g
      # ----------------------------------------------------------------------------------------------------------------
      store_M_comps = store_M_comps
      
      # ----------------------------------------------------------------------------------------------------------------
      # (ii) Sample component parameters from the priors
      # ----------------------------------------------------------------------------------------------------------------
      for (g in (store_G_plus[t + 1] + 1):store_G[t + 1]) {
        store_Zs[, , g, t + 1] = MASS::mvrnorm(N, beta, B)
      } # end g for loop
      rm(g)
      
      D[, , g] = Rfast::Dist(store_Zs[, , g, t + 1], square = T)
      
    } # end if else statement adding extras
    
    # __________________________________________________________________________________________________________________
    # (B) SAMPLE TAU
    # __________________________________________________________________________________________________________________
    zetas = (store_e[t + 1] / store_G[t + 1]) + store_M_comps[1:store_G[t + 1]]
    
    if (store_G[t + 1] == 1) {
      store_tau[1:store_G[t + 1], t + 1] = LaplacesDemon::rdirichlet(1, c(zetas, 0.001))[1]
    } else {
      store_tau[1:store_G[t + 1], t + 1] = LaplacesDemon::rdirichlet(1, zetas)
    } # end rdirichlet catch if statement
    rm(zetas)
    
    store_tau[, t + 1] = ifelse(store_tau[, t + 1] == 0, .Machine$double.eps, store_tau[, t + 1]) # prob of exactly 0 is problematic in rdirichlet
    
    # ==================================================================================================================
    # ITERATION COUNTER
    # ==================================================================================================================
    if (t %% plot_freq == 0) {
      
      if (show_info == T) {
        
        cat("\nIteration:", t, "\n")
        cat("--------------------------------------\n")
        cat("alpha:", round(store_alpha[t + 1], 3), "\n")
        cat("G+:", store_G_plus[t + 1], "\n") 
        cat("# networks in each cluster:", store_M_comps[1:store_G_plus[t + 1]], "\n")
        cat("\n\n")
        
        write(paste("\nIteration:", t), file = log_file, append = T)
        write(paste("--------------------------------------"), file = log_file, append = T)
        write(paste("alpha:", round(store_alpha[t + 1], 3)), file = log_file, append = T)
        write(paste("G+:", store_G_plus[t + 1]), file = log_file, append = T)
        write(paste(c("# networks in each cluster:", store_M_comps[1:store_G_plus[t + 1]]), collapse = " "), file = log_file, append = T)
        write(paste("\n"), file = log_file, append = T)
        
      } # end if statement
      
    } # end if statement
    
    # ==================================================================================================================
    # GET RID OF ANYTHING I NO LONGER NEED
    # ==================================================================================================================
    rm() 
    
    prog_bar$tick()
    
  } # end t loop
  
  end_time = Sys.time()
  
  write(paste("monoLaPCM fit: Done\n\n"), file = log_file, append = T)
  
  # ====================================================================================================================
  # REMOVE BURN-IN AND THIN
  # ====================================================================================================================
  {
    # ------------------------------------------------------------------------------------------------------------------
    # NETWORK-LEVEL MIXTURE
    # ------------------------------------------------------------------------------------------------------------------
    post_burn_in_indices = c(1:(iters))[-c(1:(burn_in + 1))]
    post_thinning_indices = post_burn_in_indices[seq(1, length(post_burn_in_indices), thin)]
    store_G_indices = store_G_indices[, post_thinning_indices]
    store_alpha = store_alpha[post_thinning_indices]
    acc_rate_alpha = acc_rate_alpha[post_thinning_indices]
    store_C_vector = store_C_vector[, post_thinning_indices]
    store_e = store_e[post_thinning_indices]
    acc_rate_e = acc_rate_e[post_thinning_indices]
    store_G = store_G[post_thinning_indices]
    store_G_plus = store_G_plus[post_thinning_indices]
    store_tau = store_tau[, post_thinning_indices]
    store_Zs = store_Zs[, , , post_thinning_indices]
    acc_rate_Z = acc_rate_Z[, post_thinning_indices]
  }
  
  store_ARs = list(acc_rate_alpha = acc_rate_alpha, 
                   acc_rate_Z = acc_rate_Z, 
                   acc_rate_e = acc_rate_e)
  
  initialisations = list(initial_G = initial_G,
                         initial_Gplus = initial_Gplus,
                         initial_G_indices = initial_G_indices,
                         initial_tau = initial_tau, 
                         initial_e = initial_e, 
                         initial_C_vector = initial_C_vector, 
                         initial_latent_spaces = initial_latent_spaces, 
                         initial_alpha = initial_alpha)
  
  return(list(multi = multi,
              store_G_indices = store_G_indices,
              store_alpha = store_alpha,
              store_C_vector = store_C_vector,
              store_probs_C = store_probs_C,
              store_G = store_G,
              store_G_plus = store_G_plus,
              store_e = store_e,
              store_tau = store_tau,
              store_Zs = store_Zs,
              store_ARs = store_ARs,
              time_taken = end_time - start_time,
              initialisations = initialisations))
  
} # end monoLaPCM