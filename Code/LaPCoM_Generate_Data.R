########################################################################################################################
# FUNCTION TO GENERATE THE MULTIPLEX DATA FOR LaPCoM
########################################################################################################################
simulate_data = function(net_type, M, N, alpha, G, K, tau, pi, log_file) {
  
  # ====================================================================================================================
  # ====================================================================================================================
  if ((net_type %in% c("count", "binary")) == F) {
    cat("Uh oh. Please choose the type of network you wish to simulate from (count, binary).")
  } # end network type check
  
  # ====================================================================================================================
  # STORAGE SET UP
  # ====================================================================================================================
  multi = array(NA, c(N, N, M))
  C_vector = rep(NA, M)
  C_mat = matrix(0, nrow = M, G)
  
  S_list = vector("list", G)
  S_vec_list = vector("list", G)
  mu_list = vector("list", G)
  Z_array = array(NA, c(N, 2, G))
  
  # ====================================================================================================================
  # GENERATE THE NODE-LEVEL CLUSTER ALLOCATIONS
  # ====================================================================================================================
  for (g in 1:G) {
    pi_g = pi[[g]]
    S_list[[g]] = t(rmultinom(N, 1, pi_g))
    S_vec_list[[g]] = apply(S_list[[g]], 1, function(row) { which(row == 1) })
  } # end g loop
  
  # ====================================================================================================================
  # GENERATE THE MEAN OPTIONS FOR THE NODE-LEVEL CLUSTERS
  # ====================================================================================================================
  for (g in 1:G) {
    mu_list[[g]] = matrix(NA, nrow = K[g], ncol = 2)
    if (K[g] == 1) {
      mu_list[[g]] = c(0, 0)
    } else if (K[g] == 2) {
      mu_list[[g]][1, ] = c(-0.8, 0.8)
      mu_list[[g]][2, ] = c(0.8, -0.8)
    } else if (K[g] == 3) {
      mu_list[[g]][1, ] = c(-0.9, -0.9) # used to all be 0.8s
      mu_list[[g]][2, ] = c(1.4, 0.4)
      mu_list[[g]][3, ] = c(-0.9, 1.4)
    } # end if else check
  } # end g for loop
  
  # ====================================================================================================================
  # GENERATE THE LATENT SPACES
  # ====================================================================================================================
  for (g in 1:G) {
    if (K[g] == 1) {
      Z_array[, , g] =  MASS::mvrnorm(N, mu_list[[g]], Sigma = 0.25 * Rfast::Diag.matrix(2, 1))
    } else {
      for (k in 1:K[g]) {
        Z_array[which(S_vec_list[[g]] == k), , g] = MASS::mvrnorm(length(which(S_vec_list[[g]] == k)), 
                                                                  mu_list[[g]][k, ], 
                                                                  Sigma = 0.25 * Rfast::Diag.matrix(2, 1))
      } # end k loop
    } # end if else statement about K[g]
  } # end g loop
  
    # ------------------------------------------------------------------------------------
    # LOOK AT THE LATENT SPACES
    # ------------------------------------------------------------------------------------
    par(mfrow = c(2, 1))  
    for (g in 1:G) { plot(Z_array[, , g], col = S_vec_list[[g]], pch = 16) }
  
  # ====================================================================================================================
  # FILL IN THE MULTIPLEX
  # ====================================================================================================================
  for (m in 1:M) {
    
    # ------------------------------------------------------------------------------------
    # ASSIGN NETWORK m TO A CLUSTER
    # ------------------------------------------------------------------------------------
    C_vector[m] = which(rmultinom(1, 1, tau) == 1)
    C_mat[m, C_vector[m]] = 1
    
    # ------------------------------------------------------------------------------------
    # DISTANCE MATRIX
    # ------------------------------------------------------------------------------------
    D = as.matrix(Rfast::Dist(Z_array[, , C_vector[m]], square = T))
    
    # ------------------------------------------------------------------------------------
    # GENERATE A NETWORK
    # ------------------------------------------------------------------------------------
    Y_m = matrix(0, N, N)
    eta = alpha - D
    if (net_type == "count") {
      lambda_g = exp(eta)
      Y_m[lower.tri(Y_m)] = rpois((N * (N - 1)) / 2, Rfast::lower_tri(lambda_g))
    } else if (net_type == "binary") {
      pr = (exp(eta)) / (1 + exp(eta))
      Y_m[lower.tri(Y_m)] = rbinom((N * (N - 1)) / 2, 1, Rfast::lower_tri(pr)) 
    } else {
      cat("Uh oh. We can't generate a network of this type.")
    }
    Y_m = Y_m + Rfast::transpose(Y_m)
    
    # ------------------------------------------------------------------------------------
    # STORE THE NETWORK
    # ------------------------------------------------------------------------------------
    multi[, , m] = Y_m
    
  } # end m loop
  
  # ====================================================================================================================
  # VIEW THE NETWORKS IN THEIR LATENT SPACE COORDINATE SYSTEM
  # ====================================================================================================================
  par(mfrow = c(5, M/5))
  # par(mfrow = c(2, 4))
  par(mar = rep(0, 4))
  for (g in 1:G) {
    net_inds = which(C_vector == g)

    for (m in net_inds) {
      plot(network::as.network(multi[, , m], directed = F),
           vertex.col = c("magenta", "purple", "forestgreen")[S_vec_list[[C_vector[m]]]],
           coord = Z_array[, , C_vector[m]],
           edge.col = "grey", vertex.cex = 1.5)
      
    } # end m for loop
  } # end g for loop
  
  # ====================================================================================================================
  # VIEW THE NETWORKS AS AN OVERLAY OF THE LATENT SPACE COORDINATE SYSTEM
  # ====================================================================================================================
  par(mfrow = c(5, M/5))
  par(mar = c(0.5, 0.25, 0.25, 0.5) + 0.1)
  for (g in 1:G) {
    net_inds = which(C_vector == g)
    
    for (m in net_inds) {
      z1 = Z_array[, , g][, 1]
      z2 = Z_array[, , g][, 2]
      ones = which(multi[, , m] != 0)
      
      plot(Z_array[, , g], type = "n", las = 1, xlim = c(-2.5, 2.5), ylim = c(-2.5, 2.5))
      grid()
      segments(z1[row(multi[, , m])[ones]], z2[row(multi[, , m])[ones]], 
               z1[col(multi[, , m])[ones]], z2[col(multi[, , m])[ones]], 
               lwd = 1, col = adjustcolor("grey", 0.4))
      points(Z_array[, , g], 
             col = c("magenta", "purple", "forestgreen")[S_vec_list[[g]]], pch = 19, cex = 1)
      
    } # end m for loop
  } # end g for loop
  
  # ====================================================================================================================
  # ENSURE THAT EVERY NODE HAS AT LEAST ONE CONNECTION IN EACH LATENT SPACE
  # ====================================================================================================================
  disc_flag = F
  for (g in 1:G) {
    net_inds = which(C_vector == g)
    if (any(rowSums(apply(multi[, , net_inds], 3, rowSums)) == 0)) {
      cat("Node(s)", which(rowSums(apply(multi[, , net_inds], 3, rowSums)) == 0), 
          "are disconnected across all networks in cluster", g, ".\n")
      disc_flag = T
    } # end if statement
  } # end g for loop
  
  # ====================================================================================================================
  # WRITE TO LOG FILE
  # ====================================================================================================================
  write(paste("Data Simulation: Done"), file = log_file, append = T)
  
  # ====================================================================================================================
  # STORAGE
  # ====================================================================================================================
  true_params = list(alpha = alpha, tau = tau, 
                     S_list = S_list, S_vec_list = S_vec_list,
                     mu_list = mu_list, 
                     C_mat = C_mat, C_vector = C_vector, Z_array = Z_array, K = K,
                     net_type = net_type)
  
  return(list(multi = multi, true_params = true_params, disc_flag = disc_flag))
} # end simulate_data function
