########################################################################################################################
# FULL CONDITIONAL FUNCTIONS
########################################################################################################################

# ======================================================================================================================
# GENERIC FUNCTIONS
# ======================================================================================================================
# --------------------------------------------------------------------------------------------------------------------
# FUNCTION TO CALCULATE THE LOG-LIKELIHOOD/SAMPLING DISTRIBUTION OF A NETWORK IN COMPONENT g
# --------------------------------------------------------------------------------------------------------------------
log_LPM = function(Y_m, alpha, D_g, net_type, net_mode) {
  eta = alpha - D_g
  if (net_mode == "undirected") {
    if (net_type == "count") {
      ll = sum( dpois(Y_m[upper.tri(Y_m)], lambda = exp(eta[upper.tri(eta)]), log = T) )
    } else if (net_type == "binary") {
      ll = sum( as.matrix(eta[upper.tri(eta)]) * Y_m[upper.tri(Y_m)] ) - (2 * sum( log( 1 + exp(eta[upper.tri(eta)]) ) ))
    } else {
      cat("Uh oh. There was a problem with the log_LPM function.")
    } # end net_type check
  } else if (net_mode == "directed") {
    if (net_type == "count") {
      ll = sum( dpois(Y_m, lambda = exp(eta), log = T) )
    } else if (net_type == "binary") {
      ll = sum( as.matrix(eta) * Y_m ) - (2 * sum( log( 1 + exp(eta) ) ))
    } else {
      cat("Uh oh. There was a problem with the log_LPM function.")
    } # end net_type check
  } else {
    cat("Uh oh. You must indicate whether the networks in your multiplex are 'undirected' or 'directed'.")
  } # end if else statement
  return(ll)
} # end log_LPM function

Rcpp::sourceCpp("LaPCoM_log_LPM_fast.cpp")

# --------------------------------------------------------------------------------------------------------------------
# FUNCTION TO CALCULATE THE LOG POSTERIOR OF G
# --------------------------------------------------------------------------------------------------------------------
log_post_comps_G = function(G, a, b, c, G_plus, e, M_comps) {
  term1 = ((lgamma(a + G - 1)) + (lbeta(a + b, G - 1 + c))) - ((lgamma(a)) + (lgamma(G)) + (lbeta(b, c)))
  term2 = ((G_plus * log(e)) + (lfactorial(G))) - ((G_plus * log(G)) + (lfactorial(G - G_plus)))
  term3 = 0
  for (g in 1:G_plus) {
    term3 = term3 + ((lgamma(M_comps[g] + (e / G))) - (lgamma(1 + (e / G))))
  } # end g for loop
  lpcG = term1 + term2 + term3
  return(lpcG)
} # end log_post_comps_G function

# --------------------------------------------------------------------------------------------------------------------
# FUNCTION TO CALCULATE THE LOG POSTERIOR OF e
# --------------------------------------------------------------------------------------------------------------------
log_post_e = function(e, l, r, G_plus, M_comps, G) {
  term1 = df(e, l, r, log = T)
  term2 = ((G_plus * log(e)) + (lgamma(e))) - (lgamma(sum(M_comps) + e))
  term3 = 0
  for (g in 1:G_plus) {
    term3 = term3 + ((lgamma(M_comps[g] + (e / G))) - (lgamma(1 + (e / G))))
  } # end g for loop
  lpe = term1 + term2 + term3
  return(lpe)
} # end log_post_e function

# --------------------------------------------------------------------------------------------------------------------
# FUNCTION TO CALCULATE THE LOG FULL CONDITIONAL OF ALPHA
# --------------------------------------------------------------------------------------------------------------------
log_FC_alpha = function(alpha, M, G_plus, Cmat, multi, D, m_alpha, s_alpha, net_type, net_mode) {
  ll_multi = 0
  for (m in 1:M) {
    for (g in 1:G_plus) {
      
      ll_term = log_LPM_fast(multi[, , m], alpha, D[, , g], net_type, net_mode)
      if (is.infinite(ll_term)) {
        ll_term = log_LPM(multi[, , m], alpha, D[, , g], net_type, net_mode)
      } # end if statement
      
      ll_multi = ll_multi + (Cmat[m, g] * ll_term)
      
    } # end g for loop
  } # end m for loop
  ll_prior = dnorm(alpha, m_alpha, s_alpha, log = T)
  fc = ll_multi + ll_prior
  return(fc)
} # end log_FC_alpha function

# ======================================================================================================================
# FUNCTIONS SPECIFIC TO THE LaPCoM FUNCTION
# ======================================================================================================================

# --------------------------------------------------------------------------------------------------------------------
# FUNCTION TO CALCULATE THE FULL CONDITIONAL FOR Z_g
# --------------------------------------------------------------------------------------------------------------------
log_FC_Z_g = function(multi, Cvec, g, alpha, Z_g, D_g, K_plus, Svec, mus_g, Sigmas_g, net_type, net_mode) {
  m_inds = which(Cvec == g)
  multi_subset = multi[, , m_inds]
  if (length(m_inds) == 0) {
    ll_multi = 0
  } else if (length(m_inds) == 1) {
    ll_multi = log_LPM_fast(multi_subset, alpha, D_g, net_type, net_mode)
    if (is.infinite(ll_multi)) {
      ll_multi = log_LPM(multi_subset, alpha, D_g, net_type, net_mode)
    } # end if statement
  } else {
    ll_vec = rep(NA, length(m_inds))
    for (m in 1:length(m_inds)) {
      ll_vec[m] = log_LPM_fast(multi_subset[, , m], alpha, D_g, net_type, net_mode)
      if (is.infinite(ll_vec[m])) {
        ll_vec[m] = log_LPM(multi_subset[, , m], alpha, D_g, net_type, net_mode)
      } # end if statement
    } # end m for loop
    ll_multi = sum(ll_vec)
  } # end if else statement
  
  ll_prior = 0
  for (k in 1:K_plus) {
    i_inds = which(Svec == k)
    ll_prior = ll_prior + sum( Rfast::dmvnorm(Z_g[i_inds, ], mus_g[k, ], Sigmas_g[, , k], logged = T) )
  } # end k loop
  fc = ll_multi + ll_prior
  return(fc)
} # end log_FC_Z_g function 

# --------------------------------------------------------------------------------------------------------------------
# FUNCTION TO CALCULATE THE DENSITY OF Z_g
# --------------------------------------------------------------------------------------------------------------------
log_dens_Z_g = function(Z_g, Z_g_other, K_plus, Svec, Sigma_g) {
  ll = 0
  for (k in 1:K_plus) {
    i_inds = which(Svec == k)
    for (q in 1:dim(Z_g)[2]) {
      ll = ll + sum( dnorm(Z_g[i_inds, q], Z_g_other[i_inds, q], sqrt(Sigma_g[q, q, k]), log = T) )
    } # end q for loop
  } # end k loop
  return(ll)
} # end log_dens_Z_g function

# --------------------------------------------------------------------------------------------------------------------
# FUNCTION TO CALCULATE THE LOG POSTERIOR OF K
# --------------------------------------------------------------------------------------------------------------------
log_post_comps_K = function(K_g, a, b, c, K_plus_g, w_g, N_g) {
  term1 = ((lgamma(a + K_g - 1)) + (lbeta(a + b, K_g - 1 + c))) - ((lgamma(a)) + (lgamma(K_g)) + (lbeta(b, c)))
  term2 = ((K_plus_g * log(w_g)) + (lfactorial(K_g))) - ((K_plus_g * log(K_g)) + (lfactorial(K_g - K_plus_g)))
  term3 = 0
  for (k in 1:K_plus_g) {
    term3 = term3 + ((lgamma((N_g[k]) + (w_g / K_g))) - (lgamma(1 + (w_g / K_g))))
  } # end k for loop
  lpcK = term1 + term2 + term3
  return(lpcK)
} # end log_post_comps_K function

# --------------------------------------------------------------------------------------------------------------------
# FUNCTION TO CALCULATE THE LOG POSTERIOR OF w
# --------------------------------------------------------------------------------------------------------------------
log_post_w = function(w_g, l, r, K_plus_g, N_g, K_g) {
  term1 = df(w_g, l, r, log = T)
  term2 = ((K_plus_g * log(w_g)) + (lgamma(w_g))) - ((lgamma(sum(N_g) + w_g))) 
  term3 = 0
  for (k in 1:K_plus_g) {
    term3 = term3 + ((lgamma(N_g[k] + (w_g / K_g))) - (lgamma(1 + (w_g / K_g))))
  } # end k for loop
  lpw = term1 + term2 + term3
  return(lpw)
} # end log_FC_w function

# --------------------------------------------------------------------------------------------------------------------
# FUNCTION TO CALCULATE THE LOG-POSTERIOR OF LaPCoM
# --------------------------------------------------------------------------------------------------------------------
log_post = function(multi, alpha, Z, C_vector, net_type, net_mode,
                    G_plus, tau, e, G, a_G, b_G, c_G, l_G, r_G, 
                    K_plus, K, S_vector, pi, mu, Sigma, mu_mean, mu_cov, u_sigma, v_sigma, w, a_K, b_K, c_K, l_K, r_K) {
  
  safe_log = function(x, name) {
    if (is.na(x) || x <= 0) {
      stop(paste("Error: Invalid value for", name, "- value:", x))
    } else {
      return(log(x))
    } # end if statement
  } # end safe_log function
  
  # log-likelihood
  M = dim(multi)[3]
  log_like_term = 0
  for (m in 1:M) {
    log_like_term = log_LPM_fast(Y_m = multi[, , m], alpha = alpha, D = Rfast::Dist(Z[[C_vector[m]]], square = TRUE), 
                                 net_type = net_type, net_mode = net_mode)
  } # end m for loop
  rm(m)
  
  # C
  C_log_term = sum(sapply(1:G_plus, function(g) {
    sum(sapply(1:M, function(m) {
      if (C_vector[m] == g) {
        return(safe_log(tau[g], "tau[g]"))
      } # end if statement
      return(0)
    }))
  }))
  
  # tau
  tau_log_term = sum(sapply(1:G_plus, function(g) {
    return(((e / G) - 1) * safe_log(tau[g], "tau term"))
  }))
  
  # G
  G_log_term = (lgamma(a_G + G_plus - 1) + lbeta(a_G + b_G, G_plus - 1 + c_G)) - 
    (lgamma(a_G) + lgamma(G_plus) + lbeta(b_G, c_G))
  
  # e
  e_log_term = df(e, l_G, r_G, log = TRUE)
  
  # Z
  Z_log_term = 0
  N = dim(Z[[1]])[1]
  for (g in 1:G_plus) {
    
    if (K_plus[g] == 1 && K[g] == 1) {
      # Case: Both K_g_plus and K_g are 1
      Zg_log_term = sum(sapply(1:N, function(i) {
        if (S_vector[[g]][i] == 1) {
          return(safe_log(pi[[g]][1], "pi[1]") + mvtnorm::dmvnorm(Z[[g]][i, ], mu[[g]], diag(Sigma[[g]]), log = TRUE))
        } # end if statement
      }))
    } else if (K_plus[g] == 1 && K[g] != 1) {
      # Case: K_g_plus is 1, but K_g is not 1 
      Zg_log_term = sum(sapply(1:N, function(i) {
        if (S_vector[[g]][i] == 1) {
          return(safe_log(pi[[g]][1], "pi[1]") + mvtnorm::dmvnorm(Z[[g]][i, ], mu[[g]], diag(Sigma[[g]]), log = TRUE))
        } # end if statement
      }))
    } else {
      # Case: K_plus > 1 and K > 1 (general case)
      Zg_log_term = sum(sapply(1:N, function(i) { 
        sum(sapply(1:K_plus[g], function(k) { 
          if (S_vector[[g]][i] == k) { 
            return(safe_log(pi[[g]][k], "pi[k]") + 
                     mvtnorm::dmvnorm(Z[[g]][i, ], mu[[g]][k, ], diag(Sigma[[g]][k, ]), log = TRUE)) 
          } # end if statement
          return(0)  # Ensure a return value for all cases
        })) 
      }))
    } # end if else statement
    
    Z_log_term = Z_log_term + Zg_log_term
    
  } # end g for loop
  rm(g)
  
  # mu
  mu_log_term = 0
  for (g in 1:G_plus) {
    for (k in 1:K_plus[g]) {
      if (K_plus[g] == 1) {
        mu_log_term = mu_log_term + mvtnorm::dmvnorm(mu[[g]], mu_mean, mu_cov, log = T)
      } else {
        mu_log_term = mu_log_term + mvtnorm::dmvnorm(mu[[g]][k, ], mu_mean, mu_cov, log = T)
      } # end if else statement
    } # end k for loop
    rm(k)
  } # end g for loop
  rm(g)
  
  Sigma_log_term = 0
  for (g in 1:G_plus) {
    for (k in 1:K_plus[g]) {
      if (K_plus[g] == 1) {
        Sigma_log_term = Sigma_log_term + 
          sum(sapply(1:2, function(q) { LaplacesDemon::dinvgamma(Sigma[[g]][q], u_sigma, v_sigma, log = TRUE) }))
      } else {
        Sigma_log_term = Sigma_log_term + 
          sum(sapply(1:2, function(q) { LaplacesDemon::dinvgamma(Sigma[[g]][k, q], u_sigma, v_sigma, log = TRUE) }))
      } # end if else statement
    } # end k for loop
    rm(k)
  } # end g for loop
  rm(g)
  
  # S
  S_log_term = 0
  for (g in 1:G_plus) {
    Sg_log_term = sum(sapply(1:K_plus[g], function(k) {
      sum(sapply(1:N, function(i) {
        if (S_vector[[g]][i] == k) {
          return(safe_log(pi[[g]][k], "pi[k]"))
        } # end if statement
        return(0)
      }))
    }))
    S_log_term = S_log_term + Sg_log_term
  } # end g for loop
  rm(g)
  
  # pi
  pi_log_term = 0
  for (g in 1:G_plus) {
    pi_g_log_term = sum(sapply(1:K_plus[g], function(k) {
      return(((w[g] / G) - 1) * safe_log(pi[[g]][k], "pi term"))
    }))
    pi_log_term = pi_log_term + pi_g_log_term
  } # end g for loop
  rm(g)
  
  # K
  K_log_term = 0
  for (g in 1:G_plus) {
    Kg_log_term = (lgamma(a_K + K_plus[g] - 1) + lbeta(a_K + b_K, K_plus[g] - 1 + c_K)) - 
      (lgamma(a_K) + lgamma(K_plus[g]) + lbeta(b_K, c_K))
    K_log_term = K_log_term + Kg_log_term
  } # end g for loop
  rm(g)
  
  # w
  w_log_term = 0
  for (g in 1:G_plus) {
    w_log_term = w_log_term + df(w[g], l_K, r_K, log = TRUE)
  } # end g for loop
  rm(g)
  
  # alpha
  alpha_log_term = dnorm(alpha, m_alpha, s_alpha, log = T)
  
  lp = log_like_term + C_log_term + tau_log_term + G_log_term + e_log_term + Z_log_term + mu_log_term + 
        Sigma_log_term + S_log_term + pi_log_term + K_log_term + w_log_term + alpha_log_term
  
  return(lp)
  
} # end log_post function

# ======================================================================================================================
# COVARIANCE ELLIPSE FUNCTION
# ======================================================================================================================
get_ellipse = function(mean, cov_matrix, points = 100) {
  
  angles = seq(0, 2 * base::pi, length.out = points)
  ellipse_points = t(chol(cov_matrix)) %*% rbind(cos(angles), sin(angles))
  
  ellipse_df = data.frame(
    x = ellipse_points[1, ] + mean[1],
    y = ellipse_points[2, ] + mean[2])
  
  return(ellipse_df)
  
} # end get_ellipse function