# ==================================================================================================
# FUNCTIONS SPECIFIC TO THE WITHOUT_NODE_CLUSTERING FUNCTION
# ==================================================================================================

# ------------------------------------------------------------------------------------------------
# FUNCTION TO CALCULATE THE FULL CONDITIONAL FOR Z_g
# ------------------------------------------------------------------------------------------------
log_FC_Z_g_monoLaPCM = function(multi, Cvec, g, alpha, Z_g, D_g, beta, B, net_type, net_mode) {
  m_inds = which(Cvec == g)
  multi_subset = multi[, , m_inds]
  if (length(m_inds) == 1) {
    ll_multi = log_LPM_fast(multi_subset, alpha, D_g, net_type, net_mode)
  } else {
    ll_multi = sum(apply(multi_subset, 3, log_LPM_fast, alpha, D_g, net_type, net_mode))
  } # end if else statement
  
  ll_prior = sum( Rfast::dmvnorm(Z_g, beta, B, logged = T) )
  fc = ll_multi + ll_prior
  return(fc)
} # end log_FC_Z_g_monoLaPCM function 

# ------------------------------------------------------------------------------------------------
# FUNCTION TO CALCULATE THE DENSITY OF Z_g
# ------------------------------------------------------------------------------------------------
log_dens_Z_g_monoLaPCM = function(Z_g, Z_g_other, Sigma_Z) {
  ll = sum( dnorm(Z_g, Z_g_other, sqrt(Sigma_Z[1, 1]), log = T) )
  return(ll)
} # end log_dens_Z_g_NNC function