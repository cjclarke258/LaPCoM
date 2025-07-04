####################################################################################################
# SIMULATION STUDY 2 - SCENARIO IV - VIA MANTZIOU
####################################################################################################
rm(list = ls())

setwd("/Users/cjclarke/Desktop/PhD/PhD Project 1/LaPCoM_Code")

source("Simulation_Study_2/Scripts/Seeds_Study_2.R")
source("LaPCoM_Generate_Data.R")

source("Simulation_Study_2/Scripts/Mantziou_Code/Mantziou_MCMC_SFM.R")
source("Simulation_Study_2/Scripts/Mantziou_Code/Mantziou_vec_to_graph.R")
source("Simulation_Study_2/Scripts/Mantziou_Code/Mantziou_setup.R")

# ------------------------------------------------------------------------------------------------    
# STORAGE SET UP
# ------------------------------------------------------------------------------------------------
output_SS2_SIV_Mantziou = list()

# ------------------------------------------------------------------------------------------------
# SET THE HYPERPARAMETERS FOR DATA GENERATION
# ------------------------------------------------------------------------------------------------
net_type = "binary"
M = 50
N = 60
alpha = -0.4
G = 4
K = c(1, 2, 2, 3)
tau = c(0.3, 0.3, 0.2, 0.2)
pi = list(c(1), c(0.5, 0.5), c(0.7, 0.3), c(0.4, 0.3, 0.3))

# ------------------------------------------------------------------------------------------------
# RUN THE 10 SIMULATIONS
# ------------------------------------------------------------------------------------------------
no_cores = 10 
doParallel::registerDoParallel(cores = no_cores)
cl <- parallel::makeCluster(no_cores, type = "FORK")

`%dopar%` <- foreach::`%dopar%`

tictoc::tic()

output_SS2_SIV_Mantziou = foreach::foreach(sim = 1:10) %dopar% {
  
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  # SET THE SEED
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  seed = seeds_study_2[31:40][sim]
  set.seed(seed)
  
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  # CREATE A LOG FILE
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  log_file = paste0("Simulation_Study_2/Results/SS2_SIV_Mantziou/SS2_SIV_Mantziou_Sim", sim, "_log.txt")
  file.create(log_file)
  
  write(paste0("\n######################", "\n", "Starting Simulation: ", sim, "\n", 
               "Seed: ", seed, "\n",
               "######################", "\n"), file = log_file, append = T)
  
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  # GENERATE DATA
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  sim_multi = simulate_data(net_type, M, N, alpha, G, K, tau, pi, log_file)
  
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  # FIT MODEL
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  C_max = 10
  setup_params = setup_Mantziou(sim_multi$multi, C_max)
  res_Mantziou = MCMC_mix_err_sbm_SFM(C = C_max, B = 2,
                                      iter = 1000*100, burn_in = 0.3*1000*100, thin = 100,
                                      N_data = dim(sim_multi$multi)[3], 
                                      n_nodes = dim(sim_multi$multi)[1],
                                      G_data = setup_params$multi_list_ut, # list of vectors/upper triangle matrices
                                      prob_vec = Reduce("+",setup_params$multi_list_ut) / length(setup_params$multi_list_ut),
                                      a_0 = rep(0.5, setup_params$C),
                                      b_0 = rep(0.5, setup_params$C),
                                      c_0 = rep(0.5, setup_params$C),
                                      d_0 = rep(0.5, setup_params$C),
                                      e_0 = rep(0.5, setup_params$C),
                                      f_0 = rep(0.5, setup_params$C),
                                      chi_0 = rep(1/setup_params$C, setup_params$C),
                                      a_e = 1, b_e = 400,
                                      z_init = setup_params$z_init,
                                      c_init = setup_params$c_init,
                                      p_init = 0.01, q_init = 0.01,
                                      theta_init = 0.3,
                                      e0_init = 0.2,
                                      G_m_init = setup_params$G_m_init,
                                      perturb = rep(0.001, setup_params$C),
                                      rw_steps = rep(0.001, setup_params$C),
                                      rw_steps_e0 = 0.001)
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  # STORAGE
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  output = list(sim_multi = sim_multi,
                res_Mantziou = res_Mantziou)
  
  write(paste0("\nSimulation ", sim, ": Finished!"), file = log_file, append = T)
  
  save(output, file = paste0("Simulation_Study_2/Results/SS2_SIV_Mantziou/Output_SS2_SIV_Mantziou_Simulation_", sim, ".Rdata"))
  
} # end foreach loop

parallel::stopCluster(cl)

tictoc::toc()