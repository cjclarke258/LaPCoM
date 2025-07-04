####################################################################################################
# SIMULATION STUDY 1 SCENARIO C
####################################################################################################

rm(list = ls())
# setwd("/Users/cjclarke/Desktop/PhD/PhD Project 1/LaPCoM_Code")

# ==================================================================================================
# LOAD IN ALL NECESSARY FILES AND FUNCTIONS FOR SIMULATION STUDY 1 SCENARIO C
# ==================================================================================================
source("Seeds_Study_1.R")
source("LaPCoM_Generate_Data.R")
source("LaPCoM_Initialisation_Functions.R")
source("LaPCoM_FCs.R")
source("monoLaPCM_FCs.R")
source("LaPCoM.R")
source("monoLaPCM.R")

# ==================================================================================================
# SCENARIO C
# ==================================================================================================

# ------------------------------------------------------------------------------------------------    
# STORAGE SET UP
# ------------------------------------------------------------------------------------------------
output_SS1C = list()

# ------------------------------------------------------------------------------------------------
# SET THE HYPERPARAMETERS FOR DATA GENERATION
# ------------------------------------------------------------------------------------------------
net_type = "count"
net_mode = "undirected"

M = 20
N = 30
alpha = 0.6
G = 2
K = c(1, 2)
tau = c(6/10, 4/10)
pi = list(c(1), c(0.4, 0.6))

# ------------------------------------------------------------------------------------------------
# SET THE HYPERPARAMETERS FOR INFERENCE
# ------------------------------------------------------------------------------------------------
G_max = 5
s_e = 4

K_max = 8
s_w = 3
u_sigma = 11

delta_alpha = 3
alpha_prop_sd = 0.01
delta_Z = 0.008

thin = 300
samples = 1000
show_plots = F
show_info = T

# ------------------------------------------------------------------------------------------------
# RUN THE 10 SIMULATIONS
# ------------------------------------------------------------------------------------------------
no_cores = 10
# no_cores = parallel::detectCores() - 4
doParallel::registerDoParallel(cores = no_cores)
cl = parallel::makeCluster(no_cores, type = "FORK")

`%dopar%` <- foreach::`%dopar%`

tictoc::tic()

output_SS1C = foreach::foreach(sim = 1:10) %dopar% {
  
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  # SET THE SEED
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  seed = seeds_study_1[21:30][sim]
  set.seed(seed)
  
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  # CREATE A LOG FILE
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  log_file = paste0("SS1C_Logs/Output_SS1C_Simulation", sim, ".txt")
  file.create(log_file)
  
  write(paste0("\n######################", "\n", "Starting Simulation: ", sim, "\n", 
               "Seed: ", seed, "\n",
               "######################", "\n"), file = log_file, append = T)
  
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  # GENERATE DATA
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  sim_multi = simulate_data(net_type, M, N, alpha, G, K, tau, pi, log_file)
  
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  # FIT MODEL WITH NO NODE-LEVEL CLUSTERING
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  res_without_node_clusters = monoLaPCM(multi = sim_multi$multi,
                                     net_type = net_type,
                                     net_mode = net_mode,
                                     G_max = G_max,
                                     s_e = s_e,
                                     delta_alpha = delta_alpha,
                                     alpha_prop_sd = alpha_prop_sd,
                                     delta_Z = delta_Z,
                                     thin = thin, samples = samples,
                                     show_plots = show_plots,
                                     show_info = show_info, 
                                     log_file = log_file)
  
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  # FIT MODEL WITH NODE-LEVEL CLUSTERING
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  res_with_node_clusters = LaPCoM(multi = sim_multi$multi, 
                                  net_type = net_type,
                                  net_mode = net_mode,
                                  G_max = G_max,
                                  init_method = "kmeans",
                                  add_noise_init_LS = F,
                                  s_e = s_e,
                                  delta_alpha = delta_alpha,
                                  alpha_prop_sd = alpha_prop_sd,
                                  delta_Z = delta_Z,
                                  K_max = K_max,
                                  s_w = s_w,
                                  u_sigma = u_sigma,
                                  thin = thin, samples = samples,
                                  show_plots = show_plots, 
                                  show_info = show_info, 
                                  log_file = log_file)
  
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  # STORAGE
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  output = list(sim_multi = sim_multi,
                res_without_node_clusters = res_without_node_clusters,
                res_with_node_clusters = res_with_node_clusters)
  
  write(paste0("\nSimulation ", sim, ": Finished!"), file = log_file, append = T)
  
  save(output, file = paste0("SS1C_Output/Output_SS1C_Simulation_", sim, ".Rdata"))
} # end foreach loop

parallel::stopCluster(cl)

tictoc::toc()