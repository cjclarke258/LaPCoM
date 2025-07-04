####################################################################################################
# SIMULATION STUDY 2 - SCENARIO V - VIA OUR MCMLPM
####################################################################################################
rm(list = ls())

# setwd("/Users/cjclarke/Desktop/PhD/PhD Project 1/LaPCoM_Code")

# source("Simulation_Study_2/Scripts_Results/Seeds_Study_2.R")
source("Seeds_Study_2.R")
source("LaPCoM_Generate_Data.R")
source("LaPCoM_Initialisation_Functions.R")
source("LaPCoM_FCs.R")
source("LaPCoM.R")

# ------------------------------------------------------------------------------------------------    
# STORAGE SET UP
# ------------------------------------------------------------------------------------------------
output_SS2_SV_CJ = list()

# ------------------------------------------------------------------------------------------------
# SET THE HYPERPARAMETERS FOR DATA GENERATION
# ------------------------------------------------------------------------------------------------
net_type = "binary"
net_mode = "undirected"

M = 100
N = 60
alpha = -0.4
G = 4
K = c(1, 2, 2, 3)
tau = c(0.3, 0.3, 0.2, 0.2)
pi = list(c(1), c(0.5, 0.5), c(0.7, 0.3), c(0.4, 0.3, 0.3))

# ------------------------------------------------------------------------------------------------
# SET THE HYPERPARAMETERS FOR INFERENCE
# ------------------------------------------------------------------------------------------------
G_max = 10
s_e = 4

K_max = 8
s_w = 3
u_sigma = 21

delta_alpha = 1
alpha_prop_sd = 0.01
delta_Z = 0.01

thin = 100
samples = 700
show_plots = F
show_info = T

# ------------------------------------------------------------------------------------------------
# RUN THE 10 SIMULATIONS
# ------------------------------------------------------------------------------------------------
no_cores = 10 
doParallel::registerDoParallel(cores = no_cores)
cl <- parallel::makeCluster(no_cores, type = "FORK")

`%dopar%` <- foreach::`%dopar%`

tictoc::tic()

output_SS2_SV_CJ = foreach::foreach(sim = 1:10) %dopar% {
  
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  # SET THE SEED
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  seed = seeds_study_2[41:50][sim]
  set.seed(seed)
  
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  # CREATE A LOG FILE
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  log_file = paste0("SS2_SV_CJ/SS2_SV_CJ_Sim", sim, "_log.txt")
  file.create(log_file)
  
  write(paste0("\n######################", "\n", "Starting Simulation: ", sim, "\n", 
               "Seed: ", seed, "\n",
               "######################", "\n"), file = log_file, append = T)
  
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  # GENERATE DATA
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  sim_multi = simulate_data(net_type, M, N, alpha, G, K, tau, pi, log_file)
  
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  # FIT MODEL WITH NODE-LEVEL CLUSTERING
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  res_CJ = LaPCoM(sim_multi$multi, 
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
                res_CJ = res_CJ)
  
  write(paste0("\nSimulation ", sim, ": Finished!"), file = log_file, append = T)
  
  save(output, file = paste0("SS2_SV_CJ/Output_SS2_SV_CJ_Simulation_", sim, ".Rdata"))
  
} # end foreach loop

parallel::stopCluster(cl)

tictoc::toc()