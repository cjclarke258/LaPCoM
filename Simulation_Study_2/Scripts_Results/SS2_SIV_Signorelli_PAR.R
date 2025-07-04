####################################################################################################
# SIMULATION STUDY 2 - SCENARIO IV - VIA SIGNORELLI
####################################################################################################
rm(list = ls())

setwd("/Users/cjclarke/Desktop/PhD/PhD Project 1/LaPCoM_Code")

source("Simulation_Study_2/Scripts/Seeds_Study_2.R")
source("LaPCoM_Generate_Data.R")

source("Simulation_Study_2/Scripts/Signorelli_Wit_Code/Signorelli_0-design-matrix-functions.R")
source("Simulation_Study_2/Scripts/Signorelli_Wit_Code/Signorelli_0-functions-starting-points.R")
source("Simulation_Study_2/Scripts/Signorelli_Wit_Code/Signorelli_0-EM-algorithm-glm.R")
source("Simulation_Study_2/Scripts/Signorelli_Wit_Code/Signorelli_0-glm-estimator.R")
source("Simulation_Study_2/Scripts/Signorelli_Wit_Code/extra_funcs_Signorelli.R")
source("Simulation_Study_2/Scripts/Signorelli_Wit_Code/fit_Signorelli.R")

# ------------------------------------------------------------------------------------------------    
# STORAGE SET UP
# ------------------------------------------------------------------------------------------------
output_SS2_SIV_Signorelli = list()

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

output_SS2_SIV_Signorelli = foreach::foreach(sim = 1:10) %dopar% {
  
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  # SET THE SEED
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  seed = seeds_study_2[31:40][sim]
  set.seed(seed)
  
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  # CREATE A LOG FILE
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  log_file = paste0("Simulation_Study_2/Results/SS2_SIV_Signorelli/SS2_SIV_Signorelli_Sim", sim, "_log.txt")
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
  res_Signorelli = fit_Signorelli(sim_multi$multi)
  
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  # STORAGE
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  output = list(sim_multi = sim_multi,
                res_Signorelli = res_Signorelli)
  
  write(paste0("\nSimulation ", sim, ": Finished!"), file = log_file, append = T)
  
  save(output, file = paste0("Simulation_Study_2/Results/SS2_SIV_Signorelli/Output_SS2_SIV_Signorelli_Simulation_", 
                             sim, ".Rdata"))
  
} # end foreach loop

parallel::stopCluster(cl)

tictoc::toc()