####################################################################################################
# SIMULATION STUDY 2 - SCENARIO II - VIA REBAFKA
####################################################################################################
rm(list = ls())

setwd("/Users/cjclarke/Desktop/PhD/PhD Project 1/LaPCoM_Code")

source("Simulation_Study_2/Scripts/Seeds_Study_2.R")
source("LaPCoM_Generate_Data.R")

# ------------------------------------------------------------------------------------------------    
# STORAGE SET UP
# ------------------------------------------------------------------------------------------------
output_SS2_SII_Rebafka = list()

# ------------------------------------------------------------------------------------------------
# SET THE HYPERPARAMETERS FOR DATA GENERATION
# ------------------------------------------------------------------------------------------------
net_type = "binary"
M = 50
N = 30
alpha = 0.6
G = 2
K = c(2, 3)
tau = c(0.6, 0.4)
pi = list(c(0.7, 0.3), c(0.4, 0.3, 0.3))

# ------------------------------------------------------------------------------------------------
# RUN THE 10 SIMULATIONS
# ------------------------------------------------------------------------------------------------
no_cores = 10 
doParallel::registerDoParallel(cores = no_cores)
cl <- parallel::makeCluster(no_cores, type = "FORK")

`%dopar%` <- foreach::`%dopar%`

tictoc::tic()

output_SS2_SII_Rebafka = foreach::foreach(sim = 1:10) %dopar% {
  
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  # SET THE SEED
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  seed = seeds_study_2[11:20][sim]
  set.seed(seed)
  
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  # CREATE A LOG FILE
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  log_file = paste0("Simulation_Study_2/Results/SS2_SII_Rebafka/SS2_SII_Rebafka_Sim", sim, "_log.txt")
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
  multi_list = lapply(seq(dim(sim_multi$multi)[3]), function(x) {sim_multi$multi[, , x]})

  res_Rebafka = list()
  start_time = Sys.time()
  res_Rebafka$res = graphclust::graphClustering(multi_list, nbCores = 1)
  end_time = Sys.time()
  res_Rebafka$time_taken = end_time - start_time
  
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  # STORAGE
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  output = list(sim_multi = sim_multi,
                res_Rebafka = res_Rebafka)
  
  write(paste0("\nSimulation ", sim, ": Finished!"), file = log_file, append = T)
  
  save(output, file = paste0("Simulation_Study_2/Results/SS2_SII_Rebafka/Output_SS2_SII_Rebafka_Simulation_", 
                             sim, ".Rdata"))
  
} # end foreach loop

parallel::stopCluster(cl)

tictoc::toc()