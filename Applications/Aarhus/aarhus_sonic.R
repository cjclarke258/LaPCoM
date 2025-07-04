###################################################################################################
# BABOON DATA APPLICATION (mixLPMC)
####################################################################################################

rm(list = ls())

# ==================================================================================================
# LOAD IN ALL NECESSARY FILES AND FUNCTIONS
# ==================================================================================================
# setwd("/Users/cjclarke/Desktop/PhD/PhD Project 1/LaPCoM_Code/Applications/baboons/")

source("LaPCoM_Initialisation_Functions.R")
source("LaPCoM_FCs.R")
source("LaPCoM.R")

# ================================================================================================
# LOAD AND FORMAT THE DATA
# ================================================================================================
aarhus_multi = readRDS("aarhus_adjacency_matrices_updated.rds")

# ------------------------------------------------------------------------------------------------    
# SET UP
# ------------------------------------------------------------------------------------------------
output_aarhus = list()

cols = khroma::color("light")(9)

# ------------------------------------------------------------------------------------------------
# GENERATE THE SEEDS
# ------------------------------------------------------------------------------------------------
set.seed(5678)
aarhus_seeds = sample(1:10000, 40)

# ------------------------------------------------------------------------------------------------
# SET THE HYPERPARAMETERS
# ------------------------------------------------------------------------------------------------
net_type = "binary"
net_mode = "undirected"

G_max = 5
s_e = 2

K_max = 8
s_w = 3
u_sigma = 21

delta_alpha = 2
alpha_prop_sd = 0.08

delta_Z = 0.07

# thin = 1
thin = 300
samples = 1000
show_plots = F
show_info = T

# ------------------------------------------------------------------------------------------------
# RUN THE 40 SIMULATIONS
# ------------------------------------------------------------------------------------------------
no_cores = 10
doParallel::registerDoParallel(cores = no_cores)
cl = parallel::makeCluster(no_cores, type = "FORK")

`%dopar%` = foreach::`%dopar%`

tictoc::tic()

output_aarhus = foreach::foreach(sim = 1:40) %dopar% {
  
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  # SET THE SEED
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  seed = aarhus_seeds[sim]
  set.seed(seed)
  
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  # CREATE A LOG FILE
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  log_file = paste0("aarhus_Logs/aarhus_", sim, ".txt")
  file.create(log_file)
  
  write(paste0("\n######################", "\n", "Starting Simulation: ", sim, "\n", 
               "Seed: ", seed, "\n",
               "######################", "\n"), file = log_file, append = T)
  
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  # FIT MODEL WITH NODE-LEVEL CLUSTERING
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  res_aarhus = LaPCoM(multi = aarhus_multi,
                     net_type = net_type,
                     net_mode = net_mode,
                     G_max = G_max,
                     init_method = "kmeans",
                     
                     add_noise_init_LS = T,
                     s_e = s_e,
                     delta_alpha = delta_alpha, 
                     alpha_prop_sd = alpha_prop_sd, 
                     delta_Z = delta_Z, 
                     
                     K_max = K_max, 
                     s_w = s_w, 
                     u_sigma = u_sigma,
                     
                     thin = thin, 
                     samples = samples, 
                     cols = cols, 
                     show_info = show_info, 
                     show_plots = show_plots,
                     log_file = log_file)
  
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  # STORAGE
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  output = list(aarhus_multi = aarhus_multi,
                res_aarhus = res_aarhus)
  
  write(paste0("\n Finished!"), file = log_file, append = T)
  
  save(output, file = paste0("aarhus_Output/aarhus_output_sim_", sim, ".Rdata"))
} # end foreach loop

parallel::stopCluster(cl)

tictoc::toc()