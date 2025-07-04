####################################################################################################
# KRACKHARDT DATA APPLICATION (mixLPMC)
####################################################################################################

rm(list = ls())

# ==================================================================================================
# LOAD IN ALL NECESSARY FILES AND FUNCTIONS
# ==================================================================================================

source("LaPCoM_Initialisation_Functions.R")
source("LaPCoM_FCs.R")
source("LaPCoM.R")

# ================================================================================================
# LOAD AND FORMAT THE DATA
# ================================================================================================
load("Krackhardt_Multiplex_Directed.Rdata")

# ------------------------------------------------------------------------------------------------    
# SET UP
# ------------------------------------------------------------------------------------------------
output_krack = list()

cols = khroma::color("light")(9)

# ------------------------------------------------------------------------------------------------
# GENERATE THE SEEDS
# ------------------------------------------------------------------------------------------------
set.seed(1234)
krack_seeds = sample(1:10000, 40)

# ------------------------------------------------------------------------------------------------
# SET THE HYPERPARAMETERS
# ------------------------------------------------------------------------------------------------
net_type = "binary"

G_max = 5
s_e = 4

K_max = 6
s_w = 3
u_sigma = 11

delta_alpha = 7
alpha_prop_sd = 0.01
delta_Z = 0.04

thin = 300
samples = 1000
show_plots = F
show_info = T

# ------------------------------------------------------------------------------------------------
# RUN THE 20 SIMULATIONS
# ------------------------------------------------------------------------------------------------
no_cores = 10
doParallel::registerDoParallel(cores = no_cores)
cl = parallel::makeCluster(no_cores, type = "FORK")

`%dopar%` = foreach::`%dopar%`

tictoc::tic()

output_krack = foreach::foreach(sim = 1:40) %dopar% {
    
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  # SET THE SEED
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  seed = krack_seeds[sim]
  set.seed(seed)
  
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  # CREATE A LOG FILE
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  log_file = paste0("Krack_Logs/krack_", sim, ".txt")
  file.create(log_file)
  
  write(paste0("\n######################", "\n", "Starting Simulation: ", sim, "\n", 
               "Seed: ", seed, "\n",
               "######################", "\n"), file = log_file, append = T)
  
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  # FIT MODEL WITH NODE-LEVEL CLUSTERING
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  res_krack = LaPCoM(multi = krack_multi,
                     net_type = net_type,
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
                     show_plots = show_plots, 
                     show_info = show_info, 
                     log_file = log_file)
  
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  # STORAGE
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  output = list(krack_multi = krack_multi,
                res_krack = res_krack)
  
  write(paste0("\n Finished!"), file = log_file, append = T)
  
  save(output, file = paste0("Krack_Output/krack_output_sim_", sim, ".Rdata"))
} # end foreach loop

parallel::stopCluster(cl)

tictoc::toc()