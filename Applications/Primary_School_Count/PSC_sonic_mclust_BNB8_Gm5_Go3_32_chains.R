####################################################################################################
# PRIMARY_SCHOOL DATA APPLICATION (mixLPMC)
####################################################################################################

rm(list = ls())

# ==================================================================================================
# LOAD IN ALL NECESSARY FILES AND FUNCTIONS
# ==================================================================================================
# setwd("/Volumes/My Passport/PhD/PhD Project 1/Code/Applications/Primary_School_Count/")

source("mixLPCM_Initialisation_Functions.R")
source("mixLPCM_FCs.R")
source("mixLPCM_mclust.R")

# ================================================================================================
# LOAD AND FORMAT THE DATA
# ================================================================================================
primary_school_multi = readRDS("primary_school_adj_matrices_count.rds")

# ------------------------------------------------------------------------------------------------    
# SET UP
# ------------------------------------------------------------------------------------------------
output_primary_school = list()

cols = rep(khroma::color("light")(9), 2)

# ------------------------------------------------------------------------------------------------
# GENERATE THE SEEDS
# ------------------------------------------------------------------------------------------------
set.seed(56413)
primary_school_seeds = sample(1:10000, 32)

# ------------------------------------------------------------------------------------------------
# SET THE HYPERPARAMETERS
# ------------------------------------------------------------------------------------------------
net_type = "count"
net_mode = "undirected"

G_max = 5
G_0 = NULL

a_G = 8
b_G = 18
c_G = 10

s_e = 2

K_max = 15
s_w = 1
u_sigma = 41

delta_alpha = 2
alpha_prop_sd = 0.005

delta_Z = 0.004

# thin = 1
thin = 250
samples = 1000
show_plots = F
show_info = T

# ------------------------------------------------------------------------------------------------
# RUN THE 32 SIMULATIONS
# ------------------------------------------------------------------------------------------------
no_cores = 32
doParallel::registerDoParallel(cores = no_cores)
cl = parallel::makeCluster(no_cores, type = "FORK")

`%dopar%` = foreach::`%dopar%`

tictoc::tic()

output_primary_school = foreach::foreach(sim = 1:32) %dopar% {
  
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  # SET THE SEED
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  seed = primary_school_seeds[sim]
  set.seed(seed)
  
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  # CREATE A LOG FILE
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  log_file = paste0("PSC_mclust_BNB8_Gm5_Go3_32_chains_Logs/PSC_mclust_BNB8_Gm5_Go3_", sim, ".txt")
  file.create(log_file)
  
  write(paste0("\n######################", "\n", "Starting Simulation: ", sim, "\n", 
               "Seed: ", seed, "\n",
               "######################", "\n"), file = log_file, append = T)
  
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  # FIT MODEL WITH NODE-LEVEL CLUSTERING
  # ``````````````````````````````````````````````````````````````````````````````````````````````
  res_primary_school = mixLPCM(multi = primary_school_multi,
                         net_type = net_type,
                         net_mode = net_mode,

                         init_method = "mclust",
                         
                         G_max = G_max,
                         G_0 = G_0,
                         
                         a_G = a_G,
                         b_G = b_G,
                         c_G = c_G,
                         
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
  output = list(primary_school_multi = primary_school_multi,
                res_primary_school = res_primary_school)
  
  write(paste0("\n Finished!"), file = log_file, append = T)
  
  save(output, file = paste0("PSC_mclust_BNB8_Gm5_Go3_32_chains_Output/PSC_mclust_BNB8_Gm5_Go3_sim_", sim, ".Rdata"))
} # end foreach loop

parallel::stopCluster(cl)

tictoc::toc()