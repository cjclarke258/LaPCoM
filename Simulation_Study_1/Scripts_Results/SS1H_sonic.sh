#!/bin/bash -l
#SBATCH --job-name=SS1H_CJ
# speficity number of nodes 
#SBATCH -N 1

# specify number of tasks/cores per node required
#SBATCH --ntasks-per-node 21

# specify the walltime e.g 20 mins
#SBATCH -t 330:00:00

# set to email at start,end and failed jobs
#SBATCH --mail-type=ALL
#SBATCH --mail-user=courtney.clarke@ucdconnect.ie

# run from current directory
cd $SLURM_SUBMIT_DIR
module load R/4.4.0
# module load gcc/9.3

# command to use
Rscript Simulation_Study_1_Scenario_H.R
