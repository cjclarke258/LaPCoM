# LaPCoM

Multiplexes are increasingly prevalent across various domains, motivating the development of clustering methodologies to uncover structural patterns. Existing methodologies cluster either entire networks — e.g., using latent position models or stochastic blockmodels — or nodes within a single network, such as the latent position cluster model (LPCM). We address the absence of a unified latent- space framework for simultaneous network- and node-level clustering by proposing LaPCoM, a mixture-of-mixtures model. This framework enables *co-clustering* of networks and their constituent nodes, providing both dimension reduction and two-level cluster detection. At the network level, networks are clustered based on a shared underlying latent space; at the node level, clustering arises naturally within each network cluster via the LPCM structure. Inference and automatic selection of the number of clusters at both levels are carried out via a Markov chain Monte Carlo algorithm. LaPCoM accommodates multiplex data with binary or weighted edge values. Comparisons with existing methods demonstrate good performance in recovering underlying latent spaces and network- and node- level clusters. Applications to three real-world social multiplexes reveal cohesive network-level clusters reflecting temporal patterns and node-level clusters capturing underlying social roles.

## Project Structure

- **Code**: Scripts necessary to implement LaPCoM.
  - `LaPCoM.R`: Main function to implement LaPCoM.
  - `LaPCoM_log_LPM_fast.cpp`, `LaPCoM_Initialisation_Functions.R` and `LaPCoM_FCs.R`: Functions to calculate the log-likelihood of the LPM, initialise the *LaPCoM* parameters and calculate the full conditionals in the update steps within *LaPCoM*; all required to run `LaPCoM.R`.
  - `LaPCoM_Generate_Data.R`: Function to generate data from *LaPCoM*; used in the simulation studies.
  - `LaPCoM_PP_MethodB.R`: Function to carry out the post-processing procedure for *LaPCoM*.
  - `monoLaPCM_*`: Functions required to use the single mixture version of the model for comparison purposes in our first simulation study.
- **Simulation_Study_1/**
  - `Scripts_Results/`: Scripts for running simulation scenarios.
  - `Post_Processing/`: Scripts for post-processing MCMC output.
  - `Results_Plots/`: Scripts for summarising study results.
- **Simulation_Study_2/**
  - Similar structure to Simulation_Study_1.
- **Applications/**
  - `Krackhardt/`, `Aarhus/`, and `Primary_School_Count/`: Applications of LaPCoM to illustrative datasets.

---

## How to Use

1. Clone the repository.
2. Change all working directories and path files to be your own.
3. Main model functions are in the `Code/` directory.
4. Simulation studies and applications can be run from their respective folders.

## Requirements

- R (>= 4.0.0)

*(See individual scripts for specific dependencies and ideal organisation of files.)*

## License

This project is licensed under the GNU License — see the [LICENSE](LICENSE) file for details.
