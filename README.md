# LaPCoM

Multiplex networks are increasingly common across diverse domains, motivating the development of clustering methods that uncover patterns at multiple levels. Existing approaches typically focus on clustering either entire networks or nodes within a single network. We address the lack of a unified latent space framework for simultaneous network- and node-level clustering by proposing a latent position co-clustering model (LaPCoM), based on a hierarchical mixture-of-mixtures formulation. LaPCoM enables co-clustering of networks and their constituent nodes, providing joint dimension reduction and two-level cluster detection. At the network level, it identifies global homogeneity in topological patterns by grouping networks that share similar latent representations. At the node level, it captures local connectivity and community patterns. The model adopts a Bayesian nonparametric framework using a mixture of finite mixtures, which places priors on the number of clusters at both levels and incorporates sparse priors to encourage parsimonious clustering. Inference is performed via Markov chain Monte Carlo with automatic selection of the number of clusters. LaPCoM accommodates both binary and count-valued multiplex data. Simulation studies and comparisons with existing methods demonstrate accurate recovery of latent structure and clusters. Applications to real-world social multiplexes reveal interpretable network-level clusters aligned with context-specific patterns, and node-level clusters reflecting social patterns and roles.

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
