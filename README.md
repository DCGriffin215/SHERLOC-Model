# SHERLOC Model
This repository contains the code associated with the manuscript "Development and Validation of the SHERLOC Model: A Disease Model for Non-Small Cell Lung Cancer". Where the development and validation of the "Simulation for the Health Economics of Routine NSCLC Oncology in Canada" or SHERLOC model conatined in this repository is described and a case study is performed evaluating the cost effectiveness of the FLAURA2 and MARIPOSA trials, two trials comparing Osimertinib plus platinum pemetrexed and Amivantimab plus Lazertinib against Osimertinib monotherapy in the first line EGFR positive setting.

# Project Structure
```
├── README.md
├── NSCLC microsim scripts/
│   ├── SHERLOC microsim Het Data/          # Regression model parameters and covariance matrices
│   ├── SHERLOC microsim Het Data Alt/          # Regression model parameters and covariance matrices
│   ├── SHERLOC microsim git prob.R     # Probabilistic microsimulation model 
│   ├── Parms_sample_het.R            	  # Samples parametric distribution parameters
│   ├── Parms_sample_het_alt.R            	  # Samples parametric distribution parameters for alternate model
│   ├── 
│   ├── SHERLOC microsim git.R          # deterministic microsimulation model
│   ├── Probs_determ_het.R                # Computes transition probabilities for early stages
│   ├── Probs_determ_het_alt.R                # Computes transition probabilities for early stages in alternate model
│   ├── Probs_determ_het_late.R                # Computes transition probabilities for late stages
│   ├── Probs_determ_het_late_alt.R                # Computes transition probabilities for late stages in alternate model
│   ├── Outputs/                 
│   ├── Save_outputs_prob.R  		  # Save results from 'NSCLC microsim git prob.R' to 'outputs/'
│   ├── Save_outputs_determ.R   	  # Save results from 'NSCLC microsim git.R' to 'outputs/'
│   ├── Determ_analysis.R    		      # Summary statistics and plots for the deterministic model
│   └── Prob_analysis.R                 # Summary statistics and plots for the probabilistic model
```

# Requirements
```
install.packages(c("parallel","flexsurv","survival", "MASS","ggplot2","tidyverse","scales", "patchwork","truncnorm"))
```
