# NSCLC-microsim



# Project Structure
```
├── README.md
├── NSCLC microsim scripts/
│   ├── NSCLC microsim data/          # Regression model parameters and covariance matrices
│   ├── NSCLC microsim git prob.R     # Probabilistic microsimulation model 
│   ├── Parms_sample.R            	  # Samples parametric distribution parameters
│   ├── Probs_PSA.R              	  # Computes transition probabilities for 'j' second-order iteration of probabilistic model
│   ├── NSCLC microsim git.R          # deterministic microsimulation model
│   ├── Probs_determ.R                # Computes transition probabilities for deterministic model
│   ├── optputs/                 
│   ├── Save_outputs_prob.R  		  # Save results from 'NSCLC microsim git prob.R' to 'outputs/'
│   ├── Save_outputs_determ.R   	  # Save results from 'NSCLC microsim git.R' to 'outputs/'
│   ├── Determ_analysis    		      # Summary statistics and plots for the deterministic model
│   └── Prob_analysis                 # Summary statistics and plots for the probabilistic model
```

# Requirements
```
install.packages(c("parallel","flexsurv","survival", "MASS"))
```
