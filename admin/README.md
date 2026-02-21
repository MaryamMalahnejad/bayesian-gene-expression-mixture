# Bayesian Gaussian Mixture Model for Gene Expression (GSE45827)

## Project Overview

This project implements a Bayesian two-component Gaussian mixture model in Stan to identify latent subgroups in gene expression data from the dataset **Breast_GSE45827.csv**.

The goal is to demonstrate mastery of Bayesian modeling, mixture models, and implementation using R and Stan.

---

## Data

Dataset: Breast_GSE45827.csv  
Location: data/raw/

From all numeric gene expression columns, the gene with the highest variance across samples was selected.  
The expression values were standardized (mean 0, standard deviation 1).

---

## Model

A two-component Gaussian mixture model was fitted:

x_n ~ theta * Normal(mu1, sigma1) + (1 - theta) * Normal(mu2, sigma2)

### Priors

- mu_k ~ Normal(0, 5)  
- sigma_k ~ Exponential(1)  
- theta ~ Beta(2, 2)

---

## Inference

- Implemented in Stan using cmdstanr  
- 4 MCMC chains  
- 1000 warmup iterations  
- 1000 sampling iterations per chain  
- Convergence diagnostics: Rhat ≈ 1.01 for all parameters  

---

## Results (Posterior Means)

- theta ≈ 0.44  
- mu1 ≈ 0.94  
- mu2 ≈ -0.72  
- sigma1 ≈ 0.50  
- sigma2 ≈ 0.61  

The model identifies two latent expression clusters in the selected gene.

---

## Reproducibility

To reproduce the results:

```r
source("r/run.R")
```

All outputs (model fit, summary tables, plots) are saved in:

```
output/
```

---

## Project Structure

```
data/raw/       -> Original dataset
stan/           -> Stan model file
r/run.R         -> Main execution script
output/         -> Results and plots
admin/README.md -> Project documentation
```