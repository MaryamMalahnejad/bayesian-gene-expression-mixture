# Bayesian Gaussian Mixture Model for Gene Expression (GSE45827)

## Project Overview

This project implements a Bayesian two-component Gaussian mixture model in Stan to identify latent subgroups in breast cancer gene expression data.

The objective is to demonstrate proficiency in:

- Bayesian statistical modeling  
- Finite mixture models  
- MCMC-based inference  
- Implementation using R and Stan  

---

## Dataset

The gene expression dataset used in this project is publicly available on Kaggle:

**Breast Cancer Gene Expression (CUMIDA)**  
🔗 https://www.kaggle.com/datasets/brunogrisci/breast-cancer-gene-expression-cumida  

The specific dataset used is:

`Breast_GSE45827.csv`

Due to GitHub's 100 MB file size limit, the dataset is not included in this repository.

### How to Reproduce with the Dataset

1. Download the dataset from the Kaggle link above.
2. Place the file inside:
data/raw/
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
These results demonstrate the model’s ability to identify latent structure in high-dimensional biological data using Bayesian inference.
