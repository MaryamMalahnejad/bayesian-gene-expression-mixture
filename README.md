# Bayesian Gaussian Mixture Model for Gene Expression (GSE45827)

## Project Overview

This project implements a Bayesian two-component Gaussian mixture model in Stan to identify latent expression subgroups in breast cancer gene expression data.

The objective is to demonstrate practical and theoretical understanding of:

- Bayesian statistical modeling  
- Finite mixture models  
- Markov Chain Monte Carlo (MCMC) inference  
- Probabilistic programming using Stan  
- Reproducible data science workflows in R  

---

## Dataset

The gene expression dataset used in this project is publicly available on Kaggle:

**Breast Cancer Gene Expression (CUMIDA)**  
https://www.kaggle.com/datasets/brunogrisci/breast-cancer-gene-expression-cumida  

The specific dataset used:

Breast_GSE45827.csv

Dataset characteristics:

- 151 samples  
- 54,677 gene probes  
- First columns: metadata (`samples`, `type`)  
- Remaining columns: gene expression values  

Due to GitHub's 100 MB file size limit, the dataset is not included in this repository.

### Reproducibility

To reproduce the analysis:

1. Download the dataset from Kaggle.
2. Extract the file.
3. Place it inside:
/data/raw
so that the file path becomes:
data/raw/Breast_GSE45827.csv

---

## Statistical Model

For a selected gene, standardized expression values  
x₁, x₂, …, xₙ are modeled as a two-component Gaussian mixture:

p(xₙ) = θ · Normal(μ₁, σ₁) + (1 − θ) · Normal(μ₂, σ₂)

where:

- μ₁, μ₂ are component means  
- σ₁, σ₂ > 0 are component standard deviations  
- θ ∈ (0,1) is the mixing proportion  

To prevent label switching, the component means are constrained:

μ₁ < μ₂

---

## Prior Distributions

The priors are defined as:

- μₖ ~ Normal(0, 5)  
- σₖ ~ Exponential(1)  
- θ ~ Beta(2, 2)  

Gene expression values are standardized before modeling:

x = (x_raw − mean(x_raw)) / sd(x_raw)

This improves numerical stability and makes priors weakly informative.

---

## Inference

Inference is performed using:

- Stan (via cmdstanr)
- 4 parallel MCMC chains
- 1000 warmup iterations
- 1000 sampling iterations per chain
- adapt_delta = 0.95

Convergence diagnostics:

- R-hat ≈ 1.00–1.01  
- High effective sample sizes  
- No divergent transitions  

Model evaluation:

- Leave-One-Out Cross-Validation (LOO)

---

## Example Results

Posterior mean estimates for an example probe:

- θ ≈ 0.44  
- μ₁ ≈ −0.72  
- μ₂ ≈ 0.94  
- σ₁ ≈ 0.61  
- σ₂ ≈ 0.50  

The separation between μ₁ and μ₂ suggests bimodal expression behavior for the selected gene.

Posterior predictive checks confirm that the model captures the observed distributional structure.

---

## Running the Project

From the project root directory:

```r
source("r/run.R")
Outputs are automatically generated in:
output/
```
## Generated files include:

Model fit (.rds)

Posterior summaries (.csv)

Diagnostics (R-hat, ESS)

Trace plots

Posterior density plots

Posterior predictive checks

LOO results

Run report
## Project Structure
project_mixture_geneexpr/
│
├── data/
│   └── raw/
│       └── Breast_GSE45827.csv
│
├── stan/
│   └── mixture_gaussian_2comp.stan
│
├── r/
│   └── run.R
│
├── output/
│
└── README.md
## Scientific Relevance

Finite mixture models are useful in genomics for detecting latent biological subgroups, tumor heterogeneity, or differential expression patterns.

This project demonstrates how Bayesian inference:

Quantifies uncertainty in latent structure

Avoids hard clustering

Provides full posterior distributions

Enables principled model comparison

Future Extensions

Three-component mixture model

Hierarchical Bayesian mixture

Automatic gene screening pipeline

Subtype-specific modeling

Comparison with EM-based Gaussian mixtures

Author

Maryam Malahnejad
MSc Applied Informatics
University of Duisburg-Essen
---