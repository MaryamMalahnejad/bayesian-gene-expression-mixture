# Bayesian Gaussian Mixture Model for Gene Expression (GSE45827)

## Project Overview

This project implements a Bayesian two-component Gaussian mixture model in Stan to identify latent expression subgroups in breast cancer gene expression data.

The objective is to demonstrate practical and theoretical understanding of:

- Bayesian statistical modeling  
- Finite mixture models  
- Markov Chain Monte Carlo (MCMC) inference  
- Probabilistic programming using Stan  
- Reproducible data science workflows in R  

The project is motivated by biological heterogeneity in cancer gene expression, where tumor samples may exhibit latent subgroups not directly observable from metadata.

---

## Dataset

The gene expression dataset used in this project is publicly available on Kaggle:

**Breast Cancer Gene Expression (CUMIDA)**  
https://www.kaggle.com/datasets/brunogrisci/breast-cancer-gene-expression-cumida  

Dataset used:

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

## Biological Motivation and Interpretation

For demonstration, we analyze a gene probe showing clear evidence of bimodal expression behavior across samples.

As an example, the gene ESR1 (Estrogen Receptor 1) is frequently studied in breast cancer research. ESR1 expression is strongly associated with hormone receptor status and plays a central role in tumor subtype classification (e.g., ER-positive vs ER-negative tumors).

Breast cancer is known to be biologically heterogeneous. Distinct tumor subtypes often exhibit different gene expression patterns. Modeling expression with a mixture distribution allows:

Identification of latent expression subgroups

Probabilistic classification of samples

Quantification of uncertainty

Avoidance of arbitrary thresholding

The two mixture components may correspond to biologically meaningful subpopulations, such as hormone-receptor-positive versus negative tumors.



## Statistical Model

For a selected gene, standardized expression values  
x₁, x₂, …, xₙ  
are modeled as a two-component Gaussian mixture:

p(xₙ) = θ · N(xₙ | μ₁, σ₁²) + (1 − θ) · N(xₙ | μ₂, σ₂²)

where:

- μ₁, μ₂ are component means  
- σ₁, σ₂ > 0 are component standard deviations  
- θ ∈ (0,1) is the mixing proportion  

To prevent label switching, the component means are constrained:

μ₁ < μ₂

The likelihood for independent samples is:

L(μ₁, μ₂, σ₁, σ₂, θ | x) = ∏_{n=1}^{N} p(xₙ)
---

## Why a Mixture Model?

Gene expression in cancer data is often heterogeneous.
A single Gaussian distribution may fail to capture multiple latent subpopulations.

A two-component mixture model allows:

Modeling latent biological subgroups

Capturing bimodal expression patterns

Full uncertainty quantification via posterior distributions

Interpretable probabilistic clustering

This aligns with modern probabilistic modeling approaches used in genomics.

## Prior Distributions

The priors are defined as:

- μₖ ~ Normal(0, 5)  
- σₖ ~ Exponential(1)  
- θ ~ Dirichlet(2, 2)

Gene expression values are standardized before modeling:

x = (x_raw − mean(x_raw)) / sd(x_raw)

Standardization improves numerical stability and ensures priors are weakly informative.

---
## Bayesian Workflow

The modeling process follows a principled Bayesian workflow:

Prior specification

Model implementation in Stan

MCMC sampling

Convergence diagnostics (R-hat, ESS)

Posterior predictive checks

Model comparison using LOO

Biological interpretation

This workflow ensures model validity and interpretability.

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
## Scientific Relevance

Finite mixture models are widely used in genomics to detect:

-Tumor heterogeneity

-Latent biological subgroups

-Differential expression structures

Hidden molecular subtypes

-Bayesian inference provides:

-Full posterior distributions

-Uncertainty quantification

-Principled model comparison

-Avoidance of hard clustering

This project is conceptually related to generalized linear mixture multilevel Bayesian models discussed in the course. While the present implementation is a simpler two-component Gaussian mixture, it illustrates the core principles of probabilistic modeling of biological heterogeneity.

## Limitations

-The model assumes exactly two components.

-No hierarchical structure across patients is modeled.

-Clinical covariates (tumor subtype, age) are not incorporated.

-Standardization removes absolute expression scale information.

Future extensions could include hierarchical Bayesian mixtures or multilevel modeling approaches.

## Future Extensions

-Three-component mixture model

-Hierarchical Bayesian mixture

-Automatic gene screening pipeline

-Subtype-specific modeling

-Comparison with EM-based Gaussian mixtures

## Running the Project

From the project root directory:

```r
source("r/run.R")
Outputs are automatically generated in:
output/
```
## Generated files include:

-Model fit (.rds)

-Posterior summaries (.csv)

-Diagnostics (R-hat, ESS)

-Trace plots

-Posterior density plots

-Posterior predictive checks

-LOO results

-Run report

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

# Author

Maryam Malahnejad
MSc Applied Informatics
University of Duisburg-Essen
---