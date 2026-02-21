library(readr)
library(dplyr)
library(cmdstanr)
library(posterior)


df <- read_csv("data/raw/Breast_GSE45827.csv")


num_cols <- sapply(df, is.numeric)
X <- df[, num_cols]
vars <- sapply(X, var, na.rm = TRUE)
gene_name <- names(vars)[which.max(vars)]
x <- as.vector(scale(X[[gene_name]]))

stan_data <- list(N = length(x), x = x)
cat("Using gene column:", gene_name, "\n")


mod <- cmdstan_model("stan/mixture_gaussian_2comp.stan")


fit <- mod$sample(
  data = stan_data,
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  seed = 123
)


print(fit$summary(c("theta","mu1","mu2","sigma1","sigma2"))) 

dir.create("output", showWarnings = FALSE)
fit$save_object("output/fit_mixture.rds")

# ---- Save summary table ----
sum_tbl <- fit$summary(c("theta","mu1","mu2","sigma1","sigma2"))
write.csv(sum_tbl, "output/summary_params.csv", row.names = FALSE)

writeLines(gene_name, "output/selected_gene.txt")

png("output/hist_x.png", width = 900, height = 600)
hist(x, breaks = 40, main = paste("Standardized expression:", gene_name), xlab = "x (scaled)")
dev.off()

draws <- fit$draws(c("theta","mu1","mu2","sigma1","sigma2"))
draws_mat <- posterior::as_draws_matrix(draws)

png("output/posterior_theta.png", width = 900, height = 600)
hist(draws_mat[, "theta"], breaks = 40, main = "Posterior of theta", xlab = "theta")
dev.off()

png("output/posterior_mu1.png", width = 900, height = 600)
hist(draws_mat[, "mu1"], breaks = 40, main = "Posterior of mu1", xlab = "mu1")
dev.off()

png("output/posterior_mu2.png", width = 900, height = 600)
hist(draws_mat[, "mu2"], breaks = 40, main = "Posterior of mu2", xlab = "mu2")
dev.off()