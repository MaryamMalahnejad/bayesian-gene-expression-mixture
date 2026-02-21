rm(list = ls())

pkgs <- c("cmdstanr", "posterior", "bayesplot", "ggplot2", "readr", "dplyr", "loo")
to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)

library(cmdstanr)
library(posterior)
library(bayesplot)
library(ggplot2)
library(readr)
library(dplyr)
library(loo)

set.seed(123)

stan_file <- "stan/mixture_gaussian_2comp.stan"
data_file <- "data/raw/Breast_GSE45827.csv"
out_dir <- "output"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stopifnot(file.exists(stan_file))
stopifnot(file.exists(data_file))

dat <- readr::read_csv(data_file, show_col_types = FALSE)

gene_name <- "117_at"
if (!(gene_name %in% names(dat))) gene_name <- names(dat)[3]

x_raw <- suppressWarnings(as.numeric(dat[[gene_name]]))
x_raw <- x_raw[is.finite(x_raw)]
if (length(x_raw) < 10) stop("Too few numeric values in selected gene column.")

x <- as.numeric(scale(x_raw))

stan_data <- list(N = length(x), x = x)

mod <- cmdstan_model(stan_file)

fit <- mod$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  seed = 123,
  adapt_delta = 0.95,
  refresh = 200
)

saveRDS(fit, file.path(out_dir, paste0("fit_", gene_name, ".rds")))

summ <- fit$summary()
write.csv(summ, file.path(out_dir, paste0("summary_", gene_name, ".csv")), row.names = FALSE)

draws <- fit$draws()
diag_tbl <- posterior::summarise_draws(draws, "rhat", "ess_bulk", "ess_tail")
write.csv(diag_tbl, file.path(out_dir, paste0("diagnostics_", gene_name, ".csv")), row.names = FALSE)

pars_main <- c("mu[1]", "mu[2]", "sigma[1]", "sigma[2]", "theta[1]")

png(file.path(out_dir, paste0("trace_", gene_name, ".png")), width = 1400, height = 900)
print(bayesplot::mcmc_trace(draws, pars = pars_main))
dev.off()

png(file.path(out_dir, paste0("dens_", gene_name, ".png")), width = 1400, height = 900)
print(bayesplot::mcmc_dens(draws, pars = pars_main))
dev.off()

x_rep_arr <- fit$draws("x_rep")
x_rep_mat <- posterior::as_draws_matrix(x_rep_arr)

set.seed(123)
idx <- sample(seq_len(nrow(x_rep_mat)), size = min(50, nrow(x_rep_mat)))

png(file.path(out_dir, paste0("ppc_hist_", gene_name, ".png")), width = 1400, height = 900)
print(bayesplot::ppc_hist(y = x, yrep = x_rep_mat[idx, , drop = FALSE], binwidth = 0.25))
dev.off()

log_lik <- fit$draws("log_lik")
loo_res <- loo::loo(log_lik)
capture.output(loo_res, file = file.path(out_dir, paste0("loo_", gene_name, ".txt")))

report_lines <- c(
  paste0("Gene: ", gene_name),
  paste0("N: ", length(x)),
  "Standardization: scale(x_raw)",
  paste0("Stan file: ", stan_file),
  paste0("Data file: ", data_file),
  "",
  "LOO:",
  paste(capture.output(loo_res), collapse = "\n")
)

writeLines(report_lines, con = file.path(out_dir, paste0("run_report_", gene_name, ".txt")))

message("Done. Outputs saved in: ", out_dir)