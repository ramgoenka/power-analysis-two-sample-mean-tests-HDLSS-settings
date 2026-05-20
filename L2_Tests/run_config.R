###############################################################################
# run_config.R
#
# Tests included (8 total):
#
#   ORIGINAL NORMAL-APPROXIMATION L2 TESTS:
#     BS_norm : Bai & Saranadasa 1996 via HDNRA::BS1996.TS.NABT
#     CQ_norm : Chen & Qin 2010 via HDNRA::CQ2010.TSBF.NABT
#     SD_norm : Srivastava & Du 2008 via HDNRA::SD2008.TS.NABT
#     SKK_norm : Srivastava, Katayama, Kano 2013 via HDNRA::SKK2013.TSBF.NABT
#
#   SPARSITY-AWARE L2 TEST:
#     CLZ : Chen, Li, Zhong 2014 / 2019 via highmean::apval_Chen2014
#
#   NORMAL-REFERENCE VARIANTS:
#     BS_3c : Zhang & Zhu 2022 (TS, 3-cumulant) via HDNRA::ZZ2022.TS.3cNRT
#     CQ_3c : Zhang & Zhu 2022 (TSBF, 3-cumulant) via HDNRA::ZZ2022.TSBF.3cNRT
#     ZZZ23 : Zhang et al. 2023 (TSBF, scale-inv) via HDNRA::ZZZ2023.TSBF.2cNRT
#
# Grid:
#   p : {200, 500, 1000, 2000, 10000}                         
#   n : {5, 10, 20, 50, 100}
#   cov : {identity, cs, ar1}
#   signal: {null, sparse, moderate, dense}
#   rho : {0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.8}      (cs/ar1 only; identity NA)
###############################################################################
lib_path <- Sys.getenv("R_LIBS_USER")
if (lib_path == "") lib_path <- file.path(Sys.getenv("HOME"), "R_libs")
.libPaths(c(lib_path, .libPaths()))

suppressPackageStartupMessages({
  library(HDNRA)
  library(highmean)
})

###############################################################################
# Mean-shift constructors. All three alternative signals are calibrated so
# that ||mu2||^2 = 0.0125 * p, matching the original paper.
###############################################################################
make_mu_null     <- function(p) rep(0, p)
make_mu_sparse   <- function(p) { mu <- rep(0, p); mu[1:floor(0.05 * p)] <- 0.5;  mu }
make_mu_moderate <- function(p) { mu <- rep(0, p); mu[1:floor(0.20 * p)] <- 0.25; mu }
make_mu_dense    <- function(p) rep(sqrt(0.0125), p)

###############################################################################
# Efficient covariance-structured generators. All operate at O(np) cost --
# no p-by-p matrix is ever formed or factored explicitly.
###############################################################################
generate_identity <- function(n, p, mu = rep(0, p)) {
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  if (any(mu != 0)) X <- sweep(X, 2, mu, "+")
  X
}

generate_cs <- function(n, p, rho, mu = rep(0, p)) {
  # Compound symmetry: shared common factor across columns.
  # Eigenstructure: Sigma = (1-rho)*I + rho*J has rank-1 perturbation,
  # so X = sqrt(1-rho)*Z + sqrt(rho)*W where W is broadcast across columns.
  Z <- matrix(rnorm(n * p), nrow = n, ncol = p)
  W <- rnorm(n)
  X <- sqrt(1 - rho) * Z + sqrt(rho) * W
  if (any(mu != 0)) X <- sweep(X, 2, mu, "+")
  X
}

generate_ar1 <- function(n, p, rho, mu = rep(0, p)) {
  # AR(1): explicit Cholesky recurrence on the tridiagonal precision matrix.
  # Equivalent to applying the Cholesky factor of Sigma = (rho^|i-j|) at O(np).
  Z <- matrix(rnorm(n * p), nrow = n, ncol = p)
  scale <- sqrt(1 - rho^2)
  for (k in 2:p) Z[, k] <- rho * Z[, k - 1] + scale * Z[, k]
  if (any(mu != 0)) Z <- sweep(Z, 2, mu, "+")
  Z
}

generate_two_sample <- function(n, p, cov_model, signal, rho) {
  mu2 <- switch(signal,
                "null"     = make_mu_null(p),
                "sparse"   = make_mu_sparse(p),
                "moderate" = make_mu_moderate(p),
                "dense"    = make_mu_dense(p),
                stop("Unknown signal type: ", signal))
  
  if (cov_model == "identity") {
    X <- generate_identity(n, p)
    Y <- generate_identity(n, p, mu = mu2)
  } else if (cov_model == "cs") {
    X <- generate_cs(n, p, rho = rho)
    Y <- generate_cs(n, p, rho = rho, mu = mu2)
  } else if (cov_model == "ar1") {
    X <- generate_ar1(n, p, rho = rho)
    Y <- generate_ar1(n, p, rho = rho, mu = mu2)
  } else {
    stop("Unknown cov_model: ", cov_model)
  }
  list(X = X, Y = Y)
}

###############################################################################
# Modular test runner. Returns a named numeric vector of p-values, one per
# test. NA for any test that errors out (logged but does not halt simulation).
###############################################################################
run_all_tests <- function(X, Y) {
  safe <- function(expr) tryCatch(expr, error = function(e) NA_real_)
  
  c(
    BS_norm  = safe(BS1996.TS.NABT(X, Y)$p.value),
    CQ_norm  = safe(CQ2010.TSBF.NABT(X, Y)$p.value),
    SD_norm  = safe(SD2008.TS.NABT(X, Y)$p.value),
    SKK_norm = safe(SKK2013.TSBF.NABT(X, Y)$p.value),
    BS_3c    = safe(ZZ2022.TS.3cNRT(X, Y)$p.value),
    CQ_3c    = safe(ZZ2022.TSBF.3cNRT(X, Y)$p.value),
    ZZZ23    = safe(ZZZ2023.TSBF.2cNRT(X, Y)$p.value),
    CLZ      = safe(apval_Chen2014(X, Y)$pval)
  )
}

###############################################################################
# Build the simulation grid.
###############################################################################
p_grid   <- c(200, 500, 1000)
n_grid   <- c(5, 10, 20, 50, 100)
sig_grid <- c("null", "sparse", "moderate", "dense")
rho_grid <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.8)

n_iter <- 5000
alpha  <- 0.05

# Identity has no rho (encoded as NA, single config per p/n/signal).
grid_id  <- expand.grid(
  p = p_grid, n = n_grid, cov_model = "identity",
  signal = sig_grid, rho = NA_real_,
  stringsAsFactors = FALSE
)
# CS and AR(1) sweep across the rho grid.
grid_cor <- expand.grid(
  p = p_grid, n = n_grid, cov_model = c("cs", "ar1"),
  signal = sig_grid, rho = rho_grid,
  stringsAsFactors = FALSE
)
grid <- rbind(grid_id, grid_cor)

###############################################################################
# Read SLURM_ARRAY_TASK_ID and pull this task's configuration.
###############################################################################
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (is.na(task_id) || task_id < 1 || task_id > nrow(grid)) {
  stop(sprintf("Invalid SLURM_ARRAY_TASK_ID: %s (grid has %d rows)",
               Sys.getenv("SLURM_ARRAY_TASK_ID"), nrow(grid)))
}

cfg <- grid[task_id, ]
cat(sprintf("Task %d / %d:\n", task_id, nrow(grid)))
cat(sprintf("  p=%d, n=%d, cov=%s, signal=%s, rho=%s, B=%d\n",
            cfg$p, cfg$n, cfg$cov_model, cfg$signal,
            ifelse(is.na(cfg$rho), "NA", as.character(cfg$rho)), n_iter))

###############################################################################
# Set seed for reproducibility. set.seed(2026 + task_id) ensures that
# re-running a specific task_id gives identical output.
###############################################################################
set.seed(2026 + task_id)

test_names    <- c("BS_norm", "CQ_norm", "SD_norm", "SKK_norm",
                   "BS_3c",   "CQ_3c",   "ZZZ23",   "CLZ")
reject_counts <- setNames(integer(length(test_names)), test_names)
na_counts     <- setNames(integer(length(test_names)), test_names)

t_start <- Sys.time()

for (k in 1:n_iter) {
  dat <- generate_two_sample(cfg$n, cfg$p,
                             cov_model = cfg$cov_model,
                             signal    = cfg$signal,
                             rho       = cfg$rho)
  
  pvals <- run_all_tests(dat$X, dat$Y)
  
  # Count rejections (NA does NOT count as rejection)
  rejects <- !is.na(pvals) & pvals < alpha
  reject_counts <- reject_counts + rejects
  na_counts     <- na_counts + is.na(pvals)
  
  if (k %% 500 == 0) {
    elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
    cat(sprintf("  iter %d/%d  (%.1f min elapsed)\n", k, n_iter, elapsed))
  }
}

t_end   <- Sys.time()
elapsed <- as.numeric(difftime(t_end, t_start, units = "mins"))

###############################################################################
# Build single-row result and write to disk.
###############################################################################
result <- data.frame(
  task_id   = task_id,
  p         = cfg$p,
  n         = cfg$n,
  cov_model = cfg$cov_model,
  signal    = cfg$signal,
  rho       = cfg$rho,
  n_iter    = n_iter,
  
  power_BS_norm  = reject_counts["BS_norm"]  / n_iter,
  power_CQ_norm  = reject_counts["CQ_norm"]  / n_iter,
  power_SD_norm  = reject_counts["SD_norm"]  / n_iter,
  power_SKK_norm = reject_counts["SKK_norm"] / n_iter,
  power_BS_3c    = reject_counts["BS_3c"]    / n_iter,
  power_CQ_3c    = reject_counts["CQ_3c"]    / n_iter,
  power_ZZZ23    = reject_counts["ZZZ23"]    / n_iter,
  power_CLZ      = reject_counts["CLZ"]      / n_iter,
  
  na_BS_norm  = na_counts["BS_norm"],
  na_CQ_norm  = na_counts["CQ_norm"],
  na_SD_norm  = na_counts["SD_norm"],
  na_SKK_norm = na_counts["SKK_norm"],
  na_BS_3c    = na_counts["BS_3c"],
  na_CQ_3c    = na_counts["CQ_3c"],
  na_ZZZ23    = na_counts["ZZZ23"],
  na_CLZ      = na_counts["CLZ"],
  
  time_min  = round(elapsed, 2),
  row.names = NULL
)

outfile <- sprintf("results/result_task_%04d.csv", task_id)
write.csv(result, outfile, row.names = FALSE)

cat(sprintf("\nDone. Elapsed: %.1f min\n", elapsed))
cat(sprintf("Results written to: %s\n", outfile))
print(result)
