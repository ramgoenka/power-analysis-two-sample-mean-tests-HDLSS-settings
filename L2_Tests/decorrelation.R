###############################################################################
# decorrelation.R  -- (Compound Symmetry only)
#
# decorr modes:
#   none
#   oracle  -- apply Sigma^{-1/2} with TRUE rho
#   plugin  -- apply Sigma^{-1/2} with rho_hat from pooled sample
###############################################################################
lib_path <- Sys.getenv("R_LIBS_USER")
if (lib_path == "") lib_path <- file.path(Sys.getenv("HOME"), "R_libs")
.libPaths(c(lib_path, .libPaths()))

suppressPackageStartupMessages({
  library(HDNRA)
  library(highmean)
})

###############################################################################
# Decorrelation helpers
###############################################################################
# Apply CS decorrelation with a given rho.
# Sigma^{-1/2} = a*(I - (1/p)*1*1') + b*(1/p)*1*1'
# a = 1/sqrt(1 - rho), b = 1/sqrt(1 + (p-1)*rho)
# On a row x: x' = a*(x - mean(x)) + b*mean(x). O(np) per matrix.
decorrelate_cs <- function(X, rho) {
  p <- ncol(X)
  rho <- max(min(rho, 0.999), -1/(p - 1) + 1e-6)
  a <- 1 / sqrt(1 - rho)
  b <- 1 / sqrt(1 + (p - 1) * rho)
  m_vec <- rowMeans(X)
  X_c   <- sweep(X, 1, m_vec, "-")
  sweep(a * X_c, 1, b * m_vec, "+")
}

# Plug-in rho_hat = mean off-diagonal of pooled sample correlation matrix.
# Computed via row-sum trick: sum(R_hat) = ||Z * 1_p||^2 / df. No p x p object.
estimate_cs_rho <- function(X, Y) {
  n_x <- nrow(X); n_y <- nrow(Y); p <- ncol(X)
  N   <- n_x + n_y
  df  <- N - 2
  X_c <- sweep(X, 2, colMeans(X), "-")
  Y_c <- sweep(Y, 2, colMeans(Y), "-")
  s2  <- (colSums(X_c^2) + colSums(Y_c^2)) / df
  s   <- sqrt(pmax(s2, 1e-12))
  Z   <- rbind(sweep(X_c, 2, s, "/"), sweep(Y_c, 2, s, "/"))
  row_sums <- rowSums(Z)
  sum_R    <- sum(row_sums^2) / df
  rho_hat  <- (sum_R - p) / (p * (p - 1))
  max(min(rho_hat, 0.999), -1/(p - 1) + 1e-6)
}
decorrelate_pair <- function(X, Y, mode, rho_true = NA_real_) {
  if (mode == "none") return(list(X = X, Y = Y, rho_used = NA_real_))
  rho_use <- if (mode == "oracle") rho_true else estimate_cs_rho(X, Y)
  list(X = decorrelate_cs(X, rho_use),
       Y = decorrelate_cs(Y, rho_use),
       rho_used = rho_use)
}

###############################################################################
# Mean-shift constructors  (equal-energy: ||mu2||^2 = 0.0125 * p)
###############################################################################
make_mu_null     <- function(p) rep(0, p)
make_mu_sparse   <- function(p) { mu <- rep(0, p); mu[1:floor(0.05 * p)] <- 0.5;  mu }
make_mu_moderate <- function(p) { mu <- rep(0, p); mu[1:floor(0.20 * p)] <- 0.25; mu }
make_mu_dense    <- function(p) rep(sqrt(0.0125), p)

###############################################################################
# CS data generator
###############################################################################
generate_cs <- function(n, p, rho, mu = rep(0, p)) {
  Z <- matrix(rnorm(n * p), nrow = n, ncol = p)
  W <- rnorm(n)
  X <- sqrt(1 - rho) * Z + sqrt(rho) * matrix(W, nrow = n, ncol = p)
  if (any(mu != 0)) X <- sweep(X, 2, mu, "+")
  X
}

generate_two_sample <- function(n, p, signal, rho) {
  mu2 <- switch(signal,
                "null"     = make_mu_null(p),
                "sparse"   = make_mu_sparse(p),
                "moderate" = make_mu_moderate(p),
                "dense"    = make_mu_dense(p),
                stop("Unknown signal: ", signal))
  list(X = generate_cs(n, p, rho),
       Y = generate_cs(n, p, rho, mu = mu2))
}

###############################################################################
# Test runner
###############################################################################
run_all_tests <- function(X, Y) {
  safe <- function(expr) tryCatch(expr, error = function(e) NA_real_)
  c(
    BS_norm = safe(BS1996.TS.NABT(X, Y)$p.value),
    CQ_norm = safe(CQ2010.TSBF.NABT(X, Y)$p.value),
    SD_norm = safe(SD2008.TS.NABT(X, Y)$p.value),
    BS_3c   = safe(ZZ2022.TS.3cNRT(X, Y)$p.value),
    CQ_3c   = safe(ZZ2022.TSBF.3cNRT(X, Y)$p.value)
  )
}

###############################################################################
# Pilot grid
###############################################################################
p_grid      <- c(10000)
n_grid      <- c(5, 10, 20, 50, 100)
sig_grid    <- c("null", "sparse", "moderate", "dense")
rho_grid    <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.8)
decorr_grid <- c("none", "oracle", "plugin")

n_iter <- 5000
alpha  <- 0.05

grid <- expand.grid(
  p = p_grid, n = n_grid, cov_model = "cs",
  signal = sig_grid, rho = rho_grid, decorr = decorr_grid,
  stringsAsFactors = FALSE
)

###############################################################################
# Task dispatch
###############################################################################
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (is.na(task_id) || task_id < 1 || task_id > nrow(grid)) {
  stop(sprintf("Invalid SLURM_ARRAY_TASK_ID: %s (grid has %d rows)",
               Sys.getenv("SLURM_ARRAY_TASK_ID"), nrow(grid)))
}

cfg <- grid[task_id, ]
cat(sprintf("Task %d / %d:\n", task_id, nrow(grid)))
cat(sprintf("  p=%d, n=%d, cov=%s, signal=%s, rho=%.2f, decorr=%s, B=%d\n",
            cfg$p, cfg$n, cfg$cov_model, cfg$signal, cfg$rho, cfg$decorr, n_iter))

set.seed(2026 + task_id)

test_names    <- c("BS_norm", "CQ_norm", "SD_norm",
                   "BS_3c",   "CQ_3c")
reject_counts <- setNames(integer(length(test_names)), test_names)
na_counts     <- setNames(integer(length(test_names)), test_names)

rho_used_sum   <- 0
rho_used_sumsq <- 0
rho_used_count <- 0

t_start <- Sys.time()

for (k in 1:n_iter) {
  dat <- generate_two_sample(cfg$n, cfg$p, signal = cfg$signal, rho = cfg$rho)
  dec <- decorrelate_pair(dat$X, dat$Y, mode = cfg$decorr, rho_true = cfg$rho)
  
  if (!is.na(dec$rho_used)) {
    rho_used_sum   <- rho_used_sum + dec$rho_used
    rho_used_sumsq <- rho_used_sumsq + dec$rho_used^2
    rho_used_count <- rho_used_count + 1
  }
  
  pvals <- run_all_tests(dec$X, dec$Y)
  rejects <- !is.na(pvals) & pvals < alpha
  reject_counts <- reject_counts + rejects
  na_counts     <- na_counts + is.na(pvals)
  
  if (k %% 500 == 0) {
    elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
    cat(sprintf("  iter %d/%d  (%.1f min elapsed)\n", k, n_iter, elapsed))
  }
}

elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))

rho_used_mean <- if (rho_used_count > 0) rho_used_sum / rho_used_count else NA_real_
rho_used_sd   <- if (rho_used_count > 1) {
  sqrt((rho_used_sumsq - rho_used_count * rho_used_mean^2) / (rho_used_count - 1))
} else NA_real_

###############################################################################
# Result row
###############################################################################
result <- data.frame(
  task_id        = task_id,
  p              = cfg$p,
  n              = cfg$n,
  cov_model      = cfg$cov_model,
  signal         = cfg$signal,
  rho            = cfg$rho,
  decorr         = cfg$decorr,
  n_iter         = n_iter,
  
  rho_used_mean  = rho_used_mean,
  rho_used_sd    = rho_used_sd,
  
  power_BS_norm  = reject_counts["BS_norm"]  / n_iter,
  power_CQ_norm  = reject_counts["CQ_norm"]  / n_iter,
  power_SD_norm  = reject_counts["SD_norm"]  / n_iter,
  power_BS_3c    = reject_counts["BS_3c"]    / n_iter,
  power_CQ_3c    = reject_counts["CQ_3c"]    / n_iter,
  na_BS_norm     = na_counts["BS_norm"],
  na_CQ_norm     = na_counts["CQ_norm"],
  na_SD_norm     = na_counts["SD_norm"],
  na_BS_3c       = na_counts["BS_3c"],
  na_CQ_3c       = na_counts["CQ_3c"],
  time_min       = round(elapsed, 2),
  row.names      = NULL
)

outfile <- sprintf("results/result_task_%04d.csv", task_id)
write.csv(result, outfile, row.names = FALSE)

cat(sprintf("\nDone. Elapsed: %.1f min\n", elapsed))
cat(sprintf("rho_used: mean=%.4f, sd=%.4f (over %d iters)\n",
            rho_used_mean, ifelse(is.na(rho_used_sd), 0, rho_used_sd), rho_used_count))
cat(sprintf("Results written to: %s\n", outfile))
print(result)
