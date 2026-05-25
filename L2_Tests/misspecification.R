###############################################################################
# misspecification.R -- CS-decorrelation under misspecified covariances
#
# Purpose:
#   Stress-test the CS-decorrelation preprocessing step when the true covariance
#   is not exactly compound symmetric.
#
# Covariance models:
#   identity   -- independent coordinates
#   cs         -- compound symmetry, Corr(j,k)=rho for j != k
#   ar1        -- AR(1), Corr(j,k)=rho^|j-k|
#   two_block  -- 51/49 split, two independent CS blocks,
#                 within-block Corr=rho, between=0
#
# Preprocessing modes evaluated on the SAME simulated dataset per iteration
# (common random numbers for variance reduction across decorr comparisons):
#   none       -- no preprocessing
#   cs_oracle  -- apply CS Sigma^{-1/2} using the population average off-diagonal
#                 correlation rho_pop_bar. For true CS this equals rho. For
#                 misspecified models it is the best CS-summary parameter,
#                 not the true inverse. Naming makes this explicit: this is
#                 the CS-APPROXIMATION oracle, not the true-covariance oracle.
#   cs_plugin  -- apply CS Sigma^{-1/2} using rho_hat estimated from pooled data
###############################################################################
lib_path <- Sys.getenv("R_LIBS_USER")
if (lib_path == "") lib_path <- file.path(Sys.getenv("HOME"), "R_libs")
.libPaths(c(lib_path, .libPaths()))

suppressPackageStartupMessages({
  library(HDNRA)
  library(highmean)
})

###############################################################################
# Small helpers
###############################################################################
parse_num_vec <- function(x, default) {
  if (is.na(x) || x == "") return(default)
  as.numeric(strsplit(x, ",")[[1]])
}

parse_int_vec <- function(x, default) {
  as.integer(parse_num_vec(x, default))
}

# Clamp rho to the strictly-valid CS range, with a tight upper bound to avoid
# pathological a^2 = 1/(1-rho) blow-up. At rho=0.97 the off-spike amplifier
# is a^2 = 33 and the spike compressor at p=2000 is b^2 ~ 0.0003, which is
# aggressive but bounded. At rho=0.999, a^2 = 1000, which dominates any signal.
clamp_cs_rho <- function(rho, p, log_clip = FALSE) {
  upper <- 0.97
  lower <- -1 / (p - 1) + 1e-6
  clipped <- (rho > upper) || (rho < lower)
  rho_out <- max(min(rho, upper), lower)
  if (log_clip) attr(rho_out, "clipped") <- isTRUE(clipped)
  rho_out
}

###############################################################################
# CS decorrelation helpers
###############################################################################
# Apply CS decorrelation with a given rho.
# Sigma^{-1/2} = a*(I - P_1) + b*P_1, where P_1 = 11'/p.
# On a row x: x' = a*(x - mean(x)1) + b*mean(x)1.
decorrelate_cs <- function(X, rho) {
  p <- ncol(X)
  rho <- clamp_cs_rho(rho, p)
  a <- 1 / sqrt(1 - rho)
  b <- 1 / sqrt(1 + (p - 1) * rho)
  m_vec <- rowMeans(X)
  X_c   <- sweep(X, 1, m_vec, "-")
  sweep(a * X_c, 1, b * m_vec, "+")
}

# Plug-in rho_hat = mean off-diagonal of pooled sample correlation matrix.
# Uses a row-sum trick, avoiding construction of a p x p correlation matrix.
# Returns rho_hat with attributes "clipped" (logical) and "raw" (pre-clip value)
# for downstream diagnostics.
estimate_cs_rho <- function(X, Y) {
  n_x <- nrow(X); n_y <- nrow(Y); p <- ncol(X)
  N   <- n_x + n_y
  df  <- N - 2
  
  X_c <- sweep(X, 2, colMeans(X), "-")
  Y_c <- sweep(Y, 2, colMeans(Y), "-")
  s2  <- (colSums(X_c^2) + colSums(Y_c^2)) / df
  s   <- sqrt(pmax(s2, 1e-12))
  
  Z <- rbind(sweep(X_c, 2, s, "/"), sweep(Y_c, 2, s, "/"))
  row_sums <- rowSums(Z)
  sum_R    <- sum(row_sums^2) / df
  rho_hat_raw <- (sum_R - p) / (p * (p - 1))
  
  rho_hat <- clamp_cs_rho(rho_hat_raw, p, log_clip = TRUE)
  attr(rho_hat, "raw") <- rho_hat_raw
  rho_hat
}

###############################################################################
# Population CS-summary rho under each covariance model
###############################################################################
# Population average off-diagonal correlation:
#   rho_bar = [sum_{j != k} Corr(X_j, X_k)] / [p(p-1)].
# For true CS, rho_bar = rho. For AR(1) and block covariance it is only the
# scalar CS approximation, useful for an oracle misspecification benchmark.
pop_avg_offdiag_corr <- function(cov_model, p, rho) {
  if (cov_model == "identity") return(0)
  if (cov_model == "cs") return(rho)
  
  if (cov_model == "ar1") {
    # Ordered off-diagonal sum: 2 * sum_{lag=1}^{p-1} (p-lag) rho^lag.
    lags <- seq_len(p - 1)
    return(2 * sum((p - lags) * rho^lags) / (p * (p - 1)))
  }
  
  if (cov_model == "two_block") {
    # 51/49 split (mildly imbalanced to break the 50/50 symmetry).
    b1 <- ceiling(p * 0.51)
    b2 <- p - b1
    # Ordered within-block off-diagonal entries only; between-block corr = 0.
    return(rho * (b1 * (b1 - 1) + b2 * (b2 - 1)) / (p * (p - 1)))
  }
  
  stop("Unknown cov_model: ", cov_model)
}

# Quadratic form under the CS inverse with parameter rho_cs.
# Useful for interpreting whether CS whitening should preserve or remove signal.
cs_inv_quad <- function(mu, rho_cs) {
  p <- length(mu)
  rho_cs <- clamp_cs_rho(rho_cs, p)
  norm2 <- sum(mu^2)
  s1    <- sum(mu)
  norm2 / (1 - rho_cs) -
    (rho_cs * s1^2) / ((1 - rho_cs) * (1 - rho_cs + p * rho_cs))
}

# Squared cosine of mu with 1_p/sqrt(p): the "alpha^2" alignment from the
# alignment-experiment theory. Useful when analyzing whether the alignment
# story from CS carries over to misspecified covariances.
alignment_with_ones <- function(mu) {
  p <- length(mu)
  norm2 <- sum(mu^2)
  if (norm2 <= 0) return(NA_real_)
  sum(mu)^2 / (p * norm2)
}

###############################################################################
# Mean-shift constructors: ||mu2||^2 = 0.0125 * p
###############################################################################
make_mu_null <- function(p) rep(0, p)
make_mu_sparse <- function(p) {
  mu <- rep(0, p)
  mu[1:floor(0.05 * p)] <- 0.5
  mu
}
make_mu_moderate <- function(p) {
  mu <- rep(0, p)
  mu[1:floor(0.20 * p)] <- 0.25
  mu
}
make_mu_dense <- function(p) rep(sqrt(0.0125), p)

make_mu <- function(p, signal) {
  switch(signal,
         "null"     = make_mu_null(p),
         "sparse"   = make_mu_sparse(p),
         "moderate" = make_mu_moderate(p),
         "dense"    = make_mu_dense(p),
         stop("Unknown signal: ", signal))
}

###############################################################################
# Data generators: O(np), no dense p x p covariance matrix formed
###############################################################################
generate_identity <- function(n, p, mu = rep(0, p)) {
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  if (any(mu != 0)) X <- sweep(X, 2, mu, "+")
  X
}

generate_cs <- function(n, p, rho, mu = rep(0, p)) {
  Z <- matrix(rnorm(n * p), nrow = n, ncol = p)
  W <- rnorm(n)
  X <- sqrt(1 - rho) * Z + sqrt(rho) * matrix(W, nrow = n, ncol = p)
  if (any(mu != 0)) X <- sweep(X, 2, mu, "+")
  X
}

generate_ar1 <- function(n, p, rho, mu = rep(0, p)) {
  X <- matrix(0, nrow = n, ncol = p)
  X[, 1] <- rnorm(n)
  innov_sd <- sqrt(1 - rho^2)
  if (p >= 2) {
    for (j in 2:p) {
      X[, j] <- rho * X[, j - 1] + innov_sd * rnorm(n)
    }
  }
  if (any(mu != 0)) X <- sweep(X, 2, mu, "+")
  X
}

# Two independent CS blocks with a 51/49 split (mildly imbalanced to break
# the 50/50 symmetry: at 50/50 the two block-mean correlations are
# statistically interchangeable, hiding heterogeneity. 51/49 perturbs that
# enough to expose it without changing rho_pop_bar materially.).
generate_two_block <- function(n, p, rho, mu = rep(0, p)) {
  b1 <- ceiling(p * 0.51)
  b2 <- p - b1
  
  Z1 <- matrix(rnorm(n * b1), nrow = n, ncol = b1)
  Z2 <- matrix(rnorm(n * b2), nrow = n, ncol = b2)
  W1 <- rnorm(n)
  W2 <- rnorm(n)
  
  X1 <- sqrt(1 - rho) * Z1 + sqrt(rho) * matrix(W1, nrow = n, ncol = b1)
  X2 <- sqrt(1 - rho) * Z2 + sqrt(rho) * matrix(W2, nrow = n, ncol = b2)
  X  <- cbind(X1, X2)
  
  if (any(mu != 0)) X <- sweep(X, 2, mu, "+")
  X
}

generate_sample <- function(n, p, cov_model, rho, mu = rep(0, p)) {
  switch(cov_model,
         "identity"  = generate_identity(n, p, mu),
         "cs"        = generate_cs(n, p, rho, mu),
         "ar1"       = generate_ar1(n, p, rho, mu),
         "two_block" = generate_two_block(n, p, rho, mu),
         stop("Unknown cov_model: ", cov_model))
}

generate_two_sample <- function(n, p, cov_model, rho, signal) {
  mu2 <- make_mu(p, signal)
  list(
    X   = generate_sample(n, p, cov_model, rho, mu = rep(0, p)),
    Y   = generate_sample(n, p, cov_model, rho, mu = mu2),
    mu2 = mu2
  )
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
# Simulation grid
###############################################################################
# Reduced from the original draft based on focused experimental design:
#   - n in {20, 50}: where rho_hat has interesting finite-sample variance
#       (n=5 has size issues with the underlying tests; n=100 is too stable)
#   - p in {500, 2000, 10000}: misspecification effects scale with p
#   - rho in {0.1, 0.3, 0.5, 0.7, 0.8}: endpoints + three interior points,
#       sufficient for either a sweep plot or fixed-rho summary tables.
p_grid   <- parse_int_vec(Sys.getenv("P_GRID"),   c(500, 2000, 10000))
n_grid   <- parse_int_vec(Sys.getenv("N_GRID"),   c(20, 50))
rho_grid <- parse_num_vec(Sys.getenv("RHO_GRID"), c(0.1, 0.3, 0.5, 0.7, 0.8))

cov_grid <- c("identity", "cs", "ar1", "two_block")
sig_grid <- c("null", "sparse", "moderate", "dense")

decorr_modes <- c("none", "cs_oracle", "cs_plugin")
test_names   <- c("BS_norm", "CQ_norm", "SD_norm", "BS_3c", "CQ_3c")

n_iter <- as.integer(Sys.getenv("N_ITER", unset = "5000"))
alpha  <- as.numeric(Sys.getenv("ALPHA", unset = "0.05"))

grid_nonid <- expand.grid(
  p = p_grid, n = n_grid,
  cov_model = c("cs", "ar1", "two_block"),
  signal = sig_grid, rho = rho_grid,
  stringsAsFactors = FALSE
)

grid_id <- expand.grid(
  p = p_grid, n = n_grid,
  cov_model = "identity",
  signal = sig_grid, rho = 0,
  stringsAsFactors = FALSE
)

grid <- rbind(grid_id, grid_nonid)
grid <- grid[order(grid$p, grid$n, grid$cov_model, grid$signal, grid$rho), ]
row.names(grid) <- NULL

###############################################################################
# Task dispatch
###############################################################################
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (is.na(task_id) || task_id < 1 || task_id > nrow(grid)) {
  stop(sprintf("Invalid SLURM_ARRAY_TASK_ID: %s (grid has %d rows)",
               Sys.getenv("SLURM_ARRAY_TASK_ID"), nrow(grid)))
}

cfg <- grid[task_id, ]

# Seed design:
#   - task_id is the SLURM array index
#   - SCRIPT_ID separates this script's random streams from sibling scripts
#     that may share task_id (e.g., the alignment experiment)
#   - All decorrelation modes are evaluated inside this task on the same
#     simulated data per iteration (common random numbers)
base_seed <- as.integer(Sys.getenv("BASE_SEED", unset = "2026"))
script_id <- as.integer(Sys.getenv("SCRIPT_ID", unset = "31"))
seed <- base_seed + 100000L * script_id + task_id

set.seed(seed)

cat(sprintf("Task %d / %d:\n", task_id, nrow(grid)))
cat(sprintf("  p=%d, n=%d, cov=%s, signal=%s, rho=%.3f, B=%d, alpha=%.3f\n",
            cfg$p, cfg$n, cfg$cov_model, cfg$signal, cfg$rho, n_iter, alpha))
cat(sprintf("  seed=%d\n", seed))

rho_pop_bar <- pop_avg_offdiag_corr(cfg$cov_model, cfg$p, cfg$rho)

mu2_cfg <- make_mu(cfg$p, cfg$signal)
mu2_norm2 <- sum(mu2_cfg^2)
mu2_align_ones <- alignment_with_ones(mu2_cfg)
mu2_cs_inv_quad_oracle <- if (mu2_norm2 > 0) cs_inv_quad(mu2_cfg, rho_pop_bar) else NA_real_

cat(sprintf("  population avg offdiag rho_bar=%.6f\n", rho_pop_bar))
cat(sprintf("  ||mu2||^2=%.4f, align_ones=%s, CS-inv quad=%.4f\n",
            mu2_norm2,
            ifelse(is.na(mu2_align_ones), "NA", sprintf("%.4f", mu2_align_ones)),
            ifelse(is.na(mu2_cs_inv_quad_oracle), NA_real_, mu2_cs_inv_quad_oracle)))

reject_counts <- array(0L,
                       dim = c(length(decorr_modes), length(test_names)),
                       dimnames = list(decorr_modes, test_names))
na_counts <- array(0L,
                   dim = c(length(decorr_modes), length(test_names)),
                   dimnames = list(decorr_modes, test_names))

rho_used_sum   <- setNames(rep(0, length(decorr_modes)), decorr_modes)
rho_used_sumsq <- setNames(rep(0, length(decorr_modes)), decorr_modes)
rho_used_count <- setNames(rep(0L, length(decorr_modes)), decorr_modes)

rhohat_plugin_sum   <- 0
rhohat_plugin_sumsq <- 0
rhohat_plugin_raw_sum   <- 0  # pre-clip, to detect when clipping bites
rhohat_plugin_raw_sumsq <- 0
rhohat_clip_count   <- 0L     # number of iterations where rho_hat got clipped

###############################################################################
# Monte Carlo loop
###############################################################################
t_start <- Sys.time()

for (k in 1:n_iter) {
  dat <- generate_two_sample(cfg$n, cfg$p, cfg$cov_model, cfg$rho, cfg$signal)
  
  rho_hat <- estimate_cs_rho(dat$X, dat$Y)
  rho_hat_raw <- attr(rho_hat, "raw")
  if (isTRUE(attr(rho_hat, "clipped"))) {
    rhohat_clip_count <- rhohat_clip_count + 1L
  }
  rhohat_plugin_sum     <- rhohat_plugin_sum + as.numeric(rho_hat)
  rhohat_plugin_sumsq   <- rhohat_plugin_sumsq + as.numeric(rho_hat)^2
  rhohat_plugin_raw_sum   <- rhohat_plugin_raw_sum + rho_hat_raw
  rhohat_plugin_raw_sumsq <- rhohat_plugin_raw_sumsq + rho_hat_raw^2
  
  for (mode in decorr_modes) {
    if (mode == "none") {
      X_use <- dat$X
      Y_use <- dat$Y
      rho_use <- NA_real_
    } else if (mode == "cs_oracle") {
      rho_use <- rho_pop_bar
      X_use <- decorrelate_cs(dat$X, rho_use)
      Y_use <- decorrelate_cs(dat$Y, rho_use)
    } else if (mode == "cs_plugin") {
      rho_use <- as.numeric(rho_hat)
      X_use <- decorrelate_cs(dat$X, rho_use)
      Y_use <- decorrelate_cs(dat$Y, rho_use)
    } else {
      stop("Unknown decorr mode: ", mode)
    }
    
    if (!is.na(rho_use)) {
      rho_used_sum[mode]   <- rho_used_sum[mode] + rho_use
      rho_used_sumsq[mode] <- rho_used_sumsq[mode] + rho_use^2
      rho_used_count[mode] <- rho_used_count[mode] + 1L
    }
    
    pvals <- run_all_tests(X_use, Y_use)
    rejects <- !is.na(pvals) & pvals < alpha
    reject_counts[mode, ] <- reject_counts[mode, ] + rejects
    na_counts[mode, ]     <- na_counts[mode, ] + is.na(pvals)
  }
  
  if (k %% 500 == 0) {
    elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
    cat(sprintf("  iter %d/%d  (%.1f min elapsed)\n", k, n_iter, elapsed))
  }
}

elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))

rhohat_plugin_mean <- rhohat_plugin_sum / n_iter
rhohat_plugin_sd <- if (n_iter > 1) {
  sqrt((rhohat_plugin_sumsq - n_iter * rhohat_plugin_mean^2) / (n_iter - 1))
} else NA_real_

rhohat_plugin_raw_mean <- rhohat_plugin_raw_sum / n_iter
rhohat_plugin_raw_sd <- if (n_iter > 1) {
  sqrt((rhohat_plugin_raw_sumsq - n_iter * rhohat_plugin_raw_mean^2) / (n_iter - 1))
} else NA_real_

rho_used_mean <- rep(NA_real_, length(decorr_modes))
rho_used_sd   <- rep(NA_real_, length(decorr_modes))
names(rho_used_mean) <- names(rho_used_sd) <- decorr_modes

for (mode in decorr_modes) {
  if (rho_used_count[mode] > 0) {
    rho_used_mean[mode] <- rho_used_sum[mode] / rho_used_count[mode]
    rho_used_sd[mode] <- if (rho_used_count[mode] > 1) {
      sqrt((rho_used_sumsq[mode] - rho_used_count[mode] * rho_used_mean[mode]^2) /
             (rho_used_count[mode] - 1))
    } else 0
  }
}

###############################################################################
# Result rows: one row per preprocessing mode
###############################################################################
result <- data.frame(
  task_id        = task_id,
  seed           = seed,
  p              = cfg$p,
  n              = cfg$n,
  cov_model      = cfg$cov_model,
  signal         = cfg$signal,
  rho            = cfg$rho,
  rho_pop_bar    = rho_pop_bar,
  decorr         = decorr_modes,
  n_iter         = n_iter,
  
  mu2_norm2          = mu2_norm2,
  align_ones         = mu2_align_ones,
  cs_inv_quad_oracle = mu2_cs_inv_quad_oracle,
  
  rhohat_plugin_mean     = rhohat_plugin_mean,
  rhohat_plugin_sd       = rhohat_plugin_sd,
  rhohat_plugin_raw_mean = rhohat_plugin_raw_mean,
  rhohat_plugin_raw_sd   = rhohat_plugin_raw_sd,
  rhohat_clip_count      = rhohat_clip_count,
  rho_used_mean          = as.numeric(rho_used_mean[decorr_modes]),
  rho_used_sd            = as.numeric(rho_used_sd[decorr_modes]),
  
  power_BS_norm  = as.numeric(reject_counts[, "BS_norm"] / n_iter),
  power_CQ_norm  = as.numeric(reject_counts[, "CQ_norm"] / n_iter),
  power_SD_norm  = as.numeric(reject_counts[, "SD_norm"] / n_iter),
  power_BS_3c    = as.numeric(reject_counts[, "BS_3c"]   / n_iter),
  power_CQ_3c    = as.numeric(reject_counts[, "CQ_3c"]   / n_iter),
  
  na_BS_norm     = as.integer(na_counts[, "BS_norm"]),
  na_CQ_norm     = as.integer(na_counts[, "CQ_norm"]),
  na_SD_norm     = as.integer(na_counts[, "SD_norm"]),
  na_BS_3c       = as.integer(na_counts[, "BS_3c"]),
  na_CQ_3c       = as.integer(na_counts[, "CQ_3c"]),
  
  time_min       = round(elapsed, 2),
  row.names      = NULL
)

dir.create("results", showWarnings = FALSE, recursive = TRUE)
outfile <- sprintf("results/misspec_task_%04d.csv", task_id)
write.csv(result, outfile, row.names = FALSE)

cat(sprintf("\nDone. Elapsed: %.1f min\n", elapsed))
cat(sprintf("plugin rho_hat (clipped): mean=%.6f, sd=%.6f\n",
            rhohat_plugin_mean, rhohat_plugin_sd))
cat(sprintf("plugin rho_hat (raw):     mean=%.6f, sd=%.6f\n",
            rhohat_plugin_raw_mean, rhohat_plugin_raw_sd))
cat(sprintf("rho_hat clipped in %d / %d iterations (%.2f%%)\n",
            rhohat_clip_count, n_iter, 100 * rhohat_clip_count / n_iter))
cat(sprintf("Results written to: %s\n", outfile))
print(result)