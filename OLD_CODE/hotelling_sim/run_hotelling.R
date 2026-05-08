make_mu_null <- function(p) {
  rep(0, p)
}

make_mu_sparse <- function(p, frac = 0.05, delta = 0.5) {
  mu <- rep(0, p)
  mu[1:floor(frac * p)] <- delta
  mu
}

make_mu_moderate <- function(p) {
  mu <- rep(0, p)
  mu[1:floor(0.20 * p)] <- 0.25
  mu
}

make_mu_dense <- function(p) {
  rep(sqrt(0.0125), p)
} 

generate_identity <- function(n, p, mu = rep(0, p)) {
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  if (any(mu != 0)) X <- sweep(X, 2, mu, "+")
  X
}

generate_cs <- function(n, p, rho = 0.5, mu = rep(0, p)) {
  Z <- matrix(rnorm(n * p), nrow = n, ncol = p)
  W <- rnorm(n)
  X <- sqrt(1 - rho) * Z + sqrt(rho) * W
  if (any(mu != 0)) X <- sweep(X, 2, mu, "+")
  X
}

generate_ar1 <- function(n, p, rho = 0.5, mu = rep(0, p)) {
  Z <- matrix(rnorm(n * p), nrow = n, ncol = p)
  scale <- sqrt(1 - rho^2)
  for (k in 2:p) Z[, k] <- rho * Z[, k - 1] + scale * Z[, k]
  if (any(mu != 0)) Z <- sweep(Z, 2, mu, "+")
  Z
}

generate_two_sample <- function(n, p, cov_model = "identity",
                                signal = "null", rho = 0.5) {
  mu2 <- switch(signal,
                "null"     = make_mu_null(p),
                "sparse"   = make_mu_sparse(p),
                "moderate" = make_mu_moderate(p),
                "dense"    = make_mu_dense(p)
  )
  gen <- switch(cov_model,
                "identity" = generate_identity,
                "cs"       = function(n, p, mu) generate_cs(n, p, rho = rho, mu = mu),
                "ar1"      = function(n, p, mu) generate_ar1(n, p, rho = rho, mu = mu)
  )
  list(X = gen(n, p, mu = rep(0, p)), Y = gen(n, p, mu = mu2))
}

hotelling_t2_stat <- function(X, Y) {
  n1 <- nrow(X)
  n2 <- nrow(Y)
  n_total <- n1 + n2
  denom <- n_total - 2
  
  mean_diff <- colMeans(X) - colMeans(Y)
  X_c <- sweep(X, 2, colMeans(X))
  Y_c <- sweep(Y, 2, colMeans(Y))
  resid <- rbind(X_c, Y_c)
  sv <- svd(resid, nu = 0)
  d <- sv$d
  V <- sv$v 
  tol <- max(d) * max(n_total, ncol(X)) * .Machine$double.eps
  keep <- d > tol
  r <- sum(keep)
  if (r == 0) return(0)
  d <- d[keep]
  V <- V[, keep, drop = FALSE]  
  proj <- as.numeric(crossprod(V, mean_diff))  
  T2 <- (n1 * n2) / (n1 + n2) * sum(proj^2 * denom / d^2)
  
  return(T2)
}

hotelling_perm_pval <- function(X, Y, n_perm = 1000) {
  n1 <- nrow(X)
  n2 <- nrow(Y)
  pooled <- rbind(X, Y)
  n_total <- n1 + n2
  T2_obs <- hotelling_t2_stat(X, Y)
  count <- 0
  for (b in 1:n_perm) {
    idx <- sample(n_total)
    X_perm <- pooled[idx[1:n1], , drop = FALSE]
    Y_perm <- pooled[idx[(n1 + 1):n_total], , drop = FALSE]
    T2_perm <- hotelling_t2_stat(X_perm, Y_perm)
    if (T2_perm >= T2_obs) count <- count + 1
  }
  
  return((count + 1) / (n_perm + 1))  # continuity-corrected p-value
}

# ---- Simulation grid (Tier 1 only: p <= 2000) ----
p_grid   <- c(200, 500, 2000)
n_grid   <- c(5, 10, 20, 50)
cov_grid <- c("identity", "cs", "ar1")
sig_grid <- c("null", "sparse", "moderate", "dense")

n_iter  <- 1000   # Monte Carlo iterations
n_perm  <- 1000   # permutations per iteration
alpha   <- 0.05
rho     <- 0.5

grid <- expand.grid(
  p = p_grid,
  n = n_grid,
  cov_model = cov_grid,
  signal = sig_grid,
  stringsAsFactors = FALSE
)

# ---- Read SLURM task ID ----
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (is.na(task_id) || task_id < 1 || task_id > nrow(grid)) {
  stop(sprintf("Invalid SLURM_ARRAY_TASK_ID: %s (grid has %d rows)",
               Sys.getenv("SLURM_ARRAY_TASK_ID"), nrow(grid)))
}

cfg <- grid[task_id, ]
cat(sprintf("Task %d: p=%d, n=%d, cov=%s, signal=%s, n_iter=%d, n_perm=%d\n",
            task_id, cfg$p, cfg$n, cfg$cov_model, cfg$signal, n_iter, n_perm))

set.seed(2026 + task_id)

# ---- Main simulation loop ----
reject_ht <- 0

t_start <- Sys.time()

for (k in 1:n_iter) {
  dat <- generate_two_sample(cfg$n, cfg$p,
                             cov_model = cfg$cov_model,
                             signal = cfg$signal, rho = rho)
  
  pv_ht <- tryCatch(
    hotelling_perm_pval(dat$X, dat$Y, n_perm = n_perm),
    error = function(e) NA
  )
  
  if (!is.na(pv_ht) && pv_ht < alpha) reject_ht <- reject_ht + 1
  if (k %% 50 == 0) cat(sprintf("  Iteration %d/%d\n", k, n_iter))
}

t_end <- Sys.time()
elapsed <- as.numeric(difftime(t_end, t_start, units = "mins"))

result <- data.frame(
  task_id   = task_id,
  p = cfg$p,
  n = cfg$n,
  cov_model = cfg$cov_model,
  signal = cfg$signal,
  n_iter = n_iter,
  n_perm = n_perm,
  power_HT  = reject_ht / n_iter,
  time_min  = round(elapsed, 2)
)

outfile <- sprintf("results_hotelling/result_hotelling_%03d.csv", task_id)
write.csv(result, outfile, row.names = FALSE)

cat(sprintf("\nDone! Elapsed: %.1f min\n", elapsed))
cat(sprintf("Results saved to %s\n", outfile))
print(result)
