###############################################################################
# run_config.R
###############################################################################

library(highmean)
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

make_mu_dense  <- function(p) {
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

p_grid   <- c(200, 500, 2000, 10000)
n_grid   <- c(5, 10, 20, 50)
cov_grid <- c("identity", "cs", "ar1")
sig_grid <- c("null", "sparse", "moderate", "dense")

n_iter <- 1000
alpha  <- 0.05
rho    <- 0.5

grid <- expand.grid(
  p = p_grid,
  n = n_grid,
  cov_model = cov_grid,
  signal = sig_grid,
  stringsAsFactors = FALSE
)

task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (is.na(task_id) || task_id < 1 || task_id > nrow(grid)) {
  stop(sprintf("Invalid SLURM_ARRAY_TASK_ID: %s (grid has %d rows)", 
               Sys.getenv("SLURM_ARRAY_TASK_ID"), nrow(grid)))
}

cfg <- grid[task_id, ]
cat(sprintf("Task %d: p=%d, n=%d, cov=%s, signal=%s, n_iter=%d\n",
            task_id, cfg$p, cfg$n, cfg$cov_model, cfg$signal, n_iter))

set.seed(2026 + task_id)  

reject_bs <- 0
reject_cq <- 0
reject_sd <- 0

t_start <- Sys.time()

for (k in 1:n_iter) {
  dat <- generate_two_sample(cfg$n, cfg$p, 
                             cov_model = cfg$cov_model,
                             signal = cfg$signal, rho = rho)
  
  pv_bs <- tryCatch(apval_Bai1996(dat$X, dat$Y)$pval, error = function(e) NA)
  pv_cq <- tryCatch(apval_Chen2010(dat$X, dat$Y, eq.cov = FALSE)$pval, error = function(e) NA)
  pv_sd <- tryCatch(apval_Sri2008(dat$X, dat$Y)$pval, error = function(e) NA)
  
  if (!is.na(pv_bs) && pv_bs < alpha) reject_bs <- reject_bs + 1
  if (!is.na(pv_cq) && pv_cq < alpha) reject_cq <- reject_cq + 1
  if (!is.na(pv_sd) && pv_sd < alpha) reject_sd <- reject_sd + 1
  if (k %% 100 == 0) cat(sprintf("  Iteration %d/%d\n", k, n_iter))
}

t_end <- Sys.time()
elapsed <- as.numeric(difftime(t_end, t_start, units = "mins"))

result <- data.frame(
  task_id = task_id,
  p = cfg$p,
  n = cfg$n,
  cov_model = cfg$cov_model,
  signal = cfg$signal,
  n_iter = n_iter,
  power_BS = reject_bs / n_iter,
  power_CQ = reject_cq / n_iter,
  power_SD = reject_sd / n_iter,
  time_min = round(elapsed, 2)
)

outfile <- sprintf("results/result_task_%03d.csv", task_id)
write.csv(result, outfile, row.names = FALSE)

cat(sprintf("\nDone! Elapsed: %.1f min\n", elapsed))
cat(sprintf("Results saved to %s\n", outfile))
print(result)
