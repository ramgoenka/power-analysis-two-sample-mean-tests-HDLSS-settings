###############################################################################
# install_packages.R
#
# One-time setup script. Installs all R packages required for the L2-family
# simulation:
#   - highmean: provides apval_Chen2014 (the CLZ thresholding test)
#   - HDNRA:    provides the seven other L2 tests (BS, CQ, SD, SKK, the two
#               Zhang-Zhu 2022 normal-reference 3-cumulant tests, and ZZZ2023)
#
# Run ONCE on CRCD in an interactive session, NOT as a batch job:
#
#   srun --pty --cluster=smp --partition=smp --time=01:00:00 --mem=8G bash
#   module purge
#   module load gcc r
#   Rscript install_packages.R
###############################################################################

# Set up personal library
lib_path <- Sys.getenv("R_LIBS_USER")
if (lib_path == "") lib_path <- file.path(Sys.getenv("HOME"), "R_libs")
if (!dir.exists(lib_path)) dir.create(lib_path, recursive = TRUE)
.libPaths(c(lib_path, .libPaths()))

required_pkgs <- c("highmean", "HDNRA")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("Installing %s ...\n", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org", lib = lib_path)
  } else {
    cat(sprintf("%s already installed.\n", pkg))
  }
}

# Verify each loads cleanly
for (pkg in required_pkgs) {
  library(pkg, character.only = TRUE, lib.loc = lib_path)
  cat(sprintf("  [OK] %s loaded successfully\n", pkg))
}

cat(sprintf("\nAll packages ready. Library path: %s\n", lib_path))
