###############################################################################
# install_packages.R
###############################################################################

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

for (pkg in required_pkgs) {
  library(pkg, character.only = TRUE, lib.loc = lib_path)
  cat(sprintf("  [OK] %s loaded successfully\n", pkg))
}

cat(sprintf("\nAll packages ready. Library path: %s\n", lib_path))
