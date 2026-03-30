###############################################################################
# install_packages.R
###############################################################################

lib_path <- Sys.getenv("R_LIBS_USER")
if (lib_path == "") lib_path <- file.path(Sys.getenv("HOME"), "R_libs")
if (!dir.exists(lib_path)) dir.create(lib_path, recursive = TRUE)
.libPaths(c(lib_path, .libPaths()))

install.packages("highmean", repos = "https://cloud.r-project.org", lib = lib_path)

library(highmean, lib.loc = lib_path)
cat("highmean installed and loaded successfully!\n")
cat(sprintf("Library path: %s\n", lib_path))
