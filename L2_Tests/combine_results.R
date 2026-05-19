###############################################################################
# combine_results.R
###############################################################################

result_files <- list.files("results", pattern = "^result_task_.*\\.csv$",
                          full.names = TRUE)

if (length(result_files) == 0) {
  stop("No result files found in results/ directory.")
}

cat(sprintf("Found %d result files.\n", length(result_files)))

all_results <- do.call(rbind, lapply(result_files, read.csv))
all_results <- all_results[order(all_results$task_id), ]

write.csv(all_results, "all_results_combined.csv", row.names = FALSE)
cat(sprintf("Combined results saved to all_results_combined.csv\n"))
cat(sprintf("Total configs: %d\n", nrow(all_results)))

na_cols   <- grep("^na_", names(all_results), value = TRUE)
total_nas <- sum(all_results[, na_cols])
if (total_nas > 0) {
  cat(sprintf("\n[WARNING] %d test calls returned NA across all configs.\n",
              total_nas))
  cat("Inspect the na_* columns to identify which tests/configs failed.\n")
} else {
  cat("\nNo NA p-values: all tests ran cleanly.\n")
}

cat("\n=== Type I Error (signal = null) ===\n")
null_res   <- all_results[all_results$signal == "null", ]
power_cols <- grep("^power_", names(all_results), value = TRUE)
print(null_res[, c("p", "n", "cov_model", "rho", power_cols)])

cat("\n=== Timing Summary ===\n")
cat(sprintf("Total compute (sum of tasks): %.1f hours\n",
            sum(all_results$time_min) / 60))
cat(sprintf("Mean per config: %.1f min\n", mean(all_results$time_min)))
cat(sprintf("Max per config:  %.1f min\n", max(all_results$time_min)))
