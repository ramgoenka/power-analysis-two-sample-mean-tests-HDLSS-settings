###############################################################################
# combine_hotelling.R
###############################################################################

result_files <- list.files("results_hotelling", pattern = "^result_hotelling_.*\\.csv$", full.names = TRUE)

if (length(result_files) == 0) {
  stop("No result files found in results_hotelling/ directory.")
}

cat(sprintf("Found %d result files.\n", length(result_files)))

all_results <- do.call(rbind, lapply(result_files, read.csv))
all_results <- all_results[order(all_results$task_id), ]

write.csv(all_results, "hotelling_results_combined.csv", row.names = FALSE)

cat(sprintf("Combined results saved to hotelling_results_combined.csv\n"))
cat(sprintf("Total configs: %d\n", nrow(all_results)))

# Quick summary
cat("\n=== Type I Error (null signal) ===\n")
null_res <- all_results[all_results$signal == "null", ]
print(null_res[, c("p", "n", "cov_model", "power_HT")])

cat("\n=== Power (sparse signal) ===\n")
sparse_res <- all_results[all_results$signal == "sparse", ]
print(sparse_res[, c("p", "n", "cov_model", "power_HT")])

cat("\n=== Timing Summary ===\n")
cat(sprintf("Total compute time: %.1f minutes\n", sum(all_results$time_min)))
cat(sprintf("Mean per config: %.1f minutes\n", mean(all_results$time_min)))
cat(sprintf("Max per config: %.1f minutes\n", max(all_results$time_min)))
