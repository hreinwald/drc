#!/usr/bin/env Rscript
#
# Script to identify functions with less than 10% coverage in the drc package
# Sorted by code size (largest on top)
#
# Usage: Rscript analyze_coverage.R
#
# Requirements: covr package must be installed
# Install with: install.packages("covr")

# Check if covr is installed
if (!requireNamespace("covr", quietly = TRUE)) {
  stop("Package 'covr' is required. Install it with: install.packages('covr')")
}

library(covr)

cat("Running coverage analysis on drc package...\n")
cat("This may take several minutes...\n\n")

# Run package coverage
cov <- covr::package_coverage(quiet = FALSE)

# Print overall coverage
cat("\n=== OVERALL COVERAGE ===\n")
print(cov)

# Convert to data frame for analysis
summary_df <- as.data.frame(cov)

# Filter functions with less than 10% coverage
low_cov <- summary_df[summary_df$coverage < 10, ]

if (nrow(low_cov) == 0) {
  cat("\nNo functions found with less than 10% coverage!\n")
  quit(save = "no", status = 0)
}

# Calculate code size (number of lines) for each function
low_cov$code_lines <- low_cov$last_line - low_cov$first_line + 1

# Sort by code size (largest on top)
low_cov_sorted <- low_cov[order(-low_cov$code_lines), ]

# Select and rename relevant columns for final output
result <- data.frame(
  File = low_cov_sorted$filename,
  Function = low_cov_sorted$functions,
  Coverage_Pct = round(low_cov_sorted$coverage, 2),
  Code_Lines = low_cov_sorted$code_lines,
  First_Line = low_cov_sorted$first_line,
  Last_Line = low_cov_sorted$last_line,
  stringsAsFactors = FALSE
)

# Print results
cat("\n\n=== FUNCTIONS WITH LESS THAN 10% COVERAGE ===\n")
cat("(Sorted by code size, largest on top)\n\n")
print(result, row.names = FALSE)

cat("\n\nTotal functions with <10% coverage:", nrow(result), "\n")
cat("Total functions analyzed:", nrow(summary_df), "\n")

# Save to CSV file
output_file <- "low_coverage_functions.csv"
write.csv(result, output_file, row.names = FALSE)
cat("\nResults saved to:", output_file, "\n")

# Also create a markdown table
md_file <- "low_coverage_functions.md"
cat("# Functions with Less Than 10% Coverage\n\n", file = md_file)
cat("Sorted by code size (largest functions first)\n\n", file = md_file, append = TRUE)
cat("| File | Function | Coverage % | Code Lines | Line Range |\n", file = md_file, append = TRUE)
cat("|------|----------|------------|------------|------------|\n", file = md_file, append = TRUE)

for (i in 1:nrow(result)) {
  cat(sprintf("| %s | %s | %.2f%% | %d | %d-%d |\n",
              result$File[i],
              result$Function[i],
              result$Coverage_Pct[i],
              result$Code_Lines[i],
              result$First_Line[i],
              result$Last_Line[i]),
      file = md_file, append = TRUE)
}

cat("\nMarkdown table saved to:", md_file, "\n")
