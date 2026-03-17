# Coverage Analysis for DRC Package

## Objective

Identify all functions in the drc package with less than 10% test coverage, sorted by code size (largest functions first).

## Running the Analysis

### Prerequisites

1. R must be installed on your system
2. The `covr` package must be installed:
   ```r
   install.packages("covr")
   ```

### Execute the Analysis

From the root directory of the drc package, run:

```bash
Rscript analyze_coverage.R
```

### Output Files

The script will generate two files:

1. **low_coverage_functions.csv** - CSV format with all low-coverage functions
2. **low_coverage_functions.md** - Markdown table format for easy viewing

### Expected Output Format

The analysis will produce a table with the following columns:

| Column | Description |
|--------|-------------|
| File | Source file containing the function |
| Function | Name of the function |
| Coverage_Pct | Percentage of code covered by tests (< 10%) |
| Code_Lines | Number of lines of code in the function |
| First_Line | Starting line number in the source file |
| Last_Line | Ending line number in the source file |

Results are sorted by `Code_Lines` in descending order (largest functions first).

## Alternative: Manual Analysis

If you cannot run the script, you can analyze coverage manually:

```r
library(covr)

# Run coverage
cov <- covr::package_coverage()

# Convert to data frame
summary_df <- as.data.frame(cov)

# Filter for low coverage
low_cov <- summary_df[summary_df$coverage < 10, ]

# Add code size
low_cov$code_lines <- low_cov$last_line - low_cov$first_line + 1

# Sort by size
result <- low_cov[order(-low_cov$code_lines), ]

# View results
print(result[, c("filename", "functions", "coverage", "code_lines")])
```

## Function Size Analysis

Based on static code analysis, the largest functions in the package are:

| File | Function | Approximate Lines |
|------|----------|-------------------|
| maED.R | maED | 204 |
| arandaordaz.R | arandaordaz | 128 |
| CIcompX.R | CIcompX | 119 |
| EDcomp.R | createsifct | 117 |
| findbe.R | findbe2 | 90 |
| siInner.R | siInner | 88 |
| cedergreen.R | cedergreen_edfct | 70 |
| ED_robust.R | maED_robust | 66 |
| ED_robust.R | ED_robust | 62 |

These large functions are prime candidates for having incomplete test coverage and should be prioritized when reviewing the coverage analysis results.

## Notes

- The analysis excludes functions that cannot be reliably tested (e.g., package initialization)
- Coverage percentages may vary slightly between runs due to test execution order
- Functions with 0% coverage indicate no tests exercise that code path
