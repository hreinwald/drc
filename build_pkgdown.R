#!/usr/bin/env Rscript
# Build pkgdown site
# This script safely builds the pkgdown documentation site

# Load pkgdown
if (!requireNamespace("pkgdown", quietly = TRUE)) {
  stop("pkgdown package is required. Install it with: install.packages('pkgdown')")
}

# Check if docs directory exists and is not a pkgdown site
if (dir.exists("docs")) {
  message("Found existing docs/ directory")
  message("Cleaning docs/ directory to ensure pkgdown can build properly...")
  pkgdown::clean_site(force = TRUE)
}

# Build the site
message("Building pkgdown site...")
pkgdown::build_site()

message("Done! Site built successfully in docs/")
