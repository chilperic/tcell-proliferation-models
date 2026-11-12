#!/usr/bin/env Rscript
# Install the R package dependencies needed by this repository.
# Reproducibility policy:
#   - Docker pins the R interpreter via rocker/r-ver.
#   - Package resolution uses the CRAN_REPO snapshot configured in .Rprofile
#     or provided through the CRAN_REPO environment variable.

repo <- Sys.getenv("CRAN_REPO", unset = getOption("repos")[["CRAN"]])
options(repos = c(CRAN = repo))

required <- c("deSolve", "testthat", "jsonlite", "digest")
installed <- rownames(installed.packages())
missing <- setdiff(required, installed)

cat("R version:", paste(R.version$major, R.version$minor, sep = "."), "\n")
cat("CRAN repository:", getOption("repos")[["CRAN"]], "\n")

if (length(missing) > 0L) {
  cat("Installing missing packages:", paste(missing, collapse = ", "), "\n")
  install.packages(missing, dependencies = TRUE)
} else {
  cat("All required packages already installed.\n")
}

cat("\nInstalled versions:\n")
for (pkg in required) {
  cat(sprintf("  %-10s %s\n", pkg, as.character(utils::packageVersion(pkg))))
}
