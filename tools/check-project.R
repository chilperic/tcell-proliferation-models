#!/usr/bin/env Rscript
# Smoke-test the repository after installation.
# This deliberately avoids interactive plotting and long optimization runs.

`%||%` <- function(x, y) if (is.null(x)) y else x
script_file <- sys.frame(1)$ofile %||% "tools/check-project.R"
root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = FALSE)
if (!dir.exists(file.path(root, "R"))) {
  root <- normalizePath(getwd(), mustWork = TRUE)
}
setwd(root)

message("Repository root: ", getwd())
message("R version: ", R.version.string)

source("R/validation.R")
source("R/cfse_data.R")
source("R/metrics.R")
source("R/smith_martin.R")
source("R/cyton.R")
source("R/fit.R")
source("R/plot.R")

observed <- read_cfse_counts(
  "data/synthetic_cfse/agent_based_sim_output_1.txt",
  max_generation = 20
)
times <- attr(observed, "times")
stopifnot(is.matrix(observed), length(times) == nrow(observed))

sm_params <- c(
  lambda_0 = 1.870167e-02,
  lambda   = 3.422586e-02,
  d_0      = 1.000000e-03,
  d        = 3.425349e-02,
  delta    = 1.000000e-03
)
sm <- simulate_smith_martin(sm_params, times = times, max_generation = 20)
stopifnot(is.matrix(sm), all(dim(sm) == dim(observed)))
stopifnot(is.finite(ssd(observed, sm)))

cy_params <- c(
  f0   = 1.000000e+00,
  f1   = 8.115505e-01,
  mu0d = 4.300041e+02,
  mud  = 2.341927e+01,
  mu0p = 5.316612e+01,
  mup  = 1.722716e+01
)
cy <- simulate_cyton(cy_params, times = times, max_generation = 20)
stopifnot(is.matrix(cy), all(dim(cy) == dim(observed)))
stopifnot(is.finite(ssd(observed, cy)))

if (requireNamespace("testthat", quietly = TRUE)) {
  testthat::test_dir("tests/testthat", reporter = "summary")
}

message("Installation and smoke tests passed.")
