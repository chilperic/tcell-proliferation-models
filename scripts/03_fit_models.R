# Fit Smith-Martin and Cyton models to the average of three synthetic datasets.
# This is intentionally small. Increase `n_starts` for a serious reproduction.

source("R/cfse_data.R")
source("R/metrics.R")
source("R/smith_martin.R")
source("R/cyton.R")
source("R/fit.R")

paths <- file.path("data/synthetic_cfse", paste0("agent_based_sim_output_", c(1, 3, 4), ".txt"))
observed <- read_average_cfse_counts(paths, max_generation = 20)
times <- attr(observed, "times")

set.seed(1)
sm_fit <- fit_smith_martin(observed, times = times, n_starts = 3, seed = 1)
cyton_fit <- fit_cyton(observed, times = times, n_starts = 3, seed = 1)

print(sm_fit[c("model", "value", "par", "convergence")])
print(cyton_fit[c("model", "value", "par", "convergence")])
