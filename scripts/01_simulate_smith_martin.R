# Reproduce a Smith-Martin simulation from the refactored functions.
# Run from the repository root:
#   Rscript scripts/01_simulate_smith_martin.R

source("R/cfse_data.R")
source("R/metrics.R")
source("R/smith_martin.R")
source("R/plot.R")

params <- c(
  lambda_0 = 1.870167e-02,
  lambda = 3.422586e-02,
  d_0 = 1.000000e-03,
  d = 3.425349e-02,
  delta = 1.000000e-03
)

times <- seq(0, 120, by = 12)
observed <- read_cfse_counts("data/synthetic_cfse/agent_based_sim_output_3.txt", max_generation = 20)
predicted <- simulate_smith_martin(params, times = times, max_generation = 20)

print(ssd(observed, predicted))
plot_generation_counts(observed, predicted, generations = 0:8, main = "Smith-Martin gen.")
