# Reproduce a Cyton simulation from the refactored functions.
# Run from the repository root:
#   Rscript scripts/02_simulate_cyton.R

source("R/cfse_data.R")
source("R/metrics.R")
source("R/cyton.R")
source("R/plot.R")

params <- c(
  f0 = 1.000000e+00,
  f1 = 8.115505e-01,
  mu0d = 4.300041e+02,
  mud = 2.341927e+01,
  mu0p = 5.316612e+01,
  mup = 1.722716e+01
)

times <- seq(0, 120, by = 12)
observed <- read_cfse_counts("data/synthetic_cfse/agent_based_sim_output_1.txt", max_generation = 20)
predicted <- simulate_cyton(params, times = times, max_generation = 20)

print(ssd(observed, predicted))
plot_generation_counts(observed, predicted, generations = 0:8, main = "Cyton gen.")
