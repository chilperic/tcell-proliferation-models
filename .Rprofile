# Use a fixed CRAN snapshot by default for reproducible local installs.
# Override with: CRAN_REPO=<repo-url> Rscript tools/setup-r-packages.R
local({
  repo <- Sys.getenv("CRAN_REPO", unset = "https://packagemanager.posit.co/cran/2025-04-01")
  options(repos = c(CRAN = repo))
})
