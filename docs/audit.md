# Refactor audit

## Original code state

The uploaded archive contained exploratory R scripts, generated synthetic datasets, EPS plots, hidden Git directories, `.RData`, and `.Rhistory`. The code was useful for thesis work, but not ready for direct GitHub publication.

## Problems fixed in the refactor

### 1. Absolute paths

Legacy examples used paths such as:

```r
~/Dropbox/Data/datasets/agent_based_sim_output_3.txt
~/Dropbox/Aims-project/Data/datasets/agent_based_sim_output_1.txt
```

These were replaced by repository-relative paths such as:

```r
data/synthetic_cfse/agent_based_sim_output_3.txt
```

### 2. Global variables

Legacy functions depended on objects such as `ngen`, `N`, `init.vec`, `t_pts`, and `real.data` existing in the global workspace.

Refactored functions take explicit arguments:

```r
simulate_smith_martin(params, times, initial_cells, max_generation)
simulate_cyton(params, times, initial_cells, max_generation)
fit_smith_martin(observed, times, n_starts)
fit_cyton(observed, times, n_starts)
```

### 3. Interactive commands

Commands such as `View()` and repeated plotting blocks were removed from model code. Plotting is now isolated in `plot_generation_counts()`.

### 4. Data reshaping duplication

The repeated manual summation

```r
N0d + N0p, N1d + N1p, ...
```

was replaced by `read_cfse_counts()` and `prepare_cfse_counts()`.

### 5. Mathematical ambiguity in Smith-Martin survival factors

The thesis equations use different survival factors for generation zero and later generations. Some exploratory scripts use `exp(-d_0 * delta)` for all delayed terms. The refactor uses the equation-consistent convention by default and exposes `legacy_survival = TRUE` for reproducibility checks.

### 6. Cyton formula singularity

The analytic Cyton expression divides by `theta0 - theta`. The refactor detects the near-singular case and stops with a useful error instead of silently producing unstable values.

### 7. ARE definition

The refactored `average_relative_error()` uses the target parameter in the denominator. Some exploratory scripts divided by the estimated parameter. That difference should be stated in any reproduction report.

## Remaining work before public release

1. Run all scripts in a local R environment with `deSolve` installed.
2. Compare refactored output against the legacy scripts for several parameter sets.
3. Decide whether the repository should prioritize thesis-equation consistency or exact legacy-script reproduction.
4. Clarify attribution for code comments mentioning external modification.
5. Add generated figures only if they are reproducible from scripts.
6. Consider archiving a release on Zenodo and adding a DOI to `CITATION.cff`.

## Reproducibility environment added

The repository now includes a Docker-first setup. This is the preferred GitHub publication route because it pins the R interpreter through `rocker/r-ver:4.4.3` and installs R dependencies from a fixed CRAN snapshot. Native Ubuntu installation is included for convenience, but it is not equivalent to Docker-level reproducibility because the operating system package manager controls the exact R version.
