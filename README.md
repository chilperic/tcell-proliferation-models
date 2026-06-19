# T-cell proliferation models for CFSE-like data

Refactored research code for simulating and fitting T-cell proliferation models used in the AIMS Ghana Master's thesis:

> Chilperic Armel Foko Kuate. *Reverse-engineering of T cell proliferation dynamic*. AIMS Ghana, 2016.

The repository focuses on two models used in the thesis analysis:

- **Smith-Martin model**: a delayed differential-equation model for generation-structured proliferation.
- **Cyton model**: an analytic implementation under exponential waiting-time assumptions for division and death.

The original scripts were preserved under `legacy/original_R_scripts/`. The refactored code is under `R/`.

## Why this refactor exists

The thesis appendix code was written as exploratory research code. It mixed model definitions, data loading, fitting, plotting, global variables, absolute paths, and interactive commands. That structure is normal during a thesis but weak for GitHub publication.

This version separates:

- model simulation,
- data input and reshaping,
- parameter fitting,
- metrics,
- plotting,
- reproducible scripts.

## Repository layout

```text
tcell-proliferation-models/
├── R/                         # refactored functions
├── scripts/                   # reproducible examples
├── data/synthetic_cfse/        # synthetic CFSE-like datasets from the thesis project
├── legacy/original_R_scripts/  # original unedited R scripts
├── tests/                     # lightweight tests
├── docs/                      # audit and model notes
├── DESCRIPTION
├── NAMESPACE
├── CITATION.cff
└── LICENSE
```

## Installation / setup

For reproducible local execution, use the Docker environment included in this repository:

```bash
docker build -t tcell-proliferation-models:r-4.4.3 .
docker run --rm tcell-proliferation-models:r-4.4.3
```

For native Ubuntu/Debian installation:

```bash
chmod +x tools/install_r_ubuntu_24_04.sh
./tools/install_r_ubuntu_24_04.sh
```

Detailed instructions are in [`INSTALL.md`](INSTALL.md).

To run without installing the package, source the functions from the repository root:

```r
source("R/cfse_data.R")
source("R/metrics.R")
source("R/smith_martin.R")
source("R/cyton.R")
source("R/fit.R")
source("R/plot.R")
```

## Minimal Smith-Martin example

```r
source("R/cfse_data.R")
source("R/metrics.R")
source("R/smith_martin.R")

observed <- read_cfse_counts(
  "data/synthetic_cfse/agent_based_sim_output_3.txt",
  max_generation = 20
)

times <- attr(observed, "times")

params <- c(
  lambda_0 = 1.870167e-02,
  lambda   = 3.422586e-02,
  d_0      = 1.000000e-03,
  d        = 3.425349e-02,
  delta    = 1.000000e-03
)

predicted <- simulate_smith_martin(
  params,
  times = times,
  max_generation = 20
)

ssd(observed, predicted)
```

## Minimal Cyton example

```r
source("R/cfse_data.R")
source("R/metrics.R")
source("R/cyton.R")

observed <- read_cfse_counts(
  "data/synthetic_cfse/agent_based_sim_output_1.txt",
  max_generation = 20
)

times <- attr(observed, "times")

params <- c(
  f0   = 1.000000e+00,
  f1   = 8.115505e-01,
  mu0d = 4.300041e+02,
  mud  = 2.341927e+01,
  mu0p = 5.316612e+01,
  mup  = 1.722716e+01
)

predicted <- simulate_cyton(
  params,
  times = times,
  max_generation = 20
)

ssd(observed, predicted)
```

## Reproducing example scripts

From the repository root:

```bash
Rscript scripts/01_simulate_smith_martin.R
Rscript scripts/02_simulate_cyton.R
Rscript scripts/03_fit_models.R
```

The fitting script uses only a few random starts by default so that it runs quickly. Increase `n_starts` for a serious reproduction.

## Important modelling caveat

The exploratory Smith-Martin scripts used the survival factor `exp(-d_0 * delta)` for all delayed terms. The thesis equations distinguish `exp(-d_0 * delta)` for generation zero and `exp(-d * delta)` for later generations. The refactored function uses the equation-consistent convention by default.

To reproduce the exploratory script convention, call:

```r
simulate_smith_martin(params, legacy_survival = TRUE)
```

Do not hide this distinction in a publication-quality repository. It changes the scientific interpretation.

## Data

The files in `data/synthetic_cfse/` are synthetic CFSE-like outputs generated from an individual-based simulation during the thesis project. Each row corresponds to a sampled time point. Each generation has two columns, such as `N0d` and `N0p`; `read_cfse_counts()` sums them into total generation counts `N0`, `N1`, ... .

## Citation

Use the information in `CITATION.cff`. If the repository is archived on Zenodo, replace the placeholder repository URL and add the DOI.

## License and authorship warning

A MIT license file is included as a practical default. Before publishing, verify that every non-trivial contribution is yours to license. In particular, the legacy Cyton file contains the comment `modified code by Lamin`; clarify attribution before public release.
