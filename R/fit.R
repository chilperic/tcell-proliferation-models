# Generic fitting helpers for Smith-Martin and Cyton models.

.extract_model_times <- function(observed, times = NULL) {
  if (!is.null(times)) {
    return(as.numeric(times))
  }
  attr_times <- attr(observed, "times")
  if (!is.null(attr_times)) {
    return(as.numeric(attr_times))
  }
  row_numbers <- suppressWarnings(as.numeric(rownames(observed)))
  if (!anyNA(row_numbers) && length(row_numbers) == nrow(observed)) {
    return(row_numbers)
  }
  seq_len(nrow(observed)) - 1L
}

.random_start <- function(lower, upper) {
  stats::runif(length(lower), min = lower, max = upper)
}

#' Fit the Smith-Martin model to generation-count data.
#'
#' @param observed Numeric matrix of generation counts.
#' @param times Optional output times. If omitted, a `times` attribute or row
#'   names are used.
#' @param n_starts Number of random starts for box-constrained optimization.
#' @param lower Named lower bounds.
#' @param upper Named upper bounds.
#' @param seed Optional random seed.
#' @param initial_cells Initial number of cells.
#' @param max_generation Largest generation index.
#' @param legacy_survival Passed to `simulate_smith_martin()`.
#' @return A list containing best parameters, objective value, and all starts.
fit_smith_martin <- function(observed,
                             times = NULL,
                             n_starts = 20,
                             lower = c(lambda_0 = 1e-3, lambda = 1e-3, d_0 = 1e-3, d = 1e-3, delta = 1e-3),
                             upper = c(lambda_0 = 1, lambda = 1, d_0 = 1, d = 1, delta = 100),
                             seed = NULL,
                             initial_cells = observed[1L, 1L],
                             max_generation = ncol(observed) - 1L,
                             legacy_survival = FALSE) {
  observed <- as.matrix(observed)
  times <- .extract_model_times(observed, times)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  lower <- normalize_smith_martin_params(lower)
  upper <- normalize_smith_martin_params(upper)

  objective <- function(par) {
    names(par) <- names(lower)
    predicted <- tryCatch(
      simulate_smith_martin(
        par,
        times = times,
        initial_cells = initial_cells,
        max_generation = max_generation,
        legacy_survival = legacy_survival
      ),
      error = function(e) NULL
    )
    if (is.null(predicted) || any(!is.finite(predicted))) {
      return(Inf)
    }
    ssd(observed, predicted)
  }

  starts <- replicate(n_starts, .random_start(lower, upper), simplify = FALSE)
  fits <- lapply(starts, function(start) {
    names(start) <- names(lower)
    stats::optim(start, objective, method = "L-BFGS-B", lower = lower, upper = upper)
  })

  values <- vapply(fits, `[[`, numeric(1L), "value")
  best_index <- which.min(values)

  list(
    model = "smith_martin",
    value = values[[best_index]],
    par = fits[[best_index]]$par,
    convergence = fits[[best_index]]$convergence,
    starts = fits
  )
}

#' Fit the Cyton model to generation-count data.
#'
#' @inheritParams fit_smith_martin
#' @param lower Named lower bounds for Cyton parameters.
#' @param upper Named upper bounds for Cyton parameters.
#' @return A list containing best parameters, objective value, and all starts.
fit_cyton <- function(observed,
                      times = NULL,
                      n_starts = 20,
                      lower = c(f0 = 1e-3, f1 = 1e-3, mu0d = 1e-3, mud = 1e-3, mu0p = 1e-3, mup = 1e-3),
                      upper = c(f0 = 1, f1 = 1, mu0d = 1000, mud = 1000, mu0p = 1000, mup = 1000),
                      seed = NULL,
                      initial_cells = observed[1L, 1L],
                      max_generation = ncol(observed) - 1L) {
  observed <- as.matrix(observed)
  times <- .extract_model_times(observed, times)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  lower <- normalize_cyton_params(lower)
  upper <- normalize_cyton_params(upper)

  objective <- function(par) {
    names(par) <- names(lower)
    predicted <- tryCatch(
      simulate_cyton(
        par,
        times = times,
        initial_cells = initial_cells,
        max_generation = max_generation
      ),
      error = function(e) NULL
    )
    if (is.null(predicted) || any(!is.finite(predicted))) {
      return(Inf)
    }
    ssd(observed, predicted)
  }

  starts <- replicate(n_starts, .random_start(lower, upper), simplify = FALSE)
  fits <- lapply(starts, function(start) {
    names(start) <- names(lower)
    stats::optim(start, objective, method = "L-BFGS-B", lower = lower, upper = upper)
  })

  values <- vapply(fits, `[[`, numeric(1L), "value")
  best_index <- which.min(values)

  list(
    model = "cyton",
    value = values[[best_index]],
    par = fits[[best_index]]$par,
    convergence = fits[[best_index]]$convergence,
    starts = fits
  )
}
