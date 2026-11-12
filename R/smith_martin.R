# Smith-Martin delayed differential-equation model.

#' Build the Smith-Martin initial state vector.
#'
#' The state vector is ordered as A-phase counts `NA0`, ..., `NAg`, followed by
#' total generation counts `N0`, ..., `Ng`.
#'
#' @param max_generation Largest generation index.
#' @param initial_cells Initial number of generation-zero cells.
#' @return Named numeric vector.
smith_martin_initial_state <- function(max_generation = 20, initial_cells = 50000) {
  .validate_initial(max_generation, initial_cells)
  state <- numeric(2L * (max_generation + 1L))
  names(state) <- c(paste0("NA", 0:max_generation), paste0("N", 0:max_generation))
  state[c("NA0", "N0")] <- initial_cells
  state
}

normalize_smith_martin_params <- function(params) {
  param_names <- names(params)
  params <- as.numeric(params)
  names(params) <- param_names
  if ("deltat" %in% names(params) && !("delta" %in% names(params))) {
    names(params)[names(params) == "deltat"] <- "delta"
  }
  required <- c("lambda_0", "lambda", "d_0", "d", "delta")
  missing <- setdiff(required, names(params))
  if (length(missing) > 0L) {
    stop("Missing Smith-Martin parameter(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  params <- params[required]
  if (any(!is.finite(params)) || any(params < 0) || params[["delta"]] <= 0)
    stop("Smith-Martin rates must be nonnegative and delta positive", call. = FALSE)
  params
}

#' Simulate the Smith-Martin DDE model.
#'
#' @param params Named vector with `lambda_0`, `lambda`, `d_0`, `d`, and either
#'   `delta` or legacy name `deltat`.
#' @param times Numeric vector of output times.
#' @param initial_cells Initial number of generation-zero cells.
#' @param max_generation Largest generation index.
#' @param legacy_survival If `TRUE`, reproduce the survival-factor convention
#'   used in the exploratory thesis scripts, where `exp(-d_0 * delta)` is used
#'   for all delayed terms. If `FALSE`, use `exp(-d_0 * delta)` for generation
#'   zero and `exp(-d * delta)` for later generations, matching the model equations.
#' @param history Pre-zero A-phase history: constant (legacy) or zero (synchronized initial cells).
#' @param rtol Relative solver tolerance.
#' @param atol Absolute solver tolerance.
#' @return Numeric matrix with generation counts as columns.
simulate_smith_martin <- function(params,
                                  times = seq(0, 120, by = 12),
                                  initial_cells = 50000,
                                  max_generation = 20,
                                  legacy_survival = FALSE,
                                  history = c("constant", "zero"), rtol = 1e-8, atol = 1e-8) {
  if (!requireNamespace("deSolve", quietly = TRUE)) {
    stop("Package `deSolve` is required for simulate_smith_martin().", call. = FALSE)
  }

  params <- normalize_smith_martin_params(params)
  times <- .validate_times(times)
  history <- match.arg(history)
  if (!is.finite(rtol) || rtol <= 0 || !is.finite(atol) || atol <= 0) stop("Positive finite tolerances required")
  state0 <- smith_martin_initial_state(max_generation, initial_cells)
  g <- as.integer(max_generation)

  rhs <- function(t, y, parms) {
    lambda_0 <- parms[["lambda_0"]]
    lambda <- parms[["lambda"]]
    d_0 <- parms[["d_0"]]
    d <- parms[["d"]]
    delta <- parms[["delta"]]

    lag_a <- if (t < delta) {
      if (history == "constant") state0[seq_len(g + 1L)] else numeric(g + 1L)
    } else {
      deSolve::lagvalue(t - delta)[seq_len(g + 1L)]
    }

    survival_0 <- exp(-d_0 * delta)
    survival <- if (isTRUE(legacy_survival)) survival_0 else exp(-d * delta)

    dy <- numeric(length(y))

    # A-phase compartments.
    dy[1L] <- -(lambda_0 + d_0) * y[1L]
    dy[2L] <- 2 * lambda_0 * survival_0 * lag_a[1L] - (lambda + d) * y[2L]
    if (g >= 2L) {
      dy[3L:(g + 1L)] <- 2 * lambda * survival * lag_a[2L:g] - (lambda + d) * y[3L:(g + 1L)]
    }

    # Total generation counts.
    offset <- g + 1L
    dy[offset + 1L] <- -lambda_0 * survival_0 * lag_a[1L] - d_0 * y[offset + 1L]
    dy[offset + 2L] <- 2 * lambda_0 * survival_0 * lag_a[1L] - lambda * survival * lag_a[2L] - d * y[offset + 2L]
    if (g >= 2L) {
      dy[(offset + 3L):(2L * offset)] <-
        2 * lambda * survival * lag_a[2L:g] -
        lambda * survival * lag_a[3L:(g + 1L)] -
        d * y[(offset + 3L):(2L * offset)]
    }

    list(dy)
  }

  out <- deSolve::dede(y = state0, times = times, func = rhs, parms = params, rtol = rtol, atol = atol, hmax = params[["delta"]])
  counts <- as.matrix(out[, paste0("N", 0:g), drop = FALSE])
  rownames(counts) <- as.character(out[, "time"])
  attr(counts, "times") <- out[, "time"]
  if (nrow(counts) != length(times) || any(!is.finite(counts))) stop("Incomplete or nonfinite DDE solution")
  attr(counts, "history") <- history
  attr(counts, "legacy_survival") <- legacy_survival
  counts
}
