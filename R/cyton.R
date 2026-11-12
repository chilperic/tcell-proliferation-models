# Cyton model with exponential time-to-division and time-to-death distributions.

normalize_cyton_params <- function(params) {
  param_names <- names(params)
  params <- as.numeric(params)
  names(params) <- param_names
  required <- c("f0", "f1", "mu0d", "mud", "mu0p", "mup")
  missing <- setdiff(required, names(params))
  if (length(missing) > 0L) {
    stop("Missing Cyton parameter(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  params <- params[required]
  if (any(!is.finite(params)) || any(params[c("f0", "f1")] < 0) || any(params[c("f0", "f1")] > 1) || any(params[c("mu0d", "mud", "mu0p", "mup")] <= 0))
    stop("Cyton fractions must be in [0, 1] and mean waiting times positive", call. = FALSE)
  params
}

.cyton_check_singularity <- function(theta0, theta, tolerance) {
  if (abs(theta0 - theta) < tolerance) {
    stop(
      "The analytic Cyton formula is singular or ill-conditioned when theta0 is close to theta. ",
      "Use separated parameters or a numerical convolution implementation.",
      call. = FALSE
    )
  }
}

.cyton_sum_rn_div2 <- function(i, a, n, time, theta) {
  j <- 0:i
  ((-1)^(n + i + 1) / a^(n - i)) *
    (1 / theta^(i + 1) - exp(-theta * time) * sum(time^(i - j) / factorial(i - j) / theta^(j + 1)))
}

.cyton_sum_rn_div <- function(n, time, params, initial_cells, tolerance) {
  theta0 <- 1 / params[["mu0p"]] + 1 / params[["mu0d"]]

  if (n == 0L) {
    return(params[["f0"]] * initial_cells / params[["mu0p"]] / theta0 * (1 - exp(-theta0 * time)))
  }

  theta <- 1 / params[["mup"]] + 1 / params[["mud"]]
  .cyton_check_singularity(theta0, theta, tolerance)
  a <- theta0 - theta

  prefactor <- params[["f0"]] / params[["mu0p"]] *
    (2 * params[["f1"]] / params[["mup"]])^n * initial_cells

  prefactor * (((-1 / a)^n * (1 - exp(-theta0 * time)) / theta0) +
                 sum(vapply(0:(n - 1L), .cyton_sum_rn_div2, numeric(1L), a = a, n = n, time = time, theta = theta)))
}

.cyton_sum_rn_die2 <- function(i, a, n, time, theta) {
  j <- 0:(i + 1L)
  ((-1)^(n + i) / a^(n - i - 1L)) *
    (1 / theta^(i + 2L) - exp(-theta * time) * sum(time^(i - j + 1L) / factorial(i - j + 1L) / theta^(j + 1L)))
}

.cyton_sum_rn_die <- function(n, time, params, initial_cells, tolerance) {
  theta0 <- 1 / params[["mu0p"]] + 1 / params[["mu0d"]]

  if (n == 0L) {
    return(initial_cells / params[["mu0d"]] *
             ((1 - params[["f0"]]) * params[["mu0d"]] * (1 - exp(-time / params[["mu0d"]])) +
                params[["f0"]] / theta0 * (1 - exp(-theta0 * time))))
  }

  theta <- 1 / params[["mup"]] + 1 / params[["mud"]]
  .cyton_check_singularity(theta0, theta, tolerance)
  a <- theta0 - theta

  prefactor <- 2 * initial_cells * params[["f0"]] / params[["mu0p"]] / params[["mud"]] *
    (2 * params[["f1"]] / params[["mup"]])^(n - 1L)

  convolution_sum <- if (n > 1L) {
    sum(vapply(0:(n - 2L), .cyton_sum_rn_die2, numeric(1L), a = a, n = n, time = time, theta = theta))
  } else {
    0
  }

  term <- (-1)^(n - 1L) / a^n * ((1 - exp(-theta * time)) / theta + (exp(-theta0 * time) - 1) / theta0) +
    convolution_sum

  prefactor * term
}

.cyton_generation_count <- function(n, time, params, initial_cells, tolerance) {
  if (n == 0L) {
    initial_cells -
      .cyton_sum_rn_div(0L, time, params, initial_cells, tolerance) -
      .cyton_sum_rn_die(0L, time, params, initial_cells, tolerance)
  } else {
    2 * .cyton_sum_rn_div(n - 1L, time, params, initial_cells, tolerance) -
      .cyton_sum_rn_div(n, time, params, initial_cells, tolerance) -
      .cyton_sum_rn_die(n, time, params, initial_cells, tolerance)
  }
}

#' Simulate the Cyton model under exponential waiting-time assumptions.
#'
#' @param params Named vector with `f0`, `f1`, `mu0d`, `mud`, `mu0p`, and `mup`.
#' @param times Numeric vector of output times.
#' @param initial_cells Initial number of generation-zero cells.
#' @param max_generation Largest generation index.
#' @param tolerance Numerical threshold for detecting the singular analytic case
#'   `theta0 == theta`.
#' @param clip_negative If `TRUE`, very small negative values produced by floating
#'   point cancellation are clipped to zero.
#' @param method ode uses the exponential competing-risk compartment system; legacy uses the historical analytic formula.
#' @param rtol Relative solver tolerance.
#' @param atol Absolute solver tolerance.
#' @return Numeric matrix with generation counts as columns.
simulate_cyton <- function(params,
                           times = seq(0, 120, by = 12),
                           initial_cells = 50000,
                           max_generation = 20,
                           tolerance = 1e-10,
                           clip_negative = TRUE, method = c("ode", "legacy"),
                           rtol = 1e-8, atol = 1e-8) {
  params <- normalize_cyton_params(params)
  times <- .validate_times(times)
  .validate_initial(max_generation, initial_cells)
  method <- match.arg(method)
  if (method == "ode") return(.cyton_ode(params, times, initial_cells, max_generation, rtol, atol))
  g <- as.integer(max_generation)

  values <- vapply(0:g, function(n) {
    vapply(times, function(time) {
      .cyton_generation_count(n, time, params, initial_cells, tolerance)
    }, numeric(1L))
  }, numeric(length(times)))

  values <- as.matrix(values)
  colnames(values) <- paste0("N", 0:g)
  rownames(values) <- as.character(times)

  if (isTRUE(clip_negative)) {
    values[values < 0 & values > -sqrt(.Machine$double.eps)] <- 0
  }

  attr(values, "times") <- times
  values
}


# Exponential waiting times give two Markov compartments per generation:
# committed proliferators P and nondividers Q. Death applies to both.
.cyton_ode <- function(params, times, initial_cells, g, rtol, atol) {
  if (!is.finite(rtol) || rtol <= 0 || !is.finite(atol) || atol <= 0) stop("Positive finite tolerances required")
  width <- g + 1L
  state <- numeric(2L * width)
  state[1L] <- initial_cells * params[["f0"]]
  state[width + 1L] <- initial_cells * (1 - params[["f0"]])
  division <- c(1 / params[["mu0p"]], rep(1 / params[["mup"]], g))
  death <- c(1 / params[["mu0d"]], rep(1 / params[["mud"]], g))
  rhs <- function(t, y, parms) {
    p <- y[seq_len(width)]
    q <- y[width + seq_len(width)]
    births <- c(0, 2 * division[seq_len(g)] * p[seq_len(g)])
    list(c(params[["f1"]] * births - (division + death) * p,
           (1 - params[["f1"]]) * births - death * q))
  }
  out <- deSolve::ode(state, times, rhs, parms = NULL, rtol = rtol, atol = atol)
  counts <- out[, 1L + seq_len(width), drop = FALSE] + out[, 1L + width + seq_len(width), drop = FALSE]
  if (nrow(counts) != length(times) || any(!is.finite(counts)) || any(counts < -10 * atol))
    stop("Incomplete, nonfinite or negative Cyton solution", call. = FALSE)
  colnames(counts) <- paste0("N", 0:g)
  rownames(counts) <- as.character(times)
  attr(counts, "times") <- times
  attr(counts, "method") <- "exponential-competing-risks-ode"
  counts
}
