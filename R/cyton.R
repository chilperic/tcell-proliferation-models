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
  params[required]
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
#' @return Numeric matrix with generation counts as columns.
simulate_cyton <- function(params,
                           times = seq(0, 120, by = 12),
                           initial_cells = 50000,
                           max_generation = 20,
                           tolerance = 1e-10,
                           clip_negative = TRUE) {
  params <- normalize_cyton_params(params)
  times <- sort(unique(as.numeric(times)))
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
