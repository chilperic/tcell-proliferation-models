# Model-comparison metrics.

#' Sum of squared differences.
#'
#' @param observed Numeric matrix of observed counts.
#' @param predicted Numeric matrix of predicted counts.
#' @param weights Optional numeric matrix or scalar weights.
#' @return Sum of squared weighted residuals.
ssd <- function(observed, predicted, weights = 1) {
  observed <- as.matrix(observed)
  predicted <- as.matrix(predicted)
  if (!identical(dim(observed), dim(predicted))) {
    stop("Observed and predicted matrices must have identical dimensions.", call. = FALSE)
  }
  if (!length(observed) || any(!is.finite(observed)) || any(!is.finite(predicted)))
    stop("Counts must be nonempty and finite", call. = FALSE)
  if (!(length(weights) == 1L || identical(dim(weights), dim(observed))) || any(!is.finite(weights)) || any(weights < 0))
    stop("Weights must be nonnegative finite scalar or a matrix matching counts", call. = FALSE)
  sum(weights * (observed - predicted)^2)
}

#' Average relative error.
#'
#' This uses the standard denominator `abs(true)`, not the denominator used in
#' some exploratory legacy scripts.
#'
#' @param true Numeric vector of target parameter values.
#' @param estimated Numeric vector or matrix of estimates. If a matrix is given,
#'   rows are interpreted as replicate estimates.
#' @param eps Small denominator guard.
#' @return Average relative error in percent.
average_relative_error <- function(true, estimated, eps = .Machine$double.eps) {
  true <- as.numeric(true)
  estimated <- as.matrix(estimated)

  if (ncol(estimated) != length(true)) {
    if (nrow(estimated) == length(true)) {
      estimated <- t(estimated)
    } else {
      stop("The number of estimated parameters must match `true`.", call. = FALSE)
    }
  }

  if (!length(true) || any(!is.finite(true)) || !length(estimated) || any(!is.finite(estimated)) || length(eps) != 1L || !is.finite(eps) || eps <= 0)
    stop("Finite nonempty parameters and positive finite eps required", call. = FALSE)
  denom <- pmax(abs(true), eps)
  colMeans(abs(sweep(estimated, 2L, true, `-`)) / matrix(denom, nrow = nrow(estimated), ncol = length(true), byrow = TRUE)) * 100
}
