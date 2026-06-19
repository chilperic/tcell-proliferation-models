# Lightweight base-R plotting helpers.

#' Plot observed and predicted generation counts.
#'
#' @param observed Numeric matrix of observed counts.
#' @param predicted Optional numeric matrix of predicted counts.
#' @param generations Generation indices to plot.
#' @param times Optional time vector.
#' @param main Plot title prefix.
#' @return Invisibly returns `NULL`.
plot_generation_counts <- function(observed,
                                   predicted = NULL,
                                   generations = 0:8,
                                   times = NULL,
                                   main = "Generation") {
  observed <- as.matrix(observed)
  if (!is.null(predicted)) {
    predicted <- as.matrix(predicted)
  }

  times <- .extract_model_times(observed, times)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  n <- length(generations)
  graphics::par(mfrow = c(ceiling(n / 3), min(3, n)), mar = c(4, 4, 3, 1))

  for (g in generations) {
    col <- paste0("N", g)
    if (!col %in% colnames(observed)) {
      next
    }
    y_range <- range(observed[, col], if (!is.null(predicted) && col %in% colnames(predicted)) predicted[, col] else NULL, finite = TRUE)
    graphics::plot(times, observed[, col], type = "o", xlab = "Time", ylab = "Cell count", main = paste(main, g), ylim = y_range)
    if (!is.null(predicted) && col %in% colnames(predicted)) {
      graphics::lines(times, predicted[, col], type = "o", lty = 2)
      graphics::legend("topright", legend = c("Observed", "Predicted"), lty = c(1, 2), bty = "n")
    }
    graphics::grid()
  }
  invisible(NULL)
}
