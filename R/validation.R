# Shared checks. Do not silently sort observations or recycle parameter vectors.
.validate_times <- function(times) {
  times <- as.numeric(times)
  if (length(times) < 2L || any(!is.finite(times)) || times[1L] != 0 || any(diff(times) <= 0))
    stop("times must start at zero and be finite and strictly increasing", call. = FALSE)
  times
}
.validate_initial <- function(max_generation, initial_cells) {
  if (length(max_generation) != 1L || !is.finite(max_generation) || max_generation < 1 || max_generation != as.integer(max_generation))
    stop("max_generation must be an integer >= 1", call. = FALSE)
  if (length(initial_cells) != 1L || !is.finite(initial_cells) || initial_cells < 0)
    stop("initial_cells must be finite and nonnegative", call. = FALSE)
}
.validate_fit <- function(observed, times, lower, upper, n_starts) {
  .validate_times(times)
  if (!is.numeric(observed) || any(!is.finite(observed)) || any(observed < 0) || nrow(observed) != length(times))
    stop("observed must be finite nonnegative counts aligned with times", call. = FALSE)
  if (!identical(names(lower), names(upper)) || any(lower >= upper))
    stop("Named lower bounds must be strictly below upper bounds", call. = FALSE)
  if (length(n_starts) != 1L || !is.finite(n_starts) || n_starts < 1 || n_starts != as.integer(n_starts))
    stop("n_starts must be a positive integer", call. = FALSE)
}
