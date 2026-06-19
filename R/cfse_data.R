# Utilities for reading and reshaping CFSE-like count data.

#' Read one CFSE-like simulated dataset.
#'
#' The original simulation files contain two columns per generation, for
#' example `N0d`, `N0p`, `N1d`, `N1p`, ..., `N20d`, `N20p`. Some files also
#' contain a `Time` column. This function sums the two generation-specific
#' columns and prepends the initial condition at time zero.
#'
#' @param path Path to a tabular text file.
#' @param initial_cells Initial number of cells at time zero.
#' @param times Optional time vector for rows in `path`. If the file has a
#'   `Time` column, that column is used. If neither is available, the rows are
#'   assumed to correspond to 12, 24, ..., hours.
#' @param max_generation Optional largest generation to keep.
#' @return A numeric matrix with rows as time points and columns `N0`, `N1`, ... .
read_cfse_counts <- function(path, initial_cells = 50000, times = NULL, max_generation = NULL) {
  stopifnot(length(path) == 1L)
  if (!file.exists(path)) {
    stop("File not found: ", path, call. = FALSE)
  }

  raw <- utils::read.table(path, header = TRUE, check.names = FALSE)
  prepare_cfse_counts(
    raw,
    initial_cells = initial_cells,
    times = times,
    max_generation = max_generation
  )
}

#' Reshape CFSE-like data already loaded in memory.
#'
#' @param data Data frame with generation columns such as `N0d`, `N0p`, ... .
#' @inheritParams read_cfse_counts
#' @return A numeric matrix with a `times` attribute.
prepare_cfse_counts <- function(data, initial_cells = 50000, times = NULL, max_generation = NULL) {
  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }

  if ("Time" %in% names(data)) {
    file_times <- data[["Time"]]
    data[["Time"]] <- NULL
    if (is.null(times)) {
      times <- file_times
    }
  }

  generation_ids <- sort(unique(as.integer(sub("^N([0-9]+)[dp]$", "\\1", grep("^N[0-9]+[dp]$", names(data), value = TRUE)))))
  generation_ids <- generation_ids[!is.na(generation_ids)]

  if (length(generation_ids) == 0L) {
    stop("No generation columns found. Expected names like N0d, N0p, N1d, N1p, ...", call. = FALSE)
  }

  if (!is.null(max_generation)) {
    generation_ids <- generation_ids[generation_ids <= max_generation]
  }

  counts <- vapply(generation_ids, function(g) {
    death_col <- paste0("N", g, "d")
    proliferating_col <- paste0("N", g, "p")
    missing <- setdiff(c(death_col, proliferating_col), names(data))
    if (length(missing) > 0L) {
      stop("Missing column(s): ", paste(missing, collapse = ", "), call. = FALSE)
    }
    as.numeric(data[[death_col]]) + as.numeric(data[[proliferating_col]])
  }, numeric(nrow(data)))

  counts <- as.matrix(counts)
  colnames(counts) <- paste0("N", generation_ids)

  initial <- numeric(ncol(counts))
  initial[1] <- initial_cells
  counts <- rbind(initial, counts)
  colnames(counts) <- paste0("N", generation_ids)

  if (is.null(times)) {
    times <- seq(12, by = 12, length.out = nrow(data))
  }
  times <- c(0, as.numeric(times))

  if (length(times) != nrow(counts)) {
    stop("The time vector length must match the number of data rows after prepending time zero.", call. = FALSE)
  }

  rownames(counts) <- as.character(times)
  attr(counts, "times") <- times
  counts
}

#' Read and average several CFSE-like datasets.
#'
#' @inheritParams read_cfse_counts
#' @param paths Character vector of file paths.
#' @return Numeric matrix containing element-wise average counts.
read_average_cfse_counts <- function(paths, initial_cells = 50000, max_generation = NULL) {
  if (length(paths) == 0L) {
    stop("At least one path is required.", call. = FALSE)
  }
  matrices <- lapply(paths, read_cfse_counts, initial_cells = initial_cells, max_generation = max_generation)
  template_dim <- dim(matrices[[1L]])
  ok <- vapply(matrices, function(x) identical(dim(x), template_dim), logical(1L))
  if (!all(ok)) {
    stop("All datasets must reshape to the same dimensions before averaging.", call. = FALSE)
  }
  out <- Reduce(`+`, matrices) / length(matrices)
  attr(out, "times") <- attr(matrices[[1L]], "times")
  out
}
