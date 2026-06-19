test_that("Smith-Martin simulation returns generation-count matrix", {
  skip_if_not_installed("deSolve")
  params <- c(lambda_0 = 0.018, lambda = 0.034, d_0 = 0.001, d = 0.034, delta = 0.001)
  out <- simulate_smith_martin(params, times = c(0, 12), max_generation = 5)
  expect_equal(dim(out), c(2, 6))
  expect_equal(out[1, "N0"], 50000)
})
