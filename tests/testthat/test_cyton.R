test_that("Cyton simulation preserves initial condition at time zero", {
  params <- c(f0 = 1, f1 = 0.8, mu0d = 430, mud = 23, mu0p = 53, mup = 17)
  out <- simulate_cyton(params, times = c(0, 12), max_generation = 5)
  expect_equal(out[1, "N0"], 50000)
  expect_true(all(out[1, paste0("N", 1:5)] == 0))
})
