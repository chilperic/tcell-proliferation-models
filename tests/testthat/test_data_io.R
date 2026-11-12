test_that("CFSE files reshape to one total-count matrix", {
  path <- test_path("fixtures", "agent_based_sim_output_1.txt")
  counts <- read_cfse_counts(path, max_generation = 20)
  expect_equal(dim(counts), c(11, 21))
  expect_equal(counts[1, 1], 50000)
  expect_true(all(counts[1, -1] == 0))
  expect_equal(colnames(counts)[1:3], c("N0", "N1", "N2"))
})
