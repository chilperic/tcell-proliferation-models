test_that("Cyton handles equal rates and matches the pure division solution", {
  params <- c(f0=1, f1=1, mu0d=1e15, mud=1e15, mu0p=10, mup=10)
  times <- c(0,1,5,10)
  out <- simulate_cyton(params,times=times,initial_cells=100,max_generation=12)
  for (g in 0:12) {
    expected <- 100*exp(-times/10)*(2*times/10)^g/factorial(g)
    expect_equal(unname(out[,g+1]), expected, tolerance=1e-6)
  }
})
test_that("Cyton nondividers only die and no spontaneous generations appear", {
  params <- c(f0=0, f1=.8, mu0d=20, mud=20, mu0p=10, mup=10)
  out <- simulate_cyton(params,times=c(0,5,10),initial_cells=100,max_generation=3)
  expect_equal(unname(out[,1]),100*exp(-c(0,5,10)/20),tolerance=1e-6)
  expect_true(all(out[,-1] == 0))
})
test_that("Smith-Martin zero history prevents division before the delay", {
  p <- c(lambda_0=.1,lambda=.1,d_0=0,d=0,delta=10)
  out <- simulate_smith_martin(p,times=c(0,1,5),initial_cells=100,max_generation=3,history="zero")
  expect_equal(unname(out[,1]),rep(100,3))
  expect_equal(sum(out[,-1]),0)
})
test_that("time order, fractions, missing counts and weight recycling are rejected", {
  p <- c(f0=.5,f1=.5,mu0d=20,mud=20,mu0p=10,mup=10)
  expect_error(simulate_cyton(p,times=c(0,2,1)),"increasing")
  p['f0'] <- 2
  expect_error(simulate_cyton(p),"fractions")
  expect_error(ssd(matrix(NA_real_),matrix(1)),"finite")
  expect_error(ssd(matrix(1,2,2),matrix(2,2,2),c(1,2)),"Weights")
  expect_error(prepare_cfse_counts(data.frame(N1d=0,N1p=1)),"starting at zero")
})
test_that("replicate averaging checks time coordinates", {
  a <- tempfile(); b <- tempfile()
  on.exit(unlink(c(a,b)))
  write.table(data.frame(Time=12,N0d=0,N0p=2),a,row.names=FALSE)
  write.table(data.frame(Time=24,N0d=0,N0p=2),b,row.names=FALSE)
  expect_error(read_average_cfse_counts(c(a,b)),"same generations and times")
})
