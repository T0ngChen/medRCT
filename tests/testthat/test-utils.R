
test_that("check faimly type", {
  n <- 100
  df = data.frame(x = sample(c(0,1), n, replace=T))
  df$z = sample(1:8, n, replace = T)
  df$y <- df$x + rnorm(n)
  expect_equal(family_type(df, 'x'), list(stats::binomial()))
  expect_equal(family_type(df, 'y'), list(stats::gaussian()))
  expect_error(family_type(df, 'z', unique_threshold = 10),
               "Error: The variable must be either continuous or binary."
  )
})



test_that("set_exposure assigns values correctly", {
  dt <- data.table::data.table(id = 1:100, exposure = as.factor(sample(c(0,1,2), 100, replace=T)))
  fit = lm(rnorm(100)~dt$exposure)

  # Apply the function
  result <- set_exposure(dt, "exposure", 2)

  # Check that the column has been updated correctly
  expect_equal(unique(as.numeric(result$exposure)), 1)
})



