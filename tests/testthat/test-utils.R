
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
  expect_equal(levels(result$exposure), "2")
})


test_that("test med_outcome_name", {
  expect_equal(
    med_outcome_name(l = 1, a = 0, K = 3),
    "m1_0_mmm"
  )

  expect_equal(
    med_outcome_name(l = 2, a = "1", K = 4),
    "m2_1_1mmm"
  )

  expect_equal(
    med_outcome_name(l = 1, a = 1, K = 1),
    "m1_1_m"
  )

  expect_equal(
    med_outcome_name(l = 3, a = 3, K = 5),
    "m3_3_33mmm"
  )
})



test_that("test med_outcome_all", {
  expect_equal(
    med_outcome_all(l = 2, first = 1, a = 1, K = 3),
    "m2_1_1mm"
  )

  expect_equal(
    med_outcome_all(l = 3, first = 2, a = 1, K = 4),
    "m3_1_m1mm"
  )

  expect_equal(
    med_outcome_all(l = 2, first = 2, a = 1, K = 4),
    "m2_1_mmmm"
  )

  expect_equal(
    med_outcome_all(l = 1, first = 1, a = 2, K = 1),
    "m1_2_m"
  )
})



test_that("cf_predict generates counterfactual predictions correctly", {

  data <- data.table::data.table(x = rnorm(100))
  fit_binomial <- glm(rbinom(100, 1, 0.1) ~ x, family = "binomial", data = data)

  data <- cf_predict(fit = fit_binomial, data = data, var_name = "cf_binom", n = 100, family = "binomial")

  expect_true("cf_binom" %in% colnames(data))
  expect_true(all(data$cf_binomial %in% c(0, 1)))


  fit_gaussian <- lm(rnorm(100) ~ x, data = data)

  data <- cf_predict(fit = fit_gaussian, data = data, var_name = "cf_gaussian", n = 100, family = "gaussian")

  expect_true("cf_gaussian" %in% colnames(data))
  expect_true(is.numeric(data$cf_gaussian))

})




