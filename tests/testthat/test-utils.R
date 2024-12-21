
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



test_that("test set_exposure", {
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



test_that("test cf_predict", {

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





test_that("test gen_formula", {
  result <- gen_formula(k = 1, interactions_XC = "X1:X2", marginal = TRUE)
  expect_equal(result, "M1~X+X1:X2")
  result <- gen_formula(k = 3, interactions_XC = "X1:X2", marginal = TRUE)
  expect_equal(result, "M3~X+X1:X2")
})

test_that("test gen_formula with include_all = TRUE", {
  result <- gen_formula(k = 3, K = 5, interactions_XC = "X1:X2", include_all = TRUE)
  expect_equal(result, "M3~(X+M1+M2)^2+X1:X2")

  result <- gen_formula(k = 4, first = 2, K = 5, interactions_XC = "X1:X2", include_all = TRUE)
  expect_equal(result, "M4~(X+M2+M3)^2+X1:X2")
})

test_that("test gen_formula for cases involving MM", {
  result <- gen_formula(k = 4, MM = 2, K = 5, interactions_XC = "X1:X2")
  expect_equal(result, "M4~(X+M1+M3)^2+X1:X2")

  result <- gen_formula(k = 4, first = 1, MM = 1, K = 5, interactions_XC = "X1:X2")
  expect_equal(result, "M4~(X+M2+M3)^2+X1:X2")

  result <- gen_formula(k = 2, first = 1, MM = 1, K = 5, interactions_XC = "X1:X2")
  expect_equal(result, "M2~X+X1:X2")
})



test_that("test gen_formula for some special cases", {
  result <- gen_formula(k = 1, interactions_XC = "X1:X2")
  expect_equal(result, "M1~X+X1:X2")

  result <- gen_formula(k = 1, interactions_XC = "X1:X2", include_all = TRUE)
  expect_equal(result, "M1~X+X1:X2")

  result <- gen_formula(k = 4, K = 5, interactions_XC = "X1:X2", include_all = TRUE)
  expect_equal(result, "M4~(X+M1+M2+M3)^2+X1:X2")
})



test_that("test med_joint_other", {

  k <- 4
  a <- 1
  MM <- 2
  K <- 5
  result <- med_joint_other(k = k, a = a, MM = MM, K = K)
  expect_equal(result, "m4_1_101mm")

  result <- med_joint_other(k = k, a = a, MM = MM, K = K, ordering = FALSE)
  expect_equal(result, "m4_1_1m1mm")

  k <- 3
  result <- med_joint_other(k = k, a = a, MM = MM, K = K, ordering = TRUE)
  expect_equal(result, "m3_1_10mmm")

  MM <- 1
  result <- med_joint_other(k = k, a = a, MM = MM, K = K, ordering = FALSE)
  expect_equal(result, "m3_1_m1mmm")


  K <- 1
  expect_error(med_joint_other(k = k, a = a, MM = MM, K = K, ordering = TRUE), "invalid 'times' value")
})




test_that("test med_joint_other", {
  # set mcsim and R to be very small for testing summary function
  set.seed(2024)
  result <- medRCT(
    dat = LSACdata,
    exposure = "sep",
    outcome = "child_mh",
    mediators = c("parent_mh", "preschool_att"),
    intermediate_confs = "fam_stress",
    confounders = c("child_sex", "child_atsi", "mat_cob", "mat_engl", "mat_age"),
    interactions_XC = "all",
    intervention_type = "all",
    mcsim = 10,
    bootstrap = T,
    boot_args = list(R = 10, stype = "i", ci.type = "norm")
  )
  output <- capture.output({
    results <- summary(result)
  })
  expect_true("IIE" %in% names(results))
  expect_true("expected_outcome" %in% names(results))
  expect_equal(results$sample_size, nrow(LSACdata))
  expect_equal(results$n.sim, 10)

  result_non_boot <- medRCT(
    dat = LSACdata,
    exposure = "sep",
    outcome = "child_mh",
    mediators = c("parent_mh", "preschool_att"),
    intermediate_confs = "fam_stress",
    confounders = c("child_sex", "child_atsi", "mat_cob", "mat_engl", "mat_age"),
    interactions_XC = "all",
    intervention_type = "all",
    mcsim = 10,
    bootstrap = F
  )
  output <- capture.output({
    result <- summary.medRCT(result_non_boot)
  })
  expect_false(is.null(names(result)))
  expect_length(result, 11)

  # Ensure no bootstrap-related output appears
  expect_false(any(grepl("Estimated interventional effect:", output)))
  expect_false(any(grepl("Simulations:", output)))
})
