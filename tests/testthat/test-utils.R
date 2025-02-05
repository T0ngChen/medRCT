
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
    interactions_XC = "none",
    intervention_type = "all",
    mcsim = 10,
    bootstrap = T,
    boot_args = list(R = 30, stype = "i", ci.type = "norm")
  )
  output <- capture.output({
    results <- summary(result)
  })
  expect_true("IIE" %in% names(results))
  expect_true("expected_outcome" %in% names(results))
  expect_equal(results$sample_size, 2608)
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
  expect_length(result, 15)

  # Ensure no bootstrap-related output appears
  expect_false(any(grepl("Estimated interventional effect:", output)))
  expect_false(any(grepl("Simulations:", output)))
})

######################################
# test real examples for gen_formula #
######################################

test_that("test gen_formula", {
  ######################################
  # for estimating joint distributions #
  # k in 1:K                           #
  ######################################
  expect_equal(gen_formula(k = 1, interactions_XC = "C", include_all = TRUE), "M1~X+C")
  expect_equal(gen_formula(k = 2, interactions_XC = "C", include_all = TRUE), "M2~(X+M1)^2+C")
  expect_equal(gen_formula(k = 3, interactions_XC = "C", include_all = TRUE), "M3~(X+M1+M2)^2+C")

  ######################################
  # for estimating Marginals under X=0 #
  # k in first:K                       #
  ######################################
  # if first == 1
  expect_equal(gen_formula(k = 1, interactions_XC = "C", marginal = TRUE), "M1~X+C")
  expect_equal(gen_formula(k = 2, interactions_XC = "C", marginal = TRUE), "M2~X+C")
  expect_equal(gen_formula(k = 3, interactions_XC = "C", marginal = TRUE), "M3~X+C")
  ##########################################################################
  # models for estimating joint of others under X!=0 for p_first, ..., p_K #
  # MM in first:(K-1)                                                      #
  # k in setdiff(MM:K, MM)                                                 #
  ##########################################################################
  # testing if first == 1, K = 2
  expect_equal(gen_formula(k = 2, MM = 1, first = 1, K = 2, interactions_XC = "C"), "M2~X+C")
  # testing if first == 1, K>2
  ## MM == 1
  expect_equal(gen_formula(k = 2, MM = 1, first = 1, K = 5, interactions_XC = "C"), "M2~X+C")
  expect_equal(gen_formula(k = 3, MM = 1, first = 1, K = 5, interactions_XC = "C"), "M3~(X+M2)^2+C")
  expect_equal(gen_formula(k = 4, MM = 1, first = 1, K = 5, interactions_XC = "C"), "M4~(X+M2+M3)^2+C")
  expect_equal(gen_formula(k = 5, MM = 1, first = 1, K = 5, interactions_XC = "C"), "M5~(X+M2+M3+M4)^2+C")
  ## MM == 2
  expect_equal(gen_formula(k = 3, MM = 2, first = 1, K = 5, interactions_XC = "C"), "M3~(X+M1)^2+C")
  expect_equal(gen_formula(k = 4, MM = 2, first = 1, K = 5, interactions_XC = "C"), "M4~(X+M1+M3)^2+C")
  expect_equal(gen_formula(k = 5, MM = 2, first = 1, K = 5, interactions_XC = "C"), "M5~(X+M1+M3+M4)^2+C")
  ## MM == 3
  expect_equal(gen_formula(k = 4, MM = 3, first = 1, K = 5, interactions_XC = "C"), "M4~(X+M1+M2)^2+C")
  expect_equal(gen_formula(k = 5, MM = 3, first = 1, K = 5, interactions_XC = "C"), "M5~(X+M1+M2+M4)^2+C")
  ## MM == 4
  expect_equal(gen_formula(k = 5, MM = 4, first = 1, K = 5, interactions_XC = "C"), "M5~(X+M1+M2+M3)^2+C")

  # if first != 1 & first == K-1
  expect_equal(gen_formula(k = 4, MM = 3, first = 3, K = 4, interactions_XC = "C"), "M4~(X+M1+M2)^2+C")

  # in general case that first != 0 & K-1>first
  ## MM == 3
  expect_equal(gen_formula(k = 4, MM = 3, first = 3, K = 5, interactions_XC = "C"), "M4~(X+M1+M2)^2+C")
  expect_equal(gen_formula(k = 5, MM = 3, first = 3, K = 5, interactions_XC = "C"), "M5~(X+M1+M2+M4)^2+C")
  ## MM == 4
  expect_equal(gen_formula(k = 5, MM = 3, first = 3, K = 5, interactions_XC = "C"), "M5~(X+M1+M2+M4)^2+C")

  ###################################################################################
  # models for estimating Conditionals under X!=0 for p_first_prime,...., p_K_prime #
  # MM in first:(K-1)                                                               #
  # k in (MM + 1):K                                                                 #
  ###################################################################################

  # testing if first == 1, K = 2
  expect_equal(gen_formula(k = 2, K = 2, interactions_XC = "C", include_all = TRUE), "M2~(X+M1)^2+C")
  # testing if first == 1, K>2
  ## MM == 1
  expect_equal(gen_formula(k = 2, K = 5, interactions_XC = "C", include_all = TRUE), "M2~(X+M1)^2+C")
  expect_equal(gen_formula(k = 3, K = 5, interactions_XC = "C", include_all = TRUE), "M3~(X+M1+M2)^2+C")
  expect_equal(gen_formula(k = 4, K = 5, interactions_XC = "C", include_all = TRUE), "M4~(X+M1+M2+M3)^2+C")
  expect_equal(gen_formula(k = 5, K = 5, interactions_XC = "C", include_all = TRUE), "M5~(X+M1+M2+M3+M4)^2+C")
  ## MM == 2
  expect_equal(gen_formula(k = 3, K = 5, interactions_XC = "C", include_all = TRUE), "M3~(X+M1+M2)^2+C")
  expect_equal(gen_formula(k = 4, K = 5, interactions_XC = "C", include_all = TRUE), "M4~(X+M1+M2+M3)^2+C")
  expect_equal(gen_formula(k = 5, K = 5, interactions_XC = "C", include_all = TRUE), "M5~(X+M1+M2+M3+M4)^2+C")
  ## MM == 3
  expect_equal(gen_formula(k = 4, K = 5, interactions_XC = "C", include_all = TRUE), "M4~(X+M1+M2+M3)^2+C")
  expect_equal(gen_formula(k = 5, K = 5, interactions_XC = "C", include_all = TRUE), "M5~(X+M1+M2+M3+M4)^2+C")
  ## MM == 4
  expect_equal(gen_formula(k = 5, K = 5, interactions_XC = "C", include_all = TRUE), "M5~(X+M1+M2+M3+M4)^2+C")

  # if first != 1 & first == K-1
  ## MM == 3
  expect_equal(gen_formula(k = 4, K = 4, interactions_XC = "C", include_all = TRUE), "M4~(X+M1+M2+M3)^2+C")

  # in general case that first != 0 & K-1>first
  ## MM == 3
  expect_equal(gen_formula(k = 4, K = 5, interactions_XC = "C", include_all = TRUE), "M4~(X+M1+M2+M3)^2+C")
  expect_equal(gen_formula(k = 5, K = 5, interactions_XC = "C", include_all = TRUE), "M5~(X+M1+M2+M3+M4)^2+C")
  ## MM == 4
  expect_equal(gen_formula(k = 5, K = 5, interactions_XC = "C", include_all = TRUE), "M5~(X+M1+M2+M3+M4)^2+C")

  ################################################################
  # models for estimating Joint of main ones under X=0 for p_all #
  # k in (first + 1):K                                           #
  ################################################################
  ## first == 2 && K == 5
  expect_equal(gen_formula(k = 3, interactions_XC = "C", include_all = TRUE, first = 2), "M3~(X+M2)^2+C")
  expect_equal(gen_formula(k = 4, interactions_XC = "C", include_all = TRUE, first = 2), "M4~(X+M2+M3)^2+C")
  expect_equal(gen_formula(k = 5, interactions_XC = "C", include_all = TRUE, first = 2), "M5~(X+M2+M3+M4)^2+C")
})




