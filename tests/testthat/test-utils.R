
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


###########################
# test different examples #
###########################

test_that("setting with first = 3, K = 6", {
  first = 3
  K = 6
  intervention_type = "all"
  # ESTIMATE DISTRIBUTIONS
  # Joint of M1 to MK under X=0 and X!=0 ...
  # testing function joint_dist

  joint_dist = list()
  for (k in 1:K) {
    joint_dist[[k]] = character(0)
    a = 1
    # pick up previous draws for prediction
    if (k != 1) {
      l <- 1:(k - 1)
      joint_dist[[k]] = med_outcome_name(a = a, l = l, K = K)
    }
    # random draws for the current
    joint_dist[[k]] = c(joint_dist[[k]], med_outcome_name(a = a, l = k, K = K))
  }
  expect_list = list(
    c("m1_1_mmmmmm"),
    c("m1_1_mmmmmm", "m2_1_1mmmmm"),
    c("m1_1_mmmmmm", "m2_1_1mmmmm", "m3_1_11mmmm"),
    c("m1_1_mmmmmm", "m2_1_1mmmmm", "m3_1_11mmmm", "m4_1_111mmm"),
    c("m1_1_mmmmmm", "m2_1_1mmmmm", "m3_1_11mmmm", "m4_1_111mmm", "m5_1_1111mm"),
    c("m1_1_mmmmmm", "m2_1_1mmmmm", "m3_1_11mmmm", "m4_1_111mmm", "m5_1_1111mm", "m6_1_11111m")
  )
  expect_equal(joint_dist, expect_list)

  # Estimating the target quantities
  # Marginals under X=0
  # test for function marg_dist
  marg_dist = list()
  for (k in first:K) {
    marg_dist[[k]] = character(0)
    a <- 0
    marg_dist[[k]] = paste0("m", k, "_", a, "_", strrep("m", K))
  }
  expect_list <- list(
    NULL,
    NULL,
    "m3_0_mmmmmm",
    "m4_0_mmmmmm",
    "m5_0_mmmmmm",
    "m6_0_mmmmmm"
  )
  expect_equal(marg_dist, expect_list)


  # For p_first,..., p_K
  # Joint of others under X!=0
  # test for function joint_X_nonzero
  joint_X_nonzero = list()
  if (any(intervention_type %in% c("all", "shift_k")) && first<=(K-1)) {
    for (MM in first:(K-1)) {
      index = setdiff(MM:K, MM)
      joint_X_nonzero[[MM]] = list()
      for (k in index) {
        joint_X_nonzero[[MM]][[k]] = character(0)
        a = 1
        if (first != 1) {
          l = 1:(first-1)
          joint_X_nonzero[[MM]][[k]] = med_outcome_name(a = a, l = l, K = K)
        }
        if ((k-1) > first) {
          l = setdiff(first:(k-1), MM)
          joint_X_nonzero[[MM]][[k]] = c(joint_X_nonzero[[MM]][[k]], med_joint_other(k = l, a = a, MM = MM, K = K,
                                                        ordering = FALSE))
        }

        # Perform counterfactual prediction
        joint_X_nonzero[[MM]][[k]] = c(joint_X_nonzero[[MM]][[k]], med_joint_other(k = k, a = a, MM = MM, K = K, ordering = FALSE))
      }
    }
  }

  expected_list = list(NULL,NULL,
    list(NULL, NULL, NULL,
      c("m1_1_mmmmmm", "m2_1_1mmmmm", "m4_1_11mmmm"),
      c("m1_1_mmmmmm", "m2_1_1mmmmm", "m4_1_11mmmm", "m5_1_11m1mm"),
      c("m1_1_mmmmmm", "m2_1_1mmmmm", "m4_1_11mmmm", "m5_1_11m1mm", "m6_1_11m11m")
    ),
    list(NULL, NULL, NULL, NULL,
      c("m1_1_mmmmmm", "m2_1_1mmmmm", "m3_1_11mmmm", "m5_1_111mmm"),
      c("m1_1_mmmmmm", "m2_1_1mmmmm", "m3_1_11mmmm", "m5_1_111mmm", "m6_1_111m1m")
    ),
    list(NULL, NULL, NULL, NULL, NULL,
      c("m1_1_mmmmmm", "m2_1_1mmmmm", "m3_1_11mmmm", "m4_1_111mmm", "m6_1_1111mm")
    )
  )
  expect_equal(expected_list, joint_X_nonzero)

  # For p_first_prime,...., p_K_prime
  # Conditionals under X!=0
  con_exposed = list()
  if (any(intervention_type %in% c("all", "shift_k_order"))) {
    for (MM in first:(K - 1)) {
      con_exposed[[MM]] = list()
      for (k in (MM + 1):K) {
        a = 1
        con_exposed[[MM]][[k]] = character(0)
        if (MM != 1) {
          l = 1:(MM - 1)
          con_exposed[[MM]][[k]] = med_outcome_name(a = a,
                                                    l = l,
                                                    K = K)
        }

        # Handle mediator MM
        con_exposed[[MM]][[k]] = c(con_exposed[[MM]][[k]], paste0("m", MM, "_", 0, "_", strrep("m", K)))

        # Handle mediators between MM and k
        if (k > (MM + 1)) {
          l = (MM + 1):(k - 1)
          con_exposed[[MM]][[k]] = c(con_exposed[[MM]][[k]], med_joint_other(k = l, a = a, MM = MM, K = K))
        }

        con_exposed[[MM]][[k]] = c(con_exposed[[MM]][[k]], med_joint_other(k = k, a = a, MM = MM, K = K))
      }
    }
  }

  expected_list = list(NULL, NULL,
    list(NULL, NULL, NULL,
      c("m1_1_mmmmmm", "m2_1_1mmmmm", "m3_0_mmmmmm", "m4_1_110mmm"),
      c("m1_1_mmmmmm", "m2_1_1mmmmm", "m3_0_mmmmmm", "m4_1_110mmm", "m5_1_1101mm"),
      c("m1_1_mmmmmm", "m2_1_1mmmmm", "m3_0_mmmmmm", "m4_1_110mmm", "m5_1_1101mm", "m6_1_11011m")
    ),
    list(NULL, NULL, NULL, NULL,
      c("m1_1_mmmmmm", "m2_1_1mmmmm", "m3_1_11mmmm", "m4_0_mmmmmm", "m5_1_1110mm"),
      c("m1_1_mmmmmm", "m2_1_1mmmmm", "m3_1_11mmmm", "m4_0_mmmmmm", "m5_1_1110mm", "m6_1_11101m")
    ),
    list(NULL, NULL, NULL, NULL, NULL,
      c("m1_1_mmmmmm", "m2_1_1mmmmm", "m3_1_11mmmm", "m4_1_111mmm", "m5_0_mmmmmm", "m6_1_11110m")
    )
  )
  expect_equal(con_exposed, expected_list)

  # For p_all
  # Joint of main ones under X=0
  joint_unexposed = list()
  if (any(intervention_type %in% c("all", "shift_all"))) {
    for (k in (first + 1):K) {
      a <- 0
      joint_unexposed[[k]] = character(0)
      l = first:(k - 1)
      joint_unexposed[[k]] = med_outcome_all(l = l, first = first, a = a, K = K)

      joint_unexposed[[k]] = c(joint_unexposed[[k]], med_outcome_all(l = k, first = first, a = a, K = K))
    }
  }

  expect_list <- list(NULL, NULL, NULL,
    c("m3_0_mmmmmm", "m4_0_mm0mmm"),
    c("m3_0_mmmmmm", "m4_0_mm0mmm", "m5_0_mm00mm"),
    c("m3_0_mmmmmm", "m4_0_mm0mmm", "m5_0_mm00mm", "m6_0_mm000m")
  )
  expect_equal(joint_unexposed, expect_list)
  ###################
  # check for p_all #
  ###################
  p_all = vector()
  if (first > 1) {
    l = 1:(first - 1)
    p_all = med_outcome_name(a = 1, l = l, K = K)
  }

  # all mediators of interest
  k = first:K
  p_all = c(p_all, med_outcome_all(l = k, first = first, a = 0, K = K))
  expect_equal(c(joint_dist[[K]][1:(first-1)], joint_unexposed[[K]]), p_all)

  #################
  # check for p_k #
  #################

  p_k = list()
  a = 1
  for(MM in first:K){
    p_k[[MM]] = list()
    if (first > 1) {
      l = 1:(first - 1)
      p_k[[MM]] = med_outcome_name(a = a, l = l, K = K)
    }
    p_k[[MM]] = c(p_k[[MM]], paste0("m", MM, "_", 0, "_",
                                    strrep("m", K)))
    if (length(first:K) > 1){
      k = setdiff(first:K, MM)
      p_k[[MM]] = c(p_k[[MM]], med_joint_other(k = k, a = a, MM = MM, K = K,
                                                    ordering = FALSE))
    }
  }
  expect_true(all(unlist(c(joint_X_nonzero[[3]][[6]], marg_dist[3])) %in% p_k[[3]]))
  expect_true(all(unlist(c(joint_X_nonzero[[4]][[6]], marg_dist[4])) %in% p_k[[4]]))
  expect_true(all(unlist(c(joint_X_nonzero[[5]][[6]], marg_dist[5])) %in% p_k[[5]]))
  expect_true(all(unlist(c(joint_dist[[6]][1:5], marg_dist[6])) %in% p_k[[6]]))

  #######################
  # check for p_k_prime #
  #######################
  p_k_prime = list()
  a = 1
  for(MM in first:(K - 1)){
    p_k_prime[[MM]] = list()
    if (MM != 1) {
      l = 1:(MM - 1)
      p_k_prime[[MM]] = med_outcome_name(a = a, l = l, K = K)
    }
    p_k_prime[[MM]] = c(p_k_prime[[MM]], paste0("m", MM, "_", 0, "_",
                                         strrep("m", K)))
    if ((MM + 1) <= K) {
      k = (MM + 1):K
      p_k_prime[[MM]] = c(p_k_prime[[MM]], med_joint_other(k = k, a = a, MM = MM, K = K))
    }
  }
  expect_equal(p_k_prime[[3]], unlist(con_exposed[[3]][6]))
  expect_equal(p_k_prime[[4]], unlist(con_exposed[[4]][6]))
  expect_equal(p_k_prime[[5]], unlist(con_exposed[[5]][6]))
})





test_that("setting with first = 1, K = 4", {
  first = 1
  K = 4
  intervention_type = "all"
  # ESTIMATE DISTRIBUTIONS
  # Joint of M1 to MK under X=0 and X!=0 ...
  # testing function joint_dist

  joint_dist = list()
  for (k in 1:K) {
    joint_dist[[k]] = character(0)
    a = 1
    # pick up previous draws for prediction
    if (k != 1) {
      l <- 1:(k - 1)
      joint_dist[[k]] = med_outcome_name(a = a, l = l, K = K)
    }
    # random draws for the current
    joint_dist[[k]] = c(joint_dist[[k]], med_outcome_name(a = a, l = k, K = K))
  }
  expect_list = list(
    c("m1_1_mmmm"),
    c("m1_1_mmmm", "m2_1_1mmm"),
    c("m1_1_mmmm", "m2_1_1mmm", "m3_1_11mm"),
    c("m1_1_mmmm", "m2_1_1mmm", "m3_1_11mm", "m4_1_111m")
  )
  expect_equal(joint_dist, expect_list)

  # Estimating the target quantities
  # Marginals under X=0
  # test for function marg_dist
  marg_dist = list()
  for (k in first:K) {
    marg_dist[[k]] = character(0)
    a <- 0
    marg_dist[[k]] = paste0("m", k, "_", a, "_", strrep("m", K))
  }
  expect_list <- list(
    c("m1_0_mmmm"),
    c("m2_0_mmmm"),
    c("m3_0_mmmm"),
    c("m4_0_mmmm")
  )

  expect_equal(marg_dist, expect_list)


  # For p_first,..., p_K
  # Joint of others under X!=0
  # test for function joint_X_nonzero
  joint_X_nonzero = list()
  if (any(intervention_type %in% c("all", "shift_k")) && first<=(K-1)) {
    for (MM in first:(K-1)) {
      index = setdiff(MM:K, MM)
      joint_X_nonzero[[MM]] = list()
      for (k in index) {
        joint_X_nonzero[[MM]][[k]] = character(0)
        a = 1
        if (first != 1) {
          l = 1:(first-1)
          joint_X_nonzero[[MM]][[k]] = med_outcome_name(a = a, l = l, K = K)
        }
        if ((k-1) > first) {
          l = setdiff(first:(k-1), MM)
          joint_X_nonzero[[MM]][[k]] = c(joint_X_nonzero[[MM]][[k]], med_joint_other(k = l, a = a, MM = MM, K = K,
                                                                                     ordering = FALSE))
        }

        # Perform counterfactual prediction
        joint_X_nonzero[[MM]][[k]] = c(joint_X_nonzero[[MM]][[k]], med_joint_other(k = k, a = a, MM = MM, K = K, ordering = FALSE))
      }
    }
  }

  expected_list = list(list(NULL,
      c("m2_1_mmmm"),
      c("m2_1_mmmm", "m3_1_m1mm"),
      c("m2_1_mmmm", "m3_1_m1mm", "m4_1_m11m")),
    list(NULL, NULL,
      c("m1_1_mmmm", "m3_1_1mmm"),
      c("m1_1_mmmm", "m3_1_1mmm", "m4_1_1m1m")),
    list(NULL, NULL, NULL,
      c("m1_1_mmmm", "m2_1_1mmm", "m4_1_11mm"))
  )
  expect_equal(expected_list, joint_X_nonzero)

  # For p_first_prime,...., p_K_prime
  # Conditionals under X!=0
  con_exposed = list()
  if (any(intervention_type %in% c("all", "shift_k_order"))) {
    for (MM in first:(K - 1)) {
      con_exposed[[MM]] = list()
      for (k in (MM + 1):K) {
        a = 1
        con_exposed[[MM]][[k]] = character(0)
        if (MM != 1) {
          l = 1:(MM - 1)
          con_exposed[[MM]][[k]] = med_outcome_name(a = a,
                                                    l = l,
                                                    K = K)
        }

        # Handle mediator MM
        con_exposed[[MM]][[k]] = c(con_exposed[[MM]][[k]], paste0("m", MM, "_", 0, "_", strrep("m", K)))

        # Handle mediators between MM and k
        if (k > (MM + 1)) {
          l = (MM + 1):(k - 1)
          con_exposed[[MM]][[k]] = c(con_exposed[[MM]][[k]], med_joint_other(k = l, a = a, MM = MM, K = K))
        }

        con_exposed[[MM]][[k]] = c(con_exposed[[MM]][[k]], med_joint_other(k = k, a = a, MM = MM, K = K))
      }
    }
  }

  expected_list = list(list(NULL,
      c("m1_0_mmmm", "m2_1_0mmm"),
      c("m1_0_mmmm", "m2_1_0mmm", "m3_1_01mm"),
      c("m1_0_mmmm", "m2_1_0mmm", "m3_1_01mm", "m4_1_011m")),
    list(NULL, NULL,
      c("m1_1_mmmm", "m2_0_mmmm", "m3_1_10mm"),
      c("m1_1_mmmm", "m2_0_mmmm", "m3_1_10mm", "m4_1_101m")),
    list(NULL, NULL, NULL,
      c("m1_1_mmmm", "m2_1_1mmm", "m3_0_mmmm", "m4_1_110m")))
  expect_equal(con_exposed, expected_list)

  # For p_all
  # Joint of main ones under X=0
  joint_unexposed = list()
  if (any(intervention_type %in% c("all", "shift_all"))) {
    for (k in (first + 1):K) {
      a <- 0
      joint_unexposed[[k]] = character(0)
      l = first:(k - 1)
      joint_unexposed[[k]] = med_outcome_all(l = l, first = first, a = a, K = K)

      joint_unexposed[[k]] = c(joint_unexposed[[k]], med_outcome_all(l = k, first = first, a = a, K = K))
    }
  }

  expect_list <- list(NULL,
    c("m1_0_mmmm", "m2_0_0mmm"),
    c("m1_0_mmmm", "m2_0_0mmm", "m3_0_00mm"),
    c("m1_0_mmmm", "m2_0_0mmm", "m3_0_00mm", "m4_0_000m")
  )
  expect_equal(joint_unexposed, expect_list)
  ###################
  # check for p_all #
  ###################
  p_all = vector()
  if (first > 1) {
    l = 1:(first - 1)
    p_all = med_outcome_name(a = 1, l = l, K = K)
  }

  # all mediators of interest
  k = first:K
  p_all = c(p_all, med_outcome_all(l = k, first = first, a = 0, K = K))
  expect_equal(joint_unexposed[[K]], p_all)

  #################
  # check for p_k #
  #################

  p_k = list()
  a = 1
  for(MM in first:K){
    p_k[[MM]] = list()
    if (first > 1) {
      l = 1:(first - 1)
      p_k[[MM]] = med_outcome_name(a = a, l = l, K = K)
    }
    p_k[[MM]] = paste0("m", MM, "_", 0, "_",
                                    strrep("m", K))
    if (length(first:K) > 1){
      k = setdiff(first:K, MM)
      p_k[[MM]] = c(p_k[[MM]], med_joint_other(k = k, a = a, MM = MM, K = K,
                                               ordering = FALSE))
    }
  }
  expect_true(all(unlist(c(joint_X_nonzero[[1]][[4]], marg_dist[1])) %in% p_k[[1]]))
  expect_true(all(unlist(c(joint_X_nonzero[[2]][[4]], marg_dist[2])) %in% p_k[[2]]))
  expect_true(all(unlist(c(joint_X_nonzero[[3]][[4]], marg_dist[3])) %in% p_k[[3]]))
  expect_true(all(unlist(c(joint_dist[[4]][1:3], marg_dist[4])) %in% p_k[[4]]))

  #######################
  # check for p_k_prime #
  #######################
  p_k_prime = list()
  a = 1
  for(MM in first:(K - 1)){
    p_k_prime[[MM]] = list()
    if (MM != 1) {
      l = 1:(MM - 1)
      p_k_prime[[MM]] = med_outcome_name(a = a, l = l, K = K)
    }
    p_k_prime[[MM]] = c(p_k_prime[[MM]], paste0("m", MM, "_", 0, "_",
                                                strrep("m", K)))
    if ((MM + 1) <= K) {
      k = (MM + 1):K
      p_k_prime[[MM]] = c(p_k_prime[[MM]], med_joint_other(k = k, a = a, MM = MM, K = K))
    }
  }
  expect_equal(unlist(p_k_prime[[1]]), unlist(con_exposed[[1]][4]))
  expect_equal(p_k_prime[[2]], unlist(con_exposed[[2]][4]))
  expect_equal(p_k_prime[[3]], unlist(con_exposed[[3]][4]))
})


test_that("setting with first = 1, K = 2", {
  first = 1
  K = 2
  intervention_type = "all"
  # ESTIMATE DISTRIBUTIONS
  # Joint of M1 to MK under X=0 and X!=0 ...
  # testing function joint_dist

  joint_dist = list()
  for (k in 1:K) {
    joint_dist[[k]] = character(0)
    a = 1
    # pick up previous draws for prediction
    if (k != 1) {
      l <- 1:(k - 1)
      joint_dist[[k]] = med_outcome_name(a = a, l = l, K = K)
    }
    # random draws for the current
    joint_dist[[k]] = c(joint_dist[[k]], med_outcome_name(a = a, l = k, K = K))
  }
  expect_list = list(
    c("m1_1_mm"),
    c("m1_1_mm", "m2_1_1m"))
  expect_equal(joint_dist, expect_list)

  # Estimating the target quantities
  # Marginals under X=0
  # test for function marg_dist
  marg_dist = list()
  for (k in first:K) {
    marg_dist[[k]] = character(0)
    a <- 0
    marg_dist[[k]] = paste0("m", k, "_", a, "_", strrep("m", K))
  }
  expect_list <- list(
    c("m1_0_mm"),
    c("m2_0_mm")
  )

  expect_equal(marg_dist, expect_list)


  # For p_first,..., p_K
  # Joint of others under X!=0
  # test for function joint_X_nonzero
  joint_X_nonzero = list()
  if (any(intervention_type %in% c("all", "shift_k")) && first<=(K-1)) {
    for (MM in first:(K-1)) {
      index = setdiff(MM:K, MM)
      joint_X_nonzero[[MM]] = list()
      for (k in index) {
        joint_X_nonzero[[MM]][[k]] = character(0)
        a = 1
        if (first != 1) {
          l = 1:(first-1)
          joint_X_nonzero[[MM]][[k]] = med_outcome_name(a = a, l = l, K = K)
        }
        if ((k-1) > first) {
          l = setdiff(first:(k-1), MM)
          joint_X_nonzero[[MM]][[k]] = c(joint_X_nonzero[[MM]][[k]], med_joint_other(k = l, a = a, MM = MM, K = K,
                                                                                     ordering = FALSE))
        }

        # Perform counterfactual prediction
        joint_X_nonzero[[MM]][[k]] = c(joint_X_nonzero[[MM]][[k]], med_joint_other(k = k, a = a, MM = MM, K = K, ordering = FALSE))
      }
    }
  }

  expected_list = list(list(NULL,
                            c("m2_1_mm"))
  )
  expect_equal(expected_list, joint_X_nonzero)

  # For p_first_prime,...., p_K_prime
  # Conditionals under X!=0
  con_exposed = list()
  if (any(intervention_type %in% c("all", "shift_k_order"))) {
    for (MM in first:(K - 1)) {
      con_exposed[[MM]] = list()
      for (k in (MM + 1):K) {
        a = 1
        con_exposed[[MM]][[k]] = character(0)
        if (MM != 1) {
          l = 1:(MM - 1)
          con_exposed[[MM]][[k]] = med_outcome_name(a = a,
                                                    l = l,
                                                    K = K)
        }

        # Handle mediator MM
        con_exposed[[MM]][[k]] = c(con_exposed[[MM]][[k]], paste0("m", MM, "_", 0, "_", strrep("m", K)))

        # Handle mediators between MM and k
        if (k > (MM + 1)) {
          l = (MM + 1):(k - 1)
          con_exposed[[MM]][[k]] = c(con_exposed[[MM]][[k]], med_joint_other(k = l, a = a, MM = MM, K = K))
        }

        con_exposed[[MM]][[k]] = c(con_exposed[[MM]][[k]], med_joint_other(k = k, a = a, MM = MM, K = K))
      }
    }
  }

  expected_list = list(list(NULL,
                            c("m1_0_mm", "m2_1_0m")))
  expect_equal(con_exposed, expected_list)

  # For p_all
  # Joint of main ones under X=0
  joint_unexposed = list()
  if (any(intervention_type %in% c("all", "shift_all"))) {
    for (k in (first + 1):K) {
      a <- 0
      joint_unexposed[[k]] = character(0)
      l = first:(k - 1)
      joint_unexposed[[k]] = med_outcome_all(l = l, first = first, a = a, K = K)

      joint_unexposed[[k]] = c(joint_unexposed[[k]], med_outcome_all(l = k, first = first, a = a, K = K))
    }
  }

  expect_list <- list(NULL,
                      c("m1_0_mm", "m2_0_0m"))
  expect_equal(joint_unexposed, expect_list)
  ###################
  # check for p_all #
  ###################
  p_all = vector()
  if (first > 1) {
    l = 1:(first - 1)
    p_all = med_outcome_name(a = 1, l = l, K = K)
  }

  # all mediators of interest
  k = first:K
  p_all = c(p_all, med_outcome_all(l = k, first = first, a = 0, K = K))
  expect_equal(joint_unexposed[[K]], p_all)

  #################
  # check for p_k #
  #################

  p_k = list()
  a = 1
  for(MM in first:K){
    p_k[[MM]] = list()
    if (first > 1) {
      l = 1:(first - 1)
      p_k[[MM]] = med_outcome_name(a = a, l = l, K = K)
    }
    p_k[[MM]] = paste0("m", MM, "_", 0, "_",
                       strrep("m", K))
    if (length(first:K) > 1){
      k = setdiff(first:K, MM)
      p_k[[MM]] = c(p_k[[MM]], med_joint_other(k = k, a = a, MM = MM, K = K,
                                               ordering = FALSE))
    }
  }
  expect_true(all(unlist(c(joint_X_nonzero[[1]][[2]], marg_dist[1])) %in% p_k[[1]]))
  expect_true(all(unlist(c(joint_dist[[2]][1], marg_dist[2])) %in% p_k[[2]]))

  #######################
  # check for p_k_prime #
  #######################
  p_k_prime = list()
  a = 1
  for(MM in first:(K - 1)){
    p_k_prime[[MM]] = list()
    if (MM != 1) {
      l = 1:(MM - 1)
      p_k_prime[[MM]] = med_outcome_name(a = a, l = l, K = K)
    }
    p_k_prime[[MM]] = c(p_k_prime[[MM]], paste0("m", MM, "_", 0, "_",
                                                strrep("m", K)))
    if ((MM + 1) <= K) {
      k = (MM + 1):K
      p_k_prime[[MM]] = c(p_k_prime[[MM]], med_joint_other(k = k, a = a, MM = MM, K = K))
    }
  }
  expect_equal(unlist(p_k_prime[[1]]), unlist(con_exposed[[1]][2]))
})





