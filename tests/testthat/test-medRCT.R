test_that("test for medRCT", {
  # estimate TCE using g-comp
  fit = glm(
    child_mh ~
      (sep + fam_stress + parent_mh + preschool_att)^2 +
        sep * child_sex +
        sep * child_atsi +
        sep * mat_cob +
        sep *
          mat_engl +
        sep * mat_age,
    data = medRCT::LSACdata,
    family = binomial
  )
  data = LSACdata[complete.cases(LSACdata), ]
  data$sep = 0
  predA = predict(fit, data, type = "response")
  data$sep = 1
  predB = predict(fit, data, type = "response")
  tce = mean(predB) - mean(predA)

  set.seed(2024)
  result <- medRCT(
    dat = LSACdata,
    exposure = "sep",
    outcome = "child_mh",
    mediators = c("parent_mh", "preschool_att"),
    intermediate_confs = "fam_stress",
    confounders = c(
      "child_sex",
      "child_atsi",
      "mat_cob",
      "mat_engl",
      "mat_age"
    ),
    interactions_XC = "none",
    intervention_type = "all",
    mcsim = 50,
    bootstrap = F
  )

  expect_type(result$est, "double")
  dif = result$est["TCE (p_trt - p_ctr)"] - tce
  expect_lt(dif, 0.01)

  expect_message(
    result <- medRCT(
      dat = LSACdata,
      exposure = "sep",
      outcome = "child_mh",
      mediators = c("parent_mh"),
      intermediate_confs = "fam_stress",
      confounders = c(
        "child_sex",
        "child_atsi",
        "mat_cob",
        "mat_engl",
        "mat_age"
      ),
      interactions_XC = "none",
      intervention_type = "all",
      mcsim = 10,
      bootstrap = F
    ),
    "Only able to estimate the effect type 'shift_k' with a single mediator."
  )

  expect_error(
    result <- medRCT(
      dat = LSACdata,
      exposure = "sep",
      outcome = "child_mh",
      mediators = c("parent_mh"),
      intermediate_confs = "fam_stress",
      confounders = c(
        "child_sex",
        "child_atsi",
        "mat_cob",
        "mat_engl",
        "mat_age"
      ),
      interactions_XC = "none",
      intervention_type = "all",
      mcsim = 10,
      bootstrap = F,
      effect_measure = "Diff"
    ),
    "For binary outcome, effect_measure must be one of 'RD' or 'RR'."
  )

  result1 <- medRCT(
    dat = LSACdata,
    exposure = "sep",
    outcome = "child_SDQscore",
    mediators = c("parent_mh"),
    intermediate_confs = "fam_stress",
    confounders = c(
      "child_sex",
      "child_atsi",
      "mat_cob",
      "mat_engl",
      "mat_age"
    ),
    interactions_XC = "none",
    intervention_type = "shift_k",
    mcsim = 10,
    bootstrap = T,
    boot_args = list(R = 40, ci.type = "perc")
  )
  expect_equal(result1$effect_measure, "Diff")

  expect_error(
    result2 <- medRCT(
      dat = LSACdata,
      exposure = "sep",
      outcome = "child_SDQscore",
      mediators = c("parent_mh"),
      intermediate_confs = "fam_stress",
      confounders = c(
        "child_sex",
        "child_atsi",
        "mat_cob",
        "mat_engl",
        "mat_age"
      ),
      interactions_XC = "sep:child_sex",
      intervention_type = "all",
      mcsim = 10,
      bootstrap = F,
      effect_measure = "RR"
    ),
    "For continuous outcome, effect_measure must be 'Diff'."
  )

  LSACdata$sep1 = LSACdata$sep
  LSACdata$sep1[sample(which(LSACdata$sep == 0), 1500)] = 2
  result3 <- medRCT(
    dat = LSACdata,
    exposure = "sep1",
    outcome = "child_mh",
    mediators = c("parent_mh"),
    intermediate_confs = "fam_stress",
    confounders = c(
      "child_sex",
      "child_atsi",
      "mat_cob",
      "mat_engl",
      "mat_age"
    ),
    interactions_XC = "none",
    intervention_type = "all",
    mcsim = 10,
    bootstrap = F
  )
  expect_equal(result3$effect_measure, "RD")
  result4 <- medRCT(
    dat = LSACdata,
    exposure = "sep1",
    outcome = "child_SDQscore",
    mediators = c("parent_mh"),
    intermediate_confs = "fam_stress",
    confounders = c(
      "child_sex",
      "child_atsi",
      "mat_cob",
      "mat_engl",
      "mat_age"
    ),
    interactions_XC = "none",
    intervention_type = "all",
    mcsim = 10,
    bootstrap = F
  )
  expect_equal(result4$effect_measure, "Diff")
})
