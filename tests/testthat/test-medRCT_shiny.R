library(shinytest2)
library(shiny)
test_that("test medRCT_shiny", {
  expect_error(medRCT_shiny(), "A dataset must be provided to launch the app")
  expect_error(medRCT_shiny(data = "not_a_dataframe"), "The input data must be a data frame.")
})


test_that("Reactive inputs for medRCT_shiny", {
  app <- shinytest2::AppDriver$new(medRCT_shiny(data = LSACdata))

  app$set_inputs(outcome = "child_sex", exposure = "sep")
  # Verify input values
  expect_equal(app$get_value(input = "outcome"), "child_sex")
  expect_equal(app$get_value(input = "exposure"), "sep")
  app$stop()
})



test_that("test collect_models", {
  test_data <- data.frame(
    exposure = rbinom(100, 1, 0.5),
    outcome = rbinom(100, 1, 0.5),
    mediator1 = rbinom(100, 1, 0.3),
    mediator2 = rnorm(100),
    confounder1 = rnorm(100),
    confounder2 = rbinom(100, 1, 0.5)
  )
  result <- collect_models(
    data = test_data,
    exposure = "exposure",
    outcome = "outcome",
    mediators = c("mediator1", "mediator2"),
    intermediate_confs = c("confounder1"),
    confounders = c("confounder2"),
    interactions_XC = "all",
    intervention_type = "all"
  )
  expect_type(result, "list")
  expect_true(!is.null(result$`Outcome Model`))
  expect_true(length(result$`Outcome Model`) > 0)
  expect_true(!is.null(result$`Interventional Mediator Models`))
  expect_true(length(result$`Interventional Mediator Models`) > 0)


  expect_message(
    result <- collect_models(
      data = test_data,
      exposure = "exposure",
      outcome = "outcome",
      mediators = c("mediator1"),
      intermediate_confs = c("confounder1"),
      confounders = c("confounder2"),
      interactions_XC = "all",
      intervention_type = "all"
    ),
    "Only able to estimate the effect type 'shift_k' with a single mediator."
  )

})



test_that("test gen_formula_shiny", {
  formula <- gen_formula_shiny(
    mediators = c("mediator1"),
    exposure = "exposure",
    k = 1,
    interactions_XC = "confounder1"
  )
  expect_equal(formula, "mediator1 ~ exposure + confounder1")

  # marginals
  formula <- gen_formula_shiny(
    mediators = c("mediator1", "mediator2"),
    exposure = "exposure",
    k = 2,
    interactions_XC = "confounder1",
    marginal = TRUE
  )
  expect_equal(formula, "mediator2 ~ exposure + confounder1")
  # include all
  formula <- gen_formula_shiny(
    mediators = c("mediator1", "mediator2", "mediator3"),
    exposure = "exposure",
    k = 3,
    interactions_XC = "confounder1",
    include_all = TRUE
  )
  expect_equal(formula, "mediator3 ~ (exposure+mediator1+mediator2)^2 +confounder1")

  formula <- gen_formula_shiny(
    mediators = c("mediator1", "mediator2", "mediator3"),
    exposure = "exposure",
    k = 3,
    MM = 2,
    K = 3,
    interactions_XC = "confounder1"
  )
  expect_equal(formula, "mediator3 ~ (exposure + mediator1)^2 +confounder1")

  formula <- gen_formula_shiny(
    mediators = c("mediator1", "mediator2", "mediator3"),
    exposure = "exposure",
    k = 2,
    first = 1,
    MM = 1,
    K = 3,
    interactions_XC = "confounder1"
  )
  expect_equal(formula, "mediator2 ~ exposure + confounder1")

  formula <- gen_formula_shiny(
    mediators = c("mediator1", "mediator2", "mediator3", "mediator4"),
    exposure = "exposure",
    k = 4,
    MM = 2,
    first = 1,
    K = 4,
    interactions_XC = "confounder1+confounder2",
    include_all = TRUE
  )
  expect_equal(formula, "mediator4 ~ (exposure+mediator1+mediator2+mediator3)^2+confounder1+confounder2")
})


test_that("testing clean_list removes NULL elements", {
  input <- list(a = list(1, NULL, 3), b = list(), c = list(4, 5), d=list())
  expected <- list(a = list(1, 3), c = list(4, 5))
  expect_equal(clean_list(input), expected)
  x <- list(
    list(a = 1, b = NULL, c = 3),
    list(d = NULL, e = 5),
    "non-list element",
    list()
  )

  # Expected output
  expected <- list(
    list(a = 1, c = 3), # NULL removed
    list(e = 5),        # NULL removed
    "non-list element", # Non-list elements unchanged
    list()              # Empty list unchanged
  )

  # Apply the specific part of the function
  result <- lapply(x, function(y) {
    if (is.list(y)) Filter(Negate(is.null), y) else y
  })

  # Compare result to expected output
  expect_equal(result, expected)
})


test_that("catch_model_messages captures warnings", {
  data <- data.frame(
    x = c(1, 2, 3, 4, 5),
    y = c(0, 0, 1, 1, 1)
  )
  formula <- y ~ x
  family <- binomial()

  model <- catch_model_messages(formula, data, family)

  expect_s3_class(model, "glm")
  expect_true(length(model$warnings) > 0)
  expect_null(model$errors)

})




test_that("Server logic of medRCT_shiny works as expected", {
  # Test data
  test_data <- data.frame(
    treatment = rbinom(100, 1, 0.5),
    response = rnorm(100),
    mediator1 = rnorm(100),
    mediator2 = rnorm(100),
    conf1 = rnorm(100),
    conf2 = rbinom(100, 1, 0.5)
  )

  # Test server
  testServer(medRCT_shiny(test_data),{
    # Set inputs
    session$setInputs(
      exposure = "treatment",
      outcome = "response",
      mediators = c("mediator1", "mediator2"),
      confounders = c("conf2"),
      int_confs = "conf1",
      interactions_XC = "all",
      intervention_type = "all",
      med_button = 1
    )

    # Test reactive `var_names`
    expect_equal(var_names(), c("treatment", "response"))
    expect_true("mediator1" %in% input$mediators)
    expect_true("conf1" %in% input$int_confs)
    expect_true(!is.null(model.list$model))
    expect_named(model.list$model$`Joint Mediator Models`,
                 c("Model_1", "Model_2", "Model_3"))
    expect_named(model.list$model,
                 c('Joint Mediator Models', 'Interventional Mediator Models',
                 'Outcome Model'))
  })
})
