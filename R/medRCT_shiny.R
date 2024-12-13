

#' Launch Interactive Shiny App for Model Diagnostics
#'
#' `medRCT_shiny` launches a Shiny application to provide model summaries
#' of all models fitted by the algorithm. The app serves as a user-friendly
#' interface for diagnosing potential issues with the models and facilitates adjustments
#' to improve model fit.
#'
#' @param data A \code{data.frame} containing the dataset for analysis. It should include variables for the exposure,
#'  outcome, mediators, confounders, and exposure-induced mediator-outcome confounders specified in the analysis.
#' @param ... additional arguments for shiny
#'
#' @importFrom stats formula
#' @import shiny
#'
#' @export
#'
#' @examples
#' if (interactive()) {
#'    medRCT_shiny(data=LSACdata)
#' }
medRCT_shiny <- function(data, ...){
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Please install shiny: install.packages('shiny')",
         call. = FALSE)}

  if (missing(data)) stop("A dataset must be provided to launch the app.")
  # Ensure the dataset is a data.frame
  if (!is.data.frame(data)) stop("The input data must be a data frame.")

  # UI
  ui <- fluidPage(
    tags$head(
      tags$style(HTML("
      .warn pre {
        color: #721C24;
        background-color: #FFE5B4;
        font-weight: bolder;
      }"))
    ),
    # App title
    titlePanel("medRCT: Model Diagnostics"),

    # Sidebar layout
    sidebarLayout(
      # Sidebar panel for inputs
      sidebarPanel(
        #select outcome
        selectInput("outcome", "Select the Outcome Variable:",
                    choices = c("", colnames(data))),
        #select exposure
        selectInput("exposure", "Select the Exposure Variable:",
                    choices = c("", colnames(data))),
        #select Confounders
        selectInput("confounders", "Select the Baseline Confounders:",
                    choices = c("", colnames(data)),
                    multiple = TRUE),
        #select intermediate Confounders
        selectInput("int_confs", "Select the Intermediate Confounders:",
                    choices = c("NULL", colnames(data)),
                    multiple = TRUE),
        #select mediators
        selectInput("mediators", "Select the Mediators:",
                    choices = c("", colnames(data)),
                    multiple = TRUE),
        selectizeInput(
          inputId = "interactions_XC",
          label = "Specifying the exposure-confounder or confounder-confounder interaction:",
          choices = c("all", "none"), # Predefined options
          multiple = FALSE,
          options = list(create = TRUE) # Enable free text input
        ),
        selectInput(
          inputId = "intervention_type",
          label = "Select Intervention Type:",
          choices = c("all", "shift_all", "shift_k", "shift_k_order"),
          selected = "all" # Default selected value
        ),
        actionButton("med_button", "Run medRCT models",
                     style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"
        )
      ),

      # Main panel for displaying results
      mainPanel(
        uiOutput("model_panels")
      )
    )
  )

  # the server
  server <- function(input, output, session) {

    var_names <- reactive({
      c(input$exposure, input$outcome)
      })
    # update choices when the variables change
    observe({
      var_names()
      isolate({
        updateSelectInput(
          session, "outcome",
          choices = setdiff(names(data), setdiff(c(var_names(),
                                                   input$confounders,
                                                   input$int_confs,
                                                   input$mediators),
                                                 input$outcome)),
          selected = input$outcome
        )
        updateSelectInput(
          session, "confounders",
          choices = setdiff(names(data), setdiff(c(var_names(),
                                                   input$int_confs,
                                                   input$mediators),
                                                 input$confounders)),
          selected = input$confounders
        )
        updateSelectInput(
          session, "exposure",
          choices = setdiff(names(data), setdiff(c(var_names(),
                                                   input$confounders,
                                                   input$int_confs,
                                                   input$mediators),
                                                 input$exposure)),
          selected = input$exposure
        )
        updateSelectInput(
          session, "int_confs",
          choices = setdiff(c("NULL" ,names(data)), setdiff(c(var_names(),
                                                              input$confounders,
                                                              input$mediators),
                                                            input$int_confs)),
          selected = input$int_confs
        )
        updateSelectInput(
          session, "mediators",
          choices = setdiff(names(data), setdiff(c(var_names(),
                                                   input$confounders,
                                                   input$int_confs),
                                                 input$mediators)),
          selected = input$mediators
        )
      })
    })

    observe({
      input$confounders
      isolate({
        updateSelectInput(
          session, "exposure",
          choices = setdiff(names(data), setdiff(c(var_names(),
                                                   input$confounders,
                                                   input$int_confs,
                                                   input$mediators),
                                                 input$exposure)),
          selected = input$exposure
        )
        updateSelectInput(
          session, "outcome",
          choices = setdiff(names(data), setdiff(c(var_names(),
                                                   input$confounders,
                                                   input$int_confs,
                                                   input$mediators),
                                                 input$outcome)),
          selected = input$outcome
        )
        updateSelectInput(
          session, "int_confs",
          choices = setdiff(c("NULL", names(data)),
                            setdiff(c(var_names(),
                                      input$confounders,
                                      input$mediators),
                                    input$int_confs)),
          selected = input$int_confs
        )
        updateSelectInput(
          session, "mediators",
          choices = setdiff(names(data), setdiff(c(var_names(),
                                                   input$confounders,
                                                   input$int_confs),
                                                 input$mediators)),
          selected = input$mediators
        )
      })
    })

    observe({
      input$int_confs
      isolate({
        updateSelectInput(
          session, "exposure",
          choices = setdiff(names(data), setdiff(c(var_names(),
                                                   input$int_confs,
                                                   input$confounders,
                                                   input$mediators),
                                                 input$exposure)),
          selected = input$exposure
        )
        updateSelectInput(
          session, "outcome",
          choices = setdiff(names(data), setdiff(c(var_names(),
                                                   input$int_confs,
                                                   input$confounders,
                                                   input$mediators),
                                                 input$outcome)),
          selected = input$outcome
        )
        updateSelectInput(
          session, "confounders",
          choices = setdiff(names(data), setdiff(c(var_names(),
                                                   input$int_confs,
                                                   input$mediators),
                                                 input$confounders)),
          selected = input$confounders
        )
        updateSelectInput(
          session, "mediators",
          choices = setdiff(names(data), setdiff(c(var_names(),
                                                   input$confounders,
                                                   input$int_confs),
                                                 input$mediators)),
          selected = input$mediators
        )
      })
    })

    observe({
      input$mediators
      isolate({
        updateSelectInput(
          session, "exposure",
          choices = setdiff(names(data), setdiff(c(var_names(),
                                                   input$int_confs,
                                                   input$confounders,
                                                   input$mediators),
                                                 input$exposure)),
          selected = input$exposure
        )
        updateSelectInput(
          session, "outcome",
          choices = setdiff(names(data), setdiff(c(var_names(),
                                                   input$int_confs,
                                                   input$confounders,
                                                   input$mediators),
                                                 input$outcome)),
          selected = input$outcome
        )
        updateSelectInput(
          session, "confounders",
          choices = setdiff(names(data), setdiff(c(var_names(),
                                                   input$int_confs,
                                                   input$mediators),
                                                 input$confounders)),
          selected = input$confounders
        )
        updateSelectInput(
          session, "int_confs",
          choices = setdiff(c("NULL", names(data)),
                            setdiff(c(var_names(),
                                      input$confounders,
                                      input$mediators),
                                    input$int_confs)),
          selected = input$int_confs
        )
      })
    })

    model.list <- reactiveValues()


    # generate the models when clicking the button
    observe({
      req(input$med_button)
      isolate({
        req(input$exposure,
            input$outcome,
            input$mediators,
            input$int_confs,
            input$confounders)
        int_confs <- if (is.null(input$int_confs) || all(input$int_confs == "NULL")) {
          NULL
        } else {
          input$int_confs
        }
        model.list$model = collect_models(data= data,
                                          exposure = input$exposure,
                                          outcome = input$outcome,
                                          mediators = input$mediators,
                                          intermediate_confs = int_confs,
                                          confounders = input$confounders,
                                          interactions_XC = input$interactions_XC,
                                          intervention_type = input$intervention_type)
        # remove NULL and list() from the list
        model.list$model = clean_list(model.list$model)
        # rename the list
        for (i in seq_along(model.list$model)) {
          names(model.list$model[[i]]) <- paste0("Model_", seq_along(model.list$model[[i]]))
        }
      })
    })


    # main panel
    output$model_panels <- renderUI({
      req(model.list$model)
      req(input$med_button)
      isolate({
        tagList(
          h3("Models needed to estimate the interventional effect:"),
          # Create tabs dynamically based on first-level keys in the list
          do.call(tabsetPanel, lapply(names(model.list$model), function(group_name) {
            tabPanel(
              title = group_name,
              # Add dynamic sub-tabs for each model in the group
              do.call(tabsetPanel, lapply(names(model.list$model[[group_name]]), function(model_name) {
                model <- model.list$model[[group_name]][[model_name]]

                # Build the tabPanel inline
                tabPanel(
                  title = model_name,
                  h4("Model Formula"),
                  verbatimTextOutput(outputId = paste0(group_name, "_", model_name, "_formula")),
                  # Conditionally add warnings section
                  if (!is.null(model[['warnings']]) && length(model[['warnings']]) > 0) {
                    tagList(
                      h4("Warnings"),
                      div(class = "warn",
                          verbatimTextOutput(outputId = paste0(group_name, "_", model_name, "_warnings")))
                    )
                  },
                  # Conditionally add errors section
                  if (!is.null(model[['errors']]) && length(model[['errors']]) > 0) {
                    tagList(
                      h4("Errors"),
                      div(class = "warn",
                          verbatimTextOutput(outputId = paste0(group_name, "_", model_name, "_errors")))
                    )
                  },
                  h4("Model Summary"),
                  verbatimTextOutput(outputId = paste0(group_name, "_", model_name, "_summary")),
                  tags$div(style = "height: 100px;")
                )
              }))
            )
          }))
        )
      })
    })




    observe({
      req(model.list$model)
      lapply(names(model.list$model), function(group_name) {
        lapply(names(model.list$model[[group_name]]), function(model_name) {
          model <- model.list$model[[group_name]][[model_name]]

          # Render formula
          output[[paste0(group_name, "_", model_name, "_formula")]] <- renderPrint({
            cat(paste(deparse(stats::formula(model)), collapse = "\n"))
          })

          # Render summary
          output[[paste0(group_name, "_", model_name, "_summary")]] <- renderPrint({
            summary(model)
          })
          if (!is.null(model[['warnings']]) && length(model[['warnings']]) > 0) {
            output[[paste0(group_name, "_", model_name, "_warnings")]] <- renderPrint({
              cat(paste(model[['warnings']], collapse = "\n"))
            })
          }

          # Render errors if they exist
          if (!is.null(model[['errors']]) && length(model[['errors']]) > 0) {
            output[[paste0(group_name, "_", model_name, "_errors")]] <- renderPrint({
              cat(paste(model[['errors']], collapse = "\n"))
            })
          }
        })
      })
    })

  }
  # Run the application
  shinyApp(ui = ui, server = server)
}




#' Collect All Models fitted by the algorithm
#'
#' `collect_models` fits and collects models required for the algorithm.
#'
#' @param data A \code{data.frame} containing the dataset for analysis.
#' @param exposure A \code{character} string specifying the name of the exposure variable in the dataset.
#'  The exposure variable can be categorical, with \code{0} explicitly denoting the unexposed (or control) group.
#'  Other values represent different exposure categories.
#' @param outcome A \code{character} string specifying the name of the outcome variable in the dataset.
#'  The outcome variable can be specified as either binary or continuous.
#' @param mediators A \code{character} vector specifying the names of mediator variables in the dataset. The mediators
#'  can be specified as either binary or continuous. When estimating the effect type \code{"shift_k_order"},
#'  the order of mediators in the vector is important, as it determines the causal sequence of mediators.
#' @param intermediate_confs A \code{character} vector specifying the names of intermediate confounders in the dataset.
#'  The intermediate confounders can be specified as either binary or continuous. If \code{NULL},
#'  no intermediate confounders are specified, and the natural effect will be estimated.
#' @param confounders A \code{character} vector listing the names of baseline confounders.
#' @param interactions_XC A \code{character} string specifying the exposure-confounder or confounder-confounder
#'  interaction terms to include in the regression models for confounder adjustment. The default value, \code{"all"},
#'  includes all two-way exposure-confounder interactions but excludes confounder-confounder interactions.
#'  Specify \code{"none"} to exclude all two-way exposure-confounder and confounder-confounder interactions.
#' @param intervention_type A \code{character} string indicating the type of interventional effect to be estimated.
#'
#' @importFrom stats as.formula
#'
#' @keywords internal
collect_models <- function(data,
                           exposure,
                           outcome,
                           mediators,
                           intermediate_confs,
                           confounders,
                           interactions_XC = "all",
                           intervention_type = c("all", "shift_all", "shift_k", "shift_k_order")) {
  intervention_type = sapply(intervention_type, function(arg)
    match.arg(
      arg,
      choices = c("all", "shift_all", "shift_k", "shift_k_order")
    ))

  # set intervention type to shift_k when K==1
  if (length(mediators) == 1 &
      any(intervention_type %in% c("all", "shift_all", "shift_k_order"))) {
    intervention_type = "shift_k"
    message("Only able to estimate the effect type 'shift_k' with a single mediator.")
  }

  mediators = c(intermediate_confs, mediators)

  # define the first mediator of interest
  first = length(intermediate_confs) + 1

  K <- length(mediators)

  data <- as.data.frame(data)

  no.miss = nrow(data) - sum(stats::complete.cases(data))
  message(paste0(
    "Conducting complete case analysis, ",
    no.miss,
    " observations were excluded due to missing data.\n"
  ))
  data <- data[stats::complete.cases(data), ]

  fam_type = family_type(data, mediators)

  if (interactions_XC == "all") {
    interactions_XC <- paste(paste(rep(exposure, length(confounders)), confounders, sep =
                                     "*"), collapse = "+")
  } else if (interactions_XC == "none") {
    interactions_XC <- paste(confounders, collapse = "+")
  } else {
    interactions_XC <- interactions_XC
  }

  data[[paste0(exposure)]] = as.factor(data[[paste0(exposure)]])
  exposure_level = sort(unique(as.numeric(data[[paste0(exposure)]])))
  lnzero = exposure_level[exposure_level!=0]
  res = list()
  res[1:6] <- list(list())
  names(res) = c("all mediator effects", "individual mediator effect", "shift_k effects", "shift_k_order effects", "shift_all effects", "outcome regression")

  # Joint of M1 to MK under X=0 and X!=0 ...
  for (k in 1:K) {
    res[[names(res)[1]]][[k]] <- catch_model_messages(as.formula(
      gen_formula_shiny(mediators = mediators,
                        exposure = exposure,
                        k = k,
                        interactions_XC = interactions_XC,
                        include_all = TRUE)),
      data = data,
      family = fam_type[[k]])
    }
  # Marginals under X=0
  for (k in first:K) {
    res[[names(res)[2]]][[k]] <- catch_model_messages(as.formula(
      gen_formula_shiny(mediators = mediators,
                        exposure = exposure,
                        k = k,
                        interactions_XC = interactions_XC,
                        marginal = TRUE)),
      data = data,
      family = fam_type[[k]])
  }
  # Joint of others under X!=0
  if (any(intervention_type %in% c("all", "shift_k"))) {
    for (MM in first:K) {
      index = setdiff(first:K, MM)
      for (k in index) {
        if (!(first == 1 && MM != 1 && k == index[1])) {
          res[[names(res)[3]]][[paste(MM, k, sep = "_")]] <- catch_model_messages(as.formula(
            gen_formula_shiny(mediators = mediators,
                              exposure = exposure,
                              k = k,
                              MM = MM,
                              first = first,
                              K = K,
                              interactions_XC = interactions_XC)),
            data = data,
            family = fam_type[[k]])
        }
      }
    }
  }

  # For p_first_prime,...., p_K_prime
  # Conditionals under X!=0
  if (any(intervention_type %in% c("all", "shift_k_order"))) {
    for (MM in first:(K - 1)) {
      for (k in (MM + 1):K) {
        res[[names(res)[4]]][[paste(MM, k, sep = "_")]] <- catch_model_messages(as.formula(
          gen_formula_shiny(mediators = mediators,
                            exposure = exposure,
                            k = k,
                            interactions_XC = interactions_XC,
                            include_all = TRUE)),
          data = data,
          family = fam_type[[k]])
      }
    }
  }
  # For p_all
  # Joint of main ones under X=0
  if (any(intervention_type %in% c("all", "shift_all"))) {
    for (k in (first + 1):K) {
      res[[names(res)[5]]][[k]] <- catch_model_messages(as.formula(
        gen_formula_shiny(mediators = mediators,
                          exposure = exposure,
                          k = k,
                          interactions_XC = interactions_XC,
                          first = first,
                          include_all = TRUE)),
        data = data,
        family = fam_type[[k]])
    }
  }

  # outcome
  outcome_type = family_type(data, outcome)
  res[[names(res)[6]]][[1]] <- catch_model_messages(as.formula(paste0(
    outcome, " ~ (", exposure, " + ", paste0(mediators[1:K], collapse = "+"), ")^2+",
    interactions_XC)),
    data = data,
    family = outcome_type[[1]])

  return(res)
}




#' Generate Model Formulas for Shiny-Based Model Diagnostics
#'
#' This function generates model formulas using the original variable names as they appear in the dataset.
#'
#' @param mediators A \code{character} vector specifying the names of mediator variables in the dataset. The mediators
#'  can be specified as either binary or continuous. When estimating the effect type \code{"shift_k_order"},
#'  the order of mediators in the vector is important, as it determines the causal sequence of mediators.
#' @param exposure A \code{character} string specifying the name of the exposure variable in the dataset.
#'  The exposure variable can be categorical, with \code{0} explicitly denoting the unexposed (or control) group.
#'  Other values represent different exposure categories.
#' @param k An integer specifying the index of the mediator for which the formula is being generated.
#' @param first An integer (optional) specifying the index of the first mediator. Defaults to `NULL`,
#'  indicating that the function will generate formulas without explicitly considering this parameter.
#' @param MM An integer (optional) specifying the index of the mediator whose distribution will be shifted.
#'  Defaults to `NULL`, indicating that the function will generate formulas without explicitly considering
#'  this parameter.
#' @param K An integer (optional) specifying the total number of mediators and intermediate confounders.
#'  Defaults to `NULL`, indicating that the function will generate formulas without explicitly considering this
#'  parameter.
#' @param interactions_XC A \code{character} string specifying the exposure-confounder or confounder-confounder
#'  interaction terms to include in the regression models for confounder adjustment. The default value, \code{"all"},
#'  includes all two-way exposure-confounder interactions but excludes confounder-confounder interactions.
#'  Specify \code{"none"} to exclude all two-way exposure-confounder and confounder-confounder interactions.
#' @param include_all Logical.
#' @param marginal Logical. If `TRUE`, estimating marginals under `X=0`.
#'
#' @keywords internal
gen_formula_shiny <- function(
    mediators = mediators, exposure = exposure,
    k, first=NULL, MM = NULL, K=NULL, interactions_XC,
    include_all = FALSE, marginal = FALSE) {
  if ((k == 1 || marginal) && is.null(MM))  {
    return(paste0(mediators[k], " ~ ", exposure, " + ", interactions_XC))
  } else if (include_all) {
    if (is.null(first)){
      return(paste0(
        mediators[k], " ~ (", exposure, "+",
        paste0(mediators[1:(k - 1)], collapse = "+"),
        ")^2 +", interactions_XC
      ))
    } else {
      paste0(mediators[k], " ~ (", exposure, "+",
             paste0(mediators[first:(k - 1)], collapse = "+"),
             ")^2+",
             interactions_XC)
    }
  } else if (!is.null(MM)) {
    if (first == 1 && MM == 1 && k == setdiff(first:K, MM)[1]) {
      return(paste0(mediators[k], " ~ ", exposure, " + ", interactions_XC))
    } else {
      return(paste0(
        mediators[k], " ~ (", exposure, " + ",
        paste0(mediators[setdiff(1:(k - 1), MM)], collapse = "+"),
        ")^2 +", interactions_XC
      ))
    }
  }
}


#' Clean Nested Lists
#'
#' `clean_list` removes empty lists and `NULL` elements from a nested list structure,
#' cleaning up the list for further processing.
#'
#' @param x A nested list.
#'
#' @details
#' This function performs two cleaning operations on a list:
#' 1. At the first level, it removes any empty lists (`list()`).
#' 2. At the second level, it removes `NULL` elements from sublists while preserving
#'    other elements.
#'
#'
#' @return
#' A cleaned list with empty lists and `NULL` elements removed.
#'
#' @keywords internal
clean_list <- function(x) {
  # Remove empty lists (list()) at the first level
  x <- Filter(function(y) !(is.list(y) && length(y) == 0), x)

  # Remove NULL elements in sublists at the second level
  x <- lapply(x, function(y) {
    if (is.list(y)) Filter(Negate(is.null), y) else y
  })

  return(x)
}


#' Capture Warnings and Errors from Model Fitting
#'
#' `catch_model_messages` fits a generalised linear model while capturing
#' and storing any warnings or errors generated during the fitting process.
#' This is useful for debugging.
#'
#' @param formula A formula specifying the model to be fitted.
#' @param data A data frame containing the variables referenced in the formula.
#' @param family A description of the error distribution and link function to
#'   be used in the model, as specified in `glm` (it can be either `binomial` or `gaussian`).
#'
#' @details
#' This function uses `tryCatch` and `withCallingHandlers` to handle warnings
#' and errors separately. Warnings are captured and stored in the resulting
#' model object under the `warnings` attribute, while errors are stored under
#' the `errors` attribute. If the model fits successfully without warnings or errors,
#' these attributes will be `NULL`.
#'
#' @keywords internal
catch_model_messages <- function(formula, data, family) {

  # Capture warnings and errors
  tryCatch(
    {
      warnings <- NULL
      withCallingHandlers(
        {
          # fit the model
          model <- glm(formula, data = data, family = family)
        },
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
      model$warnings <- warnings
    },
    error = function(e) {
      model$errors <- conditionMessage(e)
    }
  )

  return(model)
}
