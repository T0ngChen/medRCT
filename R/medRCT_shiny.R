

#' Launch Interactive Shiny App for Model Assessment
#'
#' `medRCT_shiny` launches a Shiny application to provide model summaries
#' of all models fitted by the algorithm. The app provides a user-friendly
#' interface for model assessment and facilitates adjustments of interaction terms
#' to improve model fit.
#'
#' @param data A \code{data.frame} containing the dataset for analysis. It should include variables for the exposure,
#'  outcome, mediators, confounders, and exposure-induced mediator-outcome confounders specified in the analysis.
#' @param ... additional arguments for shiny
#'
#' @importFrom stats formula
#' @importFrom stats setNames
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
    titlePanel("medRCT: Model Assessment"),

    # Sidebar layout
    sidebarLayout(
      # Sidebar panel for inputs
      sidebarPanel(
        helpText(HTML("For detailed guidance on specifying these arguments,
                      refer to the <a href='https://t0ngchen.github.io/medRCT/articles/intro.html' target='_blank'>medRCT vignette</a>.")),
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
        helpText(HTML("Note: For <code>interaction terms</code>, you can choose <code>all</code> or <code>none</code>.
                      If you want to specify custom interaction terms, <b>you must type them manually</b>.")),
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
            input$intervention_type,
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
        model.list$newmod = list(
          "Mediator Models" = list(),
          "Outcome Model" = model.list$model[["Outcome Model"]]
        )
        all_mediator_models <- c(
          lapply(model.list$model$L_models, function(m) list(model = m, type = "L_models")),
          lapply(model.list$model$M_models, function(m) list(model = m, type = "M_models"))
        )
        all_labels <- unique(unlist(lapply(all_mediator_models, function(m) m$model$labels)))

        model.list$newmod[["Mediator Models"]] <- setNames(
          lapply(all_labels, function(label) {
            list(
              "L_models" = lapply(
                Filter(function(m) label %in% m$model$labels && m$type == "L_models", all_mediator_models),
                function(m) m$model
              ),
              "M_models" = lapply(
                Filter(function(m) label %in% m$model$labels && m$type == "M_models", all_mediator_models),
                function(m) m$model
              )
            )
          }),
          all_labels
        )
        model.list$newmod <- flatten_models(model.list$newmod)
        model.list$newmod$`Outcome Model`[["Outcome_model"]] = model.list$newmod$`Outcome Model`[[1]]
        model.list$newmod$`Outcome Model`[[1]] = NULL
      })
    })


    # main panel
    # Main Panel
    output$model_panels <- renderUI({
      req(model.list$newmod)
      req(input$med_button)

      isolate({
        tagList(
          h3(paste0("Models needed to estimate the expected outcome in relevant target trial arms for the specified interventional effect(s) of interest: ",
                    ifelse(input$intervention_type == "all", "all interventional effects", input$intervention_type))),

          # Create first-level tabs dynamically
          do.call(tabsetPanel, lapply(names(model.list$newmod), function(group_name) {

            # Condition to check if it's "Outcome Model" (skip arm_name level)
            if (group_name == "Outcome Model") {
              tabPanel(
                title = group_name,
                do.call(tabsetPanel, lapply(names(model.list$newmod[[group_name]]), function(model_name) {

                  model <- model.list$newmod[[group_name]][[model_name]]

                  tabPanel(
                    title = model_name,
                    h4("Model Formula"),
                    verbatimTextOutput(outputId = paste0(group_name, "_", model_name, "_formula")),

                    # Warnings Section (only if available)
                    if (!is.null(model[['warnings']]) && length(model[['warnings']]) > 0) {
                      tagList(
                        h4("Warnings"),
                        div(class = "warn", verbatimTextOutput(outputId = paste0(group_name, "_", model_name, "_warnings")))
                      )
                    },

                    # Errors Section (only if available)
                    if (!is.null(model[['errors']]) && length(model[['errors']]) > 0) {
                      tagList(
                        h4("Errors"),
                        div(class = "warn", verbatimTextOutput(outputId = paste0(group_name, "_", model_name, "_errors")))
                      )
                    },

                    h4("Model Summary"),
                    verbatimTextOutput(outputId = paste0(group_name, "_", model_name, "_summary")),
                    tags$div(style = "height: 100px;")
                  )
                }))
              )

            } else {  # Keep Treatment Arm nesting for "Mediator Models"
              tabPanel(
                title = group_name,
                do.call(tabsetPanel, lapply(
                  names(model.list$newmod[[group_name]])[order(
                                              ifelse(names(model.list$newmod[[group_name]]) == "p_trt / p_ctr", 1,
                                                     ifelse(names(model.list$newmod[[group_name]]) == "p_all", 2,
                                                            ifelse(grepl("prime", names(model.list$newmod[[group_name]])), 4, 3))),
                                              # Extract the numeric part when it exists (using a regex that grabs digits)
                                              ifelse(grepl("[0-9]+", names(model.list$newmod[[group_name]])), as.numeric(gsub("[^0-9]", "", names(model.list$newmod[[group_name]]))), 0)
                                            )], function(arm_name) {
                  tabPanel(
                    h4("L-models correspond to intermediate confounder models (e.g., in L1_model, the response variable is the
                       first intermediate confounder), and M-models correspond to mediator models (e.g., in M1_model,
                       the response variable is the first mediator of interest)"),
                    title = arm_name,
                    do.call(tabsetPanel, lapply(names(model.list$newmod[[group_name]][[arm_name]]), function(model_name) {

                      model <- model.list$newmod[[group_name]][[arm_name]][[model_name]]

                      tabPanel(
                        title = model_name,
                        h4("Model Formula"),
                        verbatimTextOutput(outputId = paste0(group_name, "_", arm_name, "_", model_name, "_formula")),

                        # Warnings Section (only if available)
                        if (!is.null(model[['warnings']]) && length(model[['warnings']]) > 0) {
                          tagList(
                            h4("Warnings"),
                            div(class = "warn", verbatimTextOutput(outputId = paste0(group_name, "_", arm_name, "_", model_name, "_warnings")))
                          )
                        },

                        # Errors Section (only if available)
                        if (!is.null(model[['errors']]) && length(model[['errors']]) > 0) {
                          tagList(
                            h4("Errors"),
                            div(class = "warn", verbatimTextOutput(outputId = paste0(group_name, "_", arm_name, "_", model_name, "_errors")))
                          )
                        },

                        h4("Model Summary"),
                        verbatimTextOutput(outputId = paste0(group_name, "_", arm_name, "_", model_name, "_summary")),
                        tags$div(style = "height: 100px;")
                      )
                    }))
                  )
                }))
              )
            }
          }))
        )
      })
    })


    # Observe function for dynamic UI rendering
    observe({
      req(model.list$newmod)

      lapply(names(model.list$newmod), function(group_name) {

        if (group_name == "Outcome Model") {
          lapply(names(model.list$newmod[[group_name]]), function(model_name) {

            model <- model.list$newmod[[group_name]][[model_name]]

            # Render formula
            output[[paste0(group_name, "_", model_name, "_formula")]] <- renderPrint({
              cat(paste(deparse(stats::formula(model)), collapse = "\n"))
            })

            # Render summary
            output[[paste0(group_name, "_", model_name, "_summary")]] <- renderPrint({
              summary(model)
            })

            # Render warnings (if exist)
            if (!is.null(model[['warnings']]) && length(model[['warnings']]) > 0) {
              output[[paste0(group_name, "_", model_name, "_warnings")]] <- renderPrint({
                cat(paste(model[['warnings']], collapse = "\n"))
              })
            }

            # Render errors (if exist)
            if (!is.null(model[['errors']]) && length(model[['errors']]) > 0) {
              output[[paste0(group_name, "_", model_name, "_errors")]] <- renderPrint({
                cat(paste(model[['errors']], collapse = "\n"))
              })
            }
          })

        } else {  # For Mediator Models
          lapply(names(model.list$newmod[[group_name]]), function(arm_name) {
            lapply(names(model.list$newmod[[group_name]][[arm_name]]), function(model_name) {

              model <- model.list$newmod[[group_name]][[arm_name]][[model_name]]

              # Render formula
              output[[paste0(group_name, "_", arm_name, "_", model_name, "_formula")]] <- renderPrint({
                cat(paste(deparse(stats::formula(model)), collapse = "\n"))
              })

              # Render summary
              output[[paste0(group_name, "_", arm_name, "_", model_name, "_summary")]] <- renderPrint({
                summary(model)
              })

              # Render warnings (if exist)
              if (!is.null(model[['warnings']]) && length(model[['warnings']]) > 0) {
                output[[paste0(group_name, "_", arm_name, "_", model_name, "_warnings")]] <- renderPrint({
                  cat(paste(model[['warnings']], collapse = "\n"))
                })
              }

              # Render errors (if exist)
              if (!is.null(model[['errors']]) && length(model[['errors']]) > 0) {
                output[[paste0(group_name, "_", arm_name, "_", model_name, "_errors")]] <- renderPrint({
                  cat(paste(model[['errors']], collapse = "\n"))
                })
              }
            })
          })
        }
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
#'  The exposure variable must be categorical, with \code{0} explicitly denoting the unexposed (or control) group, which is taken as the reference group.
#'  Other values represent different, non-reference exposure categories.
#' @param outcome A \code{character} string specifying the name of the outcome variable in the dataset.
#'  The outcome variable can be either binary or continuous.
#' @param mediators A \code{character} vector specifying the names of the variables in the dataset corresponding to the mediators of interest. The mediators
#'  can be either binary or continuous. When estimating the effect type \code{"shift_k_order"},
#'  the order of mediators in the vector is important, as it interpreted as the assumed causal ordering of the mediators.
#' @param intermediate_confs A \code{character} vector specifying the names of the variables in the dataset corresponding to intermediate confounders.
#'  The intermediate confounders can be either binary or continuous. If \code{NULL},
#'  no intermediate confounders are specified.
#' @param confounders A \code{character} vector listing the names of the variables in the dataset corresponding to the baseline confounders.
#' @param interactions_XC A \code{character} string specifying the two-way interactions amongst exposure and baseline confounders
#'  to include in the regression models in the estimation procedure. The default value, \code{"all"},
#'  includes all two-way exposure-confounder interactions but excludes confounder-confounder interactions.
#'  Specify \code{"none"} to exclude all two-way interactions amongst exposure and baseline confounders.
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
  res[1:3] <- list(list())
  names(res) = c("L_models", "M_models", "Outcome Model")

  # Joint of M1 to MK under X=0 and X!=0 ...

  mod_id = 1

  for (k in 1:K) {

    if(k < first){
      res[[names(res)[1]]][[k]] <- catch_model_messages(as.formula(
        gen_formula_shiny(mediators = mediators,
                          exposure = exposure,
                          k = k,
                          interactions_XC = interactions_XC,
                          include_all = TRUE)),
        data = data,
        family = fam_type[[k]])
      res[[names(res)[1]]][[k]]$labels = c(res[[names(res)[1]]][[k]]$labels,
                                           "p_trt / p_ctr",
                                           if(first<K && any(intervention_type %in% c("all", "shift_all"))) "p_all",
                                           if(any(intervention_type %in% c("all", "shift_k"))) paste0("p_", first:K - (first-1)),
                                           if(first<K && any(intervention_type %in% c("all", "shift_k_order"))) paste0("p_", first:(K-1) - (first-1), "_prime"))
    } else if (k >= first){
      res[[names(res)[2]]][[mod_id]] <- catch_model_messages(as.formula(
        gen_formula_shiny(mediators = mediators,
                          exposure = exposure,
                          k = k,
                          interactions_XC = interactions_XC,
                          include_all = TRUE)),
        data = data,
        family = fam_type[[k]])
      res[[names(res)[2]]][[mod_id]]$labels = c(res[[names(res)[2]]][[mod_id]]$labels,
                                                "p_trt / p_ctr",
                                           if((k+1)<=K && any(intervention_type %in% c("all", "shift_k"))) paste0("p_", (k+1):K - (first - 1)),
                                           if(k < (K-1) && (k - first + 2) <= (K - first) && any(intervention_type %in% c("all", "shift_k_order"))) paste0("p_", (k-first+2):(K-first), "_prime"))
      mod_id = mod_id +1
    }
  }


  # Marginals under X=0
  for (k in first:K) {
    res[[names(res)[2]]][[mod_id]] <- catch_model_messages(as.formula(
      gen_formula_shiny(mediators = mediators,
                        exposure = exposure,
                        k = k,
                        interactions_XC = interactions_XC,
                        marginal = TRUE)),
      data = data,
      family = fam_type[[k]])
    if (k == first){
      res[[names(res)[2]]][[mod_id]]$labels = c(if(first<K && any(intervention_type %in% c("all", "shift_all"))) "p_all",
                                                if(any(intervention_type %in% c("all", "shift_k"))) paste0("p_", first - (first - 1)),
                                                if(first<K && any(intervention_type %in% c("all", "shift_k_order"))) paste0("p_", first - (first - 1), "_prime"))
    } else {
      res[[names(res)[2]]][[mod_id]]$labels = c(if(any(intervention_type %in% c("all", "shift_k"))) paste0("p_", k - (first - 1)),
                                                if (k < K && any(intervention_type %in% c("all", "shift_k_order"))) paste0("p_", k - (first - 1), "_prime"))
    }
    mod_id = mod_id +1
  }

  # Joint of others under X!=0
  if (any(intervention_type %in% c("all", "shift_k")) && first<=(K-1)) {
    for (MM in first:(K-1)) {
      index = setdiff(MM:K, MM)
      for (k in index) {
          res[[names(res)[2]]][[mod_id]] <- catch_model_messages(as.formula(
            gen_formula_shiny(mediators = mediators,
                              exposure = exposure,
                              k = k,
                              MM = MM,
                              first = first,
                              K = K,
                              interactions_XC = interactions_XC)),
            data = data,
            family = fam_type[[k]])
          res[[names(res)[2]]][[mod_id]]$labels = paste0("p_", MM - (first - 1))
          mod_id = mod_id +1
      }
    }
  }

  # p_first_prime....p_Kminus1_prime
  # Conditionals under X!=0
  if (any(intervention_type %in% c("all", "shift_k_order"))) {
    for (MM in first:(K - 1)) {
      for (k in (MM + 1):K) {
        res[[names(res)[2]]][[mod_id]] <- catch_model_messages(as.formula(
          gen_formula_shiny(mediators = mediators,
                            exposure = exposure,
                            k = k,
                            interactions_XC = interactions_XC,
                            include_all = TRUE)),
          data = data,
          family = fam_type[[k]])
        res[[names(res)[2]]][[mod_id]]$labels = paste0("p_", MM - (first - 1), "_prime")
        mod_id = mod_id +1
      }
    }
  }
  # For p_all
  # Joint of main ones under X=0
  if (any(intervention_type %in% c("all", "shift_all"))) {
    for (k in (first + 1):K) {
      res[[names(res)[2]]][[mod_id]] <- catch_model_messages(as.formula(
        gen_formula_shiny(mediators = mediators,
                          exposure = exposure,
                          k = k,
                          interactions_XC = interactions_XC,
                          first = first,
                          include_all = TRUE)),
        data = data,
        family = fam_type[[k]])
      res[[names(res)[2]]][[mod_id]]$labels = "p_all"
      mod_id = mod_id + 1
    }
  }

  # outcome
  outcome_type = family_type(data, outcome)
  res[[names(res)[3]]][[1]] <- catch_model_messages(as.formula(paste0(
    outcome, " ~ (", exposure, " + ", paste0(mediators[1:K], collapse = "+"), ")^2+",
    interactions_XC)),
    data = data,
    family = outcome_type[[1]])

  return(res)
}




#' Generate Model Formulas for Shiny-Based Model Assessment
#'
#' This function generates model formulas using the original variable names as they appear in the dataset.
#'
#' @param mediators A \code{character} vector specifying the names of the variables in the dataset corresponding to the mediators of interest. The mediators
#'  can be either binary or continuous. When estimating the effect type \code{"shift_k_order"},
#'  the order of mediators in the vector is important, as it interpreted as the assumed causal ordering of the mediators.
#' @param exposure A \code{character} string specifying the name of the exposure variable in the dataset.
#'  The exposure variable must be categorical, with \code{0} explicitly denoting the unexposed (or control) group, which is taken as the reference group.
#'  Other values represent different, non-reference exposure categories.
#' @param k An integer specifying the index of the mediator for which the formula is being generated.
#' @param first An integer (optional) specifying the index of the first mediator. Defaults to `NULL`,
#'  indicating that the function will generate formulas without explicitly considering this parameter.
#' @param MM An integer (optional) specifying the index of the mediator whose distribution will be shifted.
#'  Defaults to `NULL`, indicating that the function will generate formulas without explicitly considering
#'  this parameter.
#' @param K An integer (optional) specifying the total number of mediators and intermediate confounders.
#'  Defaults to `NULL`, indicating that the function will generate formulas without explicitly considering this
#'  parameter.
#' @param interactions_XC A \code{character} string specifying the two-way interactions amongst exposure and baseline confounders
#'  to include in the regression models in the estimation procedure. The default value, \code{"all"},
#'  includes all two-way exposure-confounder interactions but excludes confounder-confounder interactions.
#'  Specify \code{"none"} to exclude all two-way interactions amongst exposure and baseline confounders.
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



#' Flatten a Nested Model List
#' This function processes a nested list of models by extracting specific sublists (`L_models` and `M_models`),
#' renaming their elements, and merging them into the parent level.
#'
#' @param model_list A nested list containing models, with potential sublists named `L_models` and `M_models`.
#' @keywords internal
flatten_models <- function(model_list) {
  # Check if the input is a list
  if (!is.list(model_list)) return(model_list)

  # Create a copy of the list to modify
  new_list <- model_list

  for (name in names(model_list)) {
    if (is.list(model_list[[name]])) {
      if ("L_models" %in% names(model_list[[name]])) {
        l_models <- model_list[[name]]$L_models
        if (length(l_models) > 0) {
          names(l_models) <- paste0("L", seq_along(l_models), "_model")
          # Merge into the parent level
          new_list[[name]] <- c(new_list[[name]], l_models)
        }
        # Remove L_models
        new_list[[name]]$L_models <- NULL
      }

      if ("M_models" %in% names(model_list[[name]])) {
        m_models <- model_list[[name]]$M_models
        if (length(m_models) > 0) {
          names(m_models) <- paste0("M", seq_along(m_models), "_model")
          # Merge into the parent level
          new_list[[name]] <- c(new_list[[name]], m_models)
        }
        # Remove M_models
        new_list[[name]]$M_models <- NULL
      }

      # Recursively apply function to deeper levels
      new_list[[name]] <- flatten_models(new_list[[name]])
    }
  }
  return(new_list)
}


