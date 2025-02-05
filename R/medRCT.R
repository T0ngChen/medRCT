utils::globalVariables(".SD")


#' Causal Mediation Analysis Estimating Interventional Effects Mapped to A Target Trial
#'
#' `medRCT` is used to estimate interventional effects that map to a target trial evaluating hypothetical mediator interventions
#' of interest. It can handle any number of potentially correlated mediators, including mediators that are not of primary
#' interest but that are intermediate (exposure-induced) mediator-outcome confounders.
#'
#' @param dat A \code{data.frame} containing the dataset for analysis. It should include variables for the exposure,
#'  outcome, mediators, (baseline) confounders, and intermediate (exposure-induced) mediator-outcome confounders specified in the analysis.
#' @param exposure  A \code{character} string specifying the name of the exposure variable in the dataset.
#'  The exposure variable must be categorical, with \code{0} explicitly denoting the unexposed (or control) group, which is taken as the reference group.
#'  Other values represent different, non-reference exposure categories.
#' @param outcome A \code{character} string specifying the name of the outcome variable in the dataset.
#'  The outcome variable can be binary or continuous.
#' @param mediators A \code{character} vector specifying the names of the variables in the dataset corresponding to the mediators of interest. The mediators
#'  can be binary or continuous. When estimating the effect type \code{"shift_k_order"},
#'  the order of mediators in the vector is important and must correspond to the assumed causal ordering of the mediators.
#' @param intermediate_confs A \code{character} vector specifying the names of the variables in the dataset corresponding to intermediate confounders.
#'  The intermediate confounders can be binary or continuous. If \code{NULL},
#'  no intermediate confounders are specified.
#' @param confounders  A \code{character} vector listing the names of the variables in the dataset corresponding to the baseline confounders.
#' @param interactions_XC A \code{character} string specifying the two-way interactions between exposure and baseline confounders
#'  to include in the regression models in the estimation procedure. The default value, \code{"all"},
#'  includes all two-way exposure-confounder interactions but excludes confounder-confounder interactions.
#'  Specify \code{"none"} to exclude all two-way interactions between exposure and baseline confounders. See Vignette for further details.
#' @param intervention_type A \code{character} string indicating the type of interventional effect to be estimated.
#'  Options include:
#' \itemize{
#'   \item \code{"all"} (default): Estimates all three types of interventional indirect effects.
#'   \item \code{"shift_all"}: Estimates an interventional indirect effect mapped to a target trial assessing the impact of shifting the joint distribution of all
#'    mediators in the exposed, given baseline confounders, to match the corresponding distribution in the unexposed.
#'   \item \code{"shift_k"}: Estimates an interventional indirect effect mapped to a target trial assessing the impact of shifting the distribution of a specific
#'    mediator (\code{k}) in the exposed, given baseline confounders, to match the corresponding distribution in the unexposed.
#'   \item \code{"shift_k_order"}: Estimates an interventional indirect effect mapped to a target trial assessing the impact of shifting the distribution of a
#'    specific mediator (\code{k}) in the exposed, given baseline confounders, to match the corresponding distribution in the unexposed while accounting for the flow-on
#'    effects on causally descendant mediators.
#' }
#' @param mcsim An \code{integer} specifying the number of Monte Carlo simulations to perform. The default is 200.
#' It is recommended to run analysis with no fewer than 200 Monte Carlo simulations.
#' @param bootstrap A \code{logical} value indicating whether bootstrapping should be performed. If \code{TRUE}
#'  (default), bootstrapping is conducted using the \code{boot} function from the \code{boot} package.
#' @param boot_args A \code{list} of arguments for bootstrapping. The default settings are:
#' \itemize{
#'   \item \code{R}: Number of bootstrap replicates (default: 100).
#'   \item \code{stype}: Specifies the statistic type passed to the \code{boot} function (default: \code{"i"}).
#'   \item \code{ci.type}: Specifies the type of confidence interval to compute (default: \code{"norm"}).
#' }
#' @param ... Additional arguments passed to the \code{boot} function from the \code{boot} package.
#'
#' @details
#'
#' Before conducting the mediation analysis, users are encouraged to assess the models fitted by the algorithm
#' using the interactive Shiny application, which can be launched by running the function \code{medRCT_shiny}.
#' The Shiny app provides a user-friendly interface to review model summaries and identify potential warnings and errors,
#' ensuring that the models are appropriately specified before proceeding with the analysis.
#'
#' If issues with model fitting are detected, users are encouraged to adjust the exposure-confounder interaction term as needed.
#' However, \strong{mediators or confounders must not be selected based on model fitting results.}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Define confounders for the analysis
#' confounders <- c("child_sex", "child_atsi", "mat_cob", "mat_engl", "mat_age")
#'
#' # Define intermediate confounders
#' intermediate_confs <- "fam_stress"
#'
#' # Estimate interventional indirect effects for various
#' # hypothetical interventions
#' med_res <- medRCT(
#'   dat = LSACdata,
#'   exposure = "sep",
#'   outcome = "child_mh",
#'   mediators = c("parent_mh", "preschool_att"),
#'   intermediate_confs = intermediate_confs,
#'   confounders = confounders,
#'   interactions_XC = "all",
#'   intervention_type = "all",
#'   bootstrap = TRUE,
#'   boot_args = list(R = 100, stype = "i", ci.type = "norm"),
#'   mcsim = 100
#')
#' # Summarise the results
#' summary(med_res)
#' }
medRCT <- function(dat,
                   exposure,
                   outcome,
                   mediators,
                   intermediate_confs,
                   confounders,
                   interactions_XC = "all",
                   intervention_type = c("all", "shift_all", "shift_k", "shift_k_order"),
                   mcsim = 200,
                   bootstrap = TRUE,
                   boot_args = list(R = 100, stype = "i", ci.type = "norm"),
                   ...) {
  # match intervention type
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

  if(any(intervention_type %in% c("all", "shift_k_order"))) {
    message(paste0("Assumed causal order for estimating effect of type 'shift_k_order': ",  paste(mediators, collapse = ", "), "\n"))
  }


  mediators = c(intermediate_confs, mediators)

  # define the first mediator of interest
  first = length(intermediate_confs) + 1

  K <- length(mediators)

  dat <- as.data.frame(dat)
  # count the No. of missing values
  no.miss = nrow(dat) - sum(stats::complete.cases(dat))
  message(paste0(
    "Conducting complete case analysis, ",
    no.miss,
    " observations were excluded due to missing data.\n"
  ))

  if (mcsim<200) {
    message("Note: It is recommended to run analysis with no fewer than 200 Monte Carlo simulations.")
  }

  dat <- dat[stats::complete.cases(dat), ]

  fam_type = family_type(dat, mediators)

  # Rename all variables & prepare dataset
  dat$X <- dat[, exposure]
  dat$Y <- dat[, outcome]
  for (k in 1:K)
    dat[, paste0("M", k)] <- dat[, mediators[k]]
  dat <- dat[, c("X", paste0("M", 1:K), "Y", confounders)]

  # Prepare confounder terms for formulae
  # (defaults to all exposure-confounder interactions if not provided)
  if (interactions_XC == "all") {
    interactions_XC <- paste(paste(rep("X", length(confounders)), confounders, sep =
                                     "*"), collapse = "+")
  } else if (interactions_XC == "none") {
    interactions_XC <- paste(confounders, collapse = "+")
  } else {
    interactions_XC <- gsub(exposure, "X", interactions_XC)
  }

  # set R to 1 if bootstrap is not required
  if (bootstrap == FALSE) {
    boot_args$R = 1
  }

  # bootstrap
  boot.out <- boot::boot(
    data = dat,
    statistic = medRCT.fun,
    first = first,
    K = K,
    mediators = mediators,
    mcsim = mcsim,
    fam_type = fam_type,
    interactions_XC = interactions_XC,
    intervention_type = intervention_type,
    stype = boot_args$stype,
    R = boot_args$R,
    ...
  )

  # grab results from bootstrap
  if (bootstrap == TRUE) {
    est = boot.out$t0
    se = apply(boot.out$t, 2, stats::sd)
    pval <- 2 * (1 - stats::pnorm(q = abs(est / se)))
    cilow = ciupp = numeric()
    for (i in 1:length(boot.out$t0)) {
      ci.type = ifelse(is.null(boot_args$ci.type), "norm", boot_args$ci.type)
      bt <- boot::boot.ci(boot.out, index = i, type = ci.type)
      if (ci.type == "perc") {
        cilow <- c(cilow, bt$percent[4])
        ciupp <- c(ciupp, bt$percent[5])
      } else if (ci.type == "norm") {
        cilow <- c(cilow, bt$normal[2])
        ciupp <- c(ciupp, bt$normal[3])
      }
    }

    # save results
    out = list(
      est = est,
      se = se,
      pval = pval,
      cilow = cilow,
      ciupp = ciupp,
      sample.size = nrow(dat),
      mcsim = mcsim,
      bootstrap = bootstrap
    )
  } else {
    out = list(est = boot.out$t0,
               bootstrap = bootstrap)
  }

  class(out) <- "medRCT"
  out
}





