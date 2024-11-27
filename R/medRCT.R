utils::globalVariables(".SD")


#' Causal mediation analysis using interventional effects that map to a target trial
#'
#' 'medRCT' is used to estimate interventional effects that map to a target trial evaluating hypothetical mediator interventions of interest.
#' It can handle multiple mediators, including some not of primary interest but that are exposure-induced
#' mediator-outcome confounders.
#'
#' @param dat A \code{data.frame} containing the dataset for analysis. It should include variables for the exposure,
#'  outcome, mediators, confounders, and exposure-induced mediator-outcome confounders specified in the analysis.
#' @param exposure  A \code{character} string specifying the name of the exposure variable in the dataset.
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
#' @param confounders  A \code{character} vector listing the names of baseline confounders.
#' @param interactions_XC A \code{character} string specifying the exposure-confounder or confounder-confounder
#'  interaction terms to include in the regression models for confounder adjustment. The default value, \code{"all"},
#'  includes all two-way exposure-confounder interactions but excludes confounder-confounder interactions.
#'  Specify \code{"none"} to exclude all two-way exposure-confounder and confounder-confounder interactions.
#' @param intervention_type A \code{character} string indicating the type of interventional effect to be estimated.
#'  Options include:
#' \itemize{
#'   \item \code{"all"} (default): Estimates all types of interventional indirect effects.
#'   \item \code{"shift_all"}: Estimates the interventional indirect effect of shifting the joint distribution of all
#'    mediators in the exposed to match the level in the unexposed.
#'   \item \code{"shift_k"}: Estimates the interventional indirect effect of shifting the distribution of a specific
#'    mediator (\code{k}) in the exposed to match the level in the unexposed.
#'   \item \code{"shift_k_order"}: Estimates the interventional indirect effect of shifting the distribution of a
#'    specific mediator (\code{k}) in the exposed to match the level in the unexposed while accounting for the flow-on
#'    effects on its causally descendent mediators.
#' }
#' @param mcsim An \code{integer} specifying the number of Monte Carlo simulations to perform.
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
#' @export
medRCT <- function(dat,
                   exposure,
                   outcome,
                   mediators,
                   intermediate_confs,
                   confounders,
                   interactions_XC = "all",
                   intervention_type = c("all", "shift_all", "shift_k", "shift_k_order"),
                   mcsim,
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





