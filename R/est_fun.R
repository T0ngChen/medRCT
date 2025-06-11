#' Causal Mediation Analysis Estimating Interventional Effects
#'
#' This function performs the actual causal mediation analysis to estimate interventional effects mapped to a
#' hypothetical target trial.
#'
#' @param dat A \code{data.frame} containing the dataset for analysis.
#' @param ind A \code{vector} of indices specifying the subset of \code{dat} to use for the analysis.
#'  Defaults to all rows of \code{dat}. This parameter is particularly useful when using this function within the
#'  \code{boot()} function from the \code{boot} package, as it enables resampling by specifying subsets of the data.
#' @param first An \code{integer} specifying the index of the first mediator of interest in the combined list of
#'  intermediate confounders and mediators.
#' @param K An \code{integer} specifying the total number of mediators and intermediate confounders.
#'  Mediators are considered sequentially based on their order.
#' @param fam_type A \code{character} string specifying the family type for modeling. Options typically include
#'  \code{"gaussian"} for continuous variables or \code{"binomial"} for binary variables.
#' @param mediators A \code{character} vector including the variable names for mediators (including intermediate
#' confounders).
#' @param interactions_XC A \code{character} string specifying the two-way interactions amongst exposure and baseline confounders
#'  to include in the regression models in the estimation procedure. The default value, \code{"all"},
#'  includes all two-way exposure-confounder interactions but excludes confounder-confounder interactions.
#'  Specify \code{"none"} to exclude all two-way interactions amongst exposure and baseline confounders.
#' @param intervention_type A \code{character} string indicating the type of interventional effect to be estimated.
#'  Options include:
#' \itemize{
#'   \item \code{"all"} (default): Estimates all types of interventional indirect effects.
#'   \item \code{"shift_all"}: Estimates an interventional indirect effect mapped to a target trial assessing the impact of shifting the joint distribution of all
#'    mediators in the exposed to match the corresponding distribution in the unexposed.
#'   \item \code{"shift_k"}: Estimates an interventional indirect effect mapped to a target trial assessing the impact of shifting the distribution of a specific
#'    mediator (\code{k}) in the exposed to match the corresponding distribution in the unexposed.
#'   \item \code{"shift_k_order"}: Estimates an interventional indirect effect mapped to a target trial assessing the impact of shifting the distribution of a
#'    specific mediator (\code{k}) in the exposed to match the corresponding distribution in the unexposed while accounting for the flow-on
#'    effects on causally descendant mediators.
#' }
#' @param effect_measure A \code{character} string specifying the effect measure to calculate.
#'   For a continuous outcome, only \code{"Diff"} (difference in means) is allowed.
#'   For a binary outcome, must be either \code{"RD"} (risk difference) or \code{"RR"} (risk ratio).
#'   If not specified, defaults to \code{"Diff"} for continuous outcomes and \code{"RD"} for binary outcomes.
#'   Only one effect measure can be specified at a time.
#' @param mcsim An \code{integer} specifying the number of Monte Carlo simulations to perform.
#'
#' @importFrom stats as.formula binomial glm predict rbinom rnorm df.residual
#' @importFrom data.table as.data.table ":="
#'
#' @keywords internal
medRCT.fun <- function(
  dat,
  ind = 1:nrow(dat),
  first = first,
  K = K,
  fam_type = fam_type,
  mediators = mediators,
  interactions_XC = interactions_XC,
  intervention_type = intervention_type,
  effect_measure,
  mcsim
) {
  # Take bootstrap sample
  data <- dat[ind, ]

  # Replicate dataset for simulations
  dat2 <- data.table::as.data.table(data)

  dat2[, 1:(2 + K) := lapply(.SD, function(x) NA_integer_), .SDcols = 1:(2 + K)]

  dat2 <- zoo::coredata(dat2)[rep(seq(nrow(dat2)), mcsim), ]
  n <- nrow(dat2)

  # identify the exposure levels
  data$X = as.factor(data$X)
  exposure_level = sort(unique(as.numeric(dat$X)))
  lnzero = exposure_level[exposure_level != 0]

  # ESTIMATE DISTRIBUTIONS
  # Joint of M1 to MK under X=0 and X!=0 ...

  for (k in 1:K) {
    dat2 = joint_dist(
      k = k,
      K = K,
      data = data,
      dat2 = dat2,
      mediators = mediators,
      fam_type = fam_type,
      interactions_XC = interactions_XC,
      exposure_level = exposure_level,
      n = n
    )
  }

  # Estimating the target quantities
  # Marginals under X=0
  for (k in first:K) {
    dat2 <- marg_dist(
      k = k,
      first = first,
      K = K,
      data = data,
      dat2 = dat2,
      mediators = mediators,
      fam_type = fam_type,
      interactions_XC = interactions_XC,
      n = n
    )
  }

  # For p_first,..., p_K
  # Joint of others under X!=0
  if (any(intervention_type %in% c("all", "shift_k")) && first <= (K - 1)) {
    for (MM in first:(K - 1)) {
      index = setdiff(MM:K, MM)
      for (k in index) {
        dat2 <- joint_X_nonzero(
          MM = MM,
          k = k,
          first = first,
          K = K,
          data = data,
          mediators = mediators,
          dat2 = dat2,
          fam_type = fam_type,
          interactions_XC = interactions_XC,
          lnzero = lnzero,
          n = n,
          index = index
        )
      }
    }
  }

  # For p_first_prime,...., p_K_prime
  # Conditionals under X!=0
  if (any(intervention_type %in% c("all", "shift_k_order"))) {
    for (MM in first:(K - 1)) {
      for (k in (MM + 1):K) {
        dat2 <- con_exposed(
          MM = MM,
          k = k,
          K = K,
          data = data,
          dat2 = dat2,
          mediators = mediators,
          fam_type = fam_type,
          interactions_XC = interactions_XC,
          lnzero = lnzero,
          n = n
        )
      }
    }
  }

  # For p_all
  # Joint of main ones under X=0
  if (any(intervention_type %in% c("all", "shift_all"))) {
    for (k in (first + 1):K) {
      dat2 <- joint_unexposed(
        k = k,
        first = first,
        K = K,
        data = data,
        dat2 = dat2,
        mediators = mediators,
        fam_type = fam_type,
        interactions_XC = interactions_XC,
        n = n
      )
    }
  }

  # outcome
  outcome_type = family_type(data, "Y")
  fit <- glm(
    as.formula(paste0(
      "Y~(X+",
      paste0(paste0("M", 1:K), collapse = "+"),
      ")^2+",
      interactions_XC
    )),
    data = data,
    family = outcome_type[[1]]
  )

  if (!fit$converged) {
    warning(paste0(
      "Model did not converge when using the outcome Y as the response"
    ))
  }
  if (any(is.na(fit$coefficients))) {
    na_coefs <- names(coef(fit))[is.na(coef(fit))]
    warning(paste0(
      "The following coefficients were NA: ",
      paste(na_coefs, collapse = ", "),
      "when using the outcome Y as the response"
    ))
  }

  # ESTIMATE OUTCOME EXPECTATION IN EACH ARM & ESTIMATE EFFECTS
  # p_ctr

  a <- 0
  dat2 = set_exposure(data = dat2, column_name = "X", exp_val = a)
  l = 1:K
  dat2[, paste0("M", l) := mget(med_outcome_name(a = a, l = l, K = K))]
  y0 <- predict(fit, newdata = dat2, type = "response")

  results = list()
  p_ctr <- mean(y0)
  results[['p_ctr']] <- p_ctr

  # estimate causal effects
  for (a in lnzero) {
    # TCE and p_trt
    dat2 = set_exposure(data = dat2, column_name = "X", exp_val = a)
    results = compute_assign(
      dat2 = dat2,
      fit = fit,
      a = a,
      K = K,
      first = first,
      type = "trt",
      results = results,
      lnzero = lnzero,
      p_ctr = p_ctr,
      effect_measure = effect_measure
    )

    # p_all and IIE_all
    if (any(intervention_type %in% c("all", "shift_all"))) {
      results = compute_assign(
        dat2 = dat2,
        fit = fit,
        a = a,
        K = K,
        first = first,
        type = "all",
        results = results,
        lnzero = lnzero,
        p_ctr = p_ctr,
        effect_measure = effect_measure
      )
    }
    # p_first....p_K and IIE_first .... IIE_K
    if (any(intervention_type %in% c("all", "shift_k"))) {
      results = compute_assign_loop(
        dat2 = dat2,
        fit = fit,
        a = a,
        K = K,
        first = first,
        p_ctr = p_ctr,
        type = "shift_k",
        results = results,
        lnzero = lnzero,
        effect_measure = effect_measure
      )
    }
    # p_first_prime....p_Kminus1_prime
    if (any(intervention_type %in% c("all", "shift_k_order"))) {
      results = compute_assign_loop(
        dat2 = dat2,
        fit = fit,
        a = a,
        K = K,
        first = first,
        p_ctr = p_ctr,
        type = "shift_k_order",
        results = results,
        lnzero = lnzero,
        effect_measure = effect_measure
      )
    }
  }

  res = unlist(results)
  sorted_names <- sort(names(res))
  IIE_names = sorted_names[grep("IIE", sorted_names)]
  IIE_names = IIE_names[order(nchar(IIE_names))]
  IDE_names = sorted_names[grep("IDE", sorted_names)]
  IDE_names = IDE_names[order(nchar(IDE_names))]
  p_names = sorted_names[!grepl("IIE|TCE|IDE", sorted_names)]
  p_names = p_names[order(nchar(p_names))]

  # Move TCE-related names after IIE-related names
  final_order <- c(
    IIE_names, # All IIE-related names
    sorted_names[grep("TCE", sorted_names)], # All TCE-related names
    IDE_names,
    p_names # Remaining names
  )
  res = res[final_order]

  return(res)
}


medRCT.fun.safe <- function(
  dat,
  ind = 1:nrow(dat),
  first = first,
  K = K,
  fam_type = fam_type,
  mediators = mediators,
  interactions_XC = interactions_XC,
  intervention_type = intervention_type,
  effect_measure,
  mcsim
) {
  res <- tryCatch(
    withCallingHandlers(
      {
        medRCT.fun(
          dat = dat,
          ind = ind,
          first = first,
          K = K,
          fam_type = fam_type,
          mediators = mediators,
          interactions_XC = interactions_XC,
          intervention_type = intervention_type,
          effect_measure = effect_measure,
          mcsim = mcsim
        )
      },
      warning = function(w) {
        stop(w$message)
      }
    ),
    error = function(e) {
      message(e$message)
      NA
    }
  )
  res
}
