#' Estimation of Joint Distribution of Mediators and Intermediate Confounders under
#' Exposed and Unexposed
#'
#' @param k An integer specifying the index of the current mediator being processed.
#' @param K An integer specifying the total number of mediators and intermediate confounders.
#' @param data A \code{data.frame} containing the dataset for analysis. The dataset is
#'  used to fit the model for the mediators.
#' @param dat2 A \code{data.frame} used to perform counterfactual predictions and
#'  random draws.
#' @param fam_type A \code{character} string specifying the family type for modeling.
#'  Options typically include \code{"gaussian"} for continuous variables or \code{"binomial"}
#'  for binary variables.
#' @param mediators A \code{character} vector including the variable names for mediators
#' (including intermediate confounders).
#' @param interactions_XC  A \code{character} string specifying the two-way interactions amongst exposure and baseline confounders
#'  to include in the regression models in the estimation procedure.
#' @param exposure_level A numeric vector specifying the levels of the exposure
#'  (e.g., \code{c(0, 1)}) for which counterfactual predictions are performed.
#' @param separation_method Method to handle separation, only relevant for binomial (binary outcome) models.
#'   Options are \code{"brglm"} (Logistic regression models are fitted using bias reduction methods for generalised linear models implemented in the \code{brglm2} package)
#'   or \code{"discard"} (if separation is detected, the function returns \code{NA}. If this occurs during the main estimation,
#'   the program stops; if it occurs during bootstrapping, the affected bootstrap samples are discarded).
#' @param n An integer specifying the number of observations for `dat2`.
#'
#' @importFrom stats coef as.formula
#'
#' @keywords internal
joint_dist <- function(
  k,
  K,
  data,
  dat2,
  fam_type,
  mediators,
  interactions_XC,
  exposure_level,
  separation_method,
  n
) {
  # Fit the model
  fit <- fit_model(
    stats::as.formula(gen_formula(
      k = k,
      interactions_XC = interactions_XC,
      include_all = TRUE
    )),
    data = data,
    family = fam_type[[k]],
    separation_method = separation_method
  )

  # Check convergence and coefficients
  if (!fit$converged) {
    warning(paste0(
      "Model did not converge when using variable ",
      mediators[k],
      " as the response"
    ))
  }
  if (any(is.na(fit$coefficients))) {
    na_coefs <- names(stats::coef(fit))[is.na(stats::coef(fit))]
    warning(paste0(
      "The following coefficients were NA: ",
      paste(na_coefs, collapse = ", "),
      "when using variable ",
      mediators[k],
      " as the response"
    ))
  }

  # Loop through exposure levels
  for (a in exposure_level) {
    # Set exposure
    dat2 <- set_exposure(data = dat2, column_name = "X", exp_val = a)

    # Handle mediators M1 to M_(k-1) if k > 1
    if (k != 1) {
      l <- 1:(k - 1)
      dat2[, paste0("M", l) := mget(med_outcome_name(a = a, l = l, K = K))]
    }

    # Perform counterfactual prediction and random draw
    dat2 <- cf_predict(
      fit = fit,
      data = dat2,
      var_name = med_outcome_name(a = a, l = k, K = K),
      n = n,
      family = fam_type[[k]]$family
    )
  }
  dat2
}


#' Estimation of Marginal Distributions of Mediators under Unexposed
#'
#' @param k An integer specifying the index of the current mediator being processed.
#' @param first An integer specifying the index of the first mediator.
#' @param K An integer specifying the total number of mediators and intermediate confounders.
#' @param data A \code{data.frame} containing the dataset for analysis. The dataset is
#'  used to fit the model for the mediators.
#' @param dat2 A \code{data.frame} used to perform counterfactual predictions and
#'  random draws.
#' @param fam_type A \code{character} string specifying the family type for modeling.
#'  Options typically include \code{"gaussian"} for continuous variables or \code{"binomial"}
#'  for binary variables.
#' @param mediators A \code{character} vector including the variable names for mediators
#'  (including intermediate confounders).
#' @param interactions_XC A \code{character} string specifying the two-way interactions amongst exposure and baseline confounders
#'  to include in the regression models in the estimation procedure.
#' @param separation_method Method to handle separation, only relevant for binomial (binary outcome) models.
#'   Options are \code{"brglm"} (Logistic regression models are fitted using bias reduction methods for generalised linear models implemented in the \code{brglm2} package)
#'   or \code{"discard"} (if separation is detected, the function returns \code{NA}. If this occurs during the main estimation,
#'   the program stops; if it occurs during bootstrapping, the affected bootstrap samples are discarded).
#' @param n An integer specifying the number of observations for `dat2`.
#'
#' @importFrom stats coef as.formula
#'
#' @keywords internal
marg_dist <- function(
  k,
  first,
  K,
  data,
  dat2,
  fam_type,
  mediators,
  interactions_XC,
  separation_method,
  n
) {
  # Fit the model
  fit <- fit_model(
    stats::as.formula(gen_formula(
      k = k,
      interactions_XC = interactions_XC,
      marginal = TRUE
    )),
    data = data,
    family = fam_type[[k]],
    separation_method = separation_method
  )

  if (!fit$converged) {
    warning(paste0(
      "Model did not converge when using variable ",
      mediators[k],
      " as the response"
    ))
  }
  if (any(is.na(fit$coefficients))) {
    na_coefs <- names(stats::coef(fit))[is.na(stats::coef(fit))]
    warning(paste0(
      "The following coefficients were NA: ",
      paste(na_coefs, collapse = ", "),
      "when using variable ",
      mediators[k],
      " as the response"
    ))
  }

  # Set exposure to 0
  a <- 0
  dat2 <- set_exposure(data = dat2, column_name = "X", exp_val = a)

  # Perform counterfactual prediction and random draw
  dat2 <- cf_predict(
    fit = fit,
    data = dat2,
    var_name = paste0("m", k, "_", a, "_", strrep("m", K)),
    n = n,
    family = fam_type[[k]]$family
  )

  dat2 = dat2
}


#' Estimation of Joint Distribution of All Other Mediators under Exposed
#'
#' @param MM An integer specifying the index of the mediator whose distribution
#'  will be shifted.
#' @param k An integer specifying the index of the current mediator being processed.
#' @param first An integer specifying the index of the first mediator.
#' @param K An integer specifying the total number of mediators and intermediate confounders.
#' @param data A \code{data.frame} containing the dataset for analysis. The dataset is
#'  used to fit the model for the mediators.
#' @param dat2 A \code{data.frame} used to perform counterfactual predictions and
#'  random draws.
#' @param fam_type A \code{character} string specifying the family type for modeling.
#'  Options typically include \code{"gaussian"} for continuous variables or \code{"binomial"}
#'  for binary variables.
#' @param mediators A \code{character} vector including the variable names for mediators
#'  (including intermediate confounders).
#' @param interactions_XC A \code{character} string specifying the two-way interactions amongst exposure and baseline confounders
#'  to include in the regression models in the estimation procedure.
#' @param lnzero A numeric vector specifying the non-zero levels of the exposure.
#' @param separation_method Method to handle separation, only relevant for binomial (binary outcome) models.
#'   Options are \code{"brglm"} (Logistic regression models are fitted using bias reduction methods for generalised linear models implemented in the \code{brglm2} package)
#'   or \code{"discard"} (if separation is detected, the function returns \code{NA}. If this occurs during the main estimation,
#'   the program stops; if it occurs during bootstrapping, the affected bootstrap samples are discarded).
#' @param n An integer specifying the number of observations for `dat2`.
#' @param index An integer vector specifying the indices of all mediators, excluding the mediator
#'  specified by `MM`.
#'
#' @importFrom stats coef as.formula
#'
#' @keywords internal
joint_X_nonzero <- function(
  MM,
  k,
  first,
  K,
  data,
  dat2,
  fam_type,
  mediators,
  interactions_XC,
  lnzero,
  separation_method,
  n,
  index
) {
  # Check for intermediate confounders
  fit <- fit_model(
    stats::as.formula(gen_formula(
      k = k,
      MM = MM,
      first = first,
      K = K,
      interactions_XC = interactions_XC
    )),
    data = data,
    family = fam_type[[k]],
    separation_method = separation_method
  )

  # Check model convergence and coefficients
  if (!fit$converged) {
    warning(paste0(
      "Model did not converge when using variable ",
      mediators[k],
      " as the response"
    ))
  }
  if (any(is.na(fit$coefficients))) {
    na_coefs <- names(stats::coef(fit))[is.na(stats::coef(fit))]
    warning(paste0(
      "The following coefficients were NA: ",
      paste(na_coefs, collapse = ", "),
      "when using variable ",
      mediators[k],
      " as the response"
    ))
  }

  # loop over exposure levels
  for (a in lnzero) {
    dat2 <- set_exposure(data = dat2, column_name = "X", exp_val = a)

    # Update mediators
    if (first != 1) {
      l = 1:(first - 1)
      dat2[, paste0("M", l) := mget(med_outcome_name(a = a, l = l, K = K))]
    }
    if ((k - 1) > first) {
      l = setdiff(first:(k - 1), MM)
      dat2[,
        paste0("M", l) := mget(med_joint_other(
          k = l,
          a = a,
          MM = MM,
          K = K,
          ordering = FALSE
        ))
      ]
    }

    # Perform counterfactual prediction
    dat2 <- cf_predict(
      fit = fit,
      data = dat2,
      var_name = med_joint_other(
        k = k,
        a = a,
        MM = MM,
        K = K,
        ordering = FALSE
      ),
      n = n,
      family = fam_type[[k]]$family
    )
  }
  dat2
}


#' Estimation of Conditionals of Mediators under Exposed
#'
#' @param MM An integer specifying the index of the mediator whose distribution
#'  will be shifted.
#' @param k An integer specifying the index of the current mediator being processed.
#' @param K An integer specifying the total number of mediators and intermediate confounders.
#' @param data A \code{data.frame} containing the dataset for analysis. The dataset is
#'  used to fit the model for the mediators.
#' @param dat2 A \code{data.frame} used to perform counterfactual predictions and
#'  random draws.
#' @param fam_type A \code{character} string specifying the family type for modeling.
#'  Options typically include \code{"gaussian"} for continuous variables or \code{"binomial"}
#'  for binary variables.
#' @param mediators A \code{character} vector including the variable names for mediators
#'  (including intermediate confounders).
#' @param interactions_XC A \code{character} string specifying the two-way interactions amongst exposure and baseline confounders
#'  to include in the regression models in the estimation procedure.
#' @param lnzero A numeric vector specifying the non-zero levels of the exposure.
#' @param separation_method Method to handle separation, only relevant for binomial (binary outcome) models.
#'   Options are \code{"brglm"} (Logistic regression models are fitted using bias reduction methods for generalised linear models implemented in the \code{brglm2} package)
#'   or \code{"discard"} (if separation is detected, the function returns \code{NA}. If this occurs during the main estimation,
#'   the program stops; if it occurs during bootstrapping, the affected bootstrap samples are discarded).
#' @param n An integer specifying the number of observations for `dat2`.
#'
#' @importFrom stats coef as.formula
#'
#' @keywords internal
con_exposed <- function(
  MM,
  k,
  K,
  data,
  dat2,
  fam_type,
  mediators,
  interactions_XC,
  lnzero,
  separation_method,
  n
) {
  # Fit the model for the mediator k
  fit <- fit_model(
    stats::as.formula(gen_formula(
      k = k,
      interactions_XC = interactions_XC,
      include_all = TRUE
    )),
    data = data,
    family = fam_type[[k]],
    separation_method = separation_method
  )

  # Check for convergence or NA coefficients
  if (!fit$converged) {
    warning(paste0(
      "Model did not converge when using variable ",
      mediators[k],
      " as the response"
    ))
  }
  if (any(is.na(fit$coefficients))) {
    na_coefs <- names(stats::coef(fit))[is.na(stats::coef(fit))]
    warning(paste0(
      "The following coefficients were NA: ",
      paste(na_coefs, collapse = ", "),
      "when using variable ",
      mediators[k],
      " as the response"
    ))
  }

  # Iterate over exposure levels
  for (a in lnzero) {
    # Set exposure
    dat2 <- set_exposure(data = dat2, column_name = "X", exp_val = a)

    # Handle mediators before MM
    if (MM != 1) {
      l = 1:(MM - 1)
      dat2[, paste0("M", l) := mget(med_outcome_name(a = a, l = l, K = K))]
    }

    # Handle mediator MM
    dat2[, paste0("M", MM) := get(paste0("m", MM, "_", 0, "_", strrep("m", K)))]

    # Handle mediators between MM and k
    if (k > (MM + 1)) {
      l = (MM + 1):(k - 1)
      dat2[,
        paste0("M", l) := mget(med_joint_other(k = l, a = a, MM = MM, K = K))
      ]
    }

    # Perform counterfactual prediction

    dat2 <- cf_predict(
      fit = fit,
      data = dat2,
      var_name = med_joint_other(k = k, a = a, MM = MM, K = K),
      n = n,
      family = fam_type[[k]]$family
    )
  }
  dat2
}


#' Estimation of Joint Distribution of All Mediators under Unexposed
#'
#' @param k An integer specifying the index of the current mediator being processed.
#' @param first An integer specifying the index of the first mediator.
#' @param K An integer specifying the total number of mediators and intermediate confounders.
#' @param data A \code{data.frame} containing the dataset for analysis. The dataset is
#'  used to fit the model for the mediators.
#' @param dat2 A \code{data.frame} used to perform counterfactual predictions and
#'  random draws.
#' @param fam_type A \code{character} string specifying the family type for modeling.
#'  Options typically include \code{"gaussian"} for continuous variables or \code{"binomial"}
#'  for binary variables.
#' @param mediators A \code{character} vector including the variable names for mediators
#'  (including intermediate confounders).
#' @param interactions_XC A \code{character} string specifying the two-way interactions amongst exposure and baseline confounders
#'  to include in the regression models in the estimation procedure.
#' @param separation_method Method to handle separation, only relevant for binomial (binary outcome) models.
#'   Options are \code{"brglm"} (Logistic regression models are fitted using bias reduction methods for generalised linear models implemented in the \code{brglm2} package)
#'   or \code{"discard"} (if separation is detected, the function returns \code{NA}. If this occurs during the main estimation,
#'   the program stops; if it occurs during bootstrapping, the affected bootstrap samples are discarded).
#' @param n An integer specifying the number of observations for `dat2`.
#'
#' @importFrom stats coef as.formula
#'
#' @keywords internal
joint_unexposed <- function(
  k,
  first,
  K,
  data,
  dat2,
  fam_type,
  mediators,
  interactions_XC,
  separation_method,
  n
) {
  # Fit the model for the mediator k
  fit <- fit_model(
    stats::as.formula(gen_formula(
      k = k,
      interactions_XC = interactions_XC,
      include_all = TRUE,
      first = first
    )),
    data = data,
    family = fam_type[[k]],
    separation_method = separation_method
  )

  # Check for convergence or NA coefficients
  if (!fit$converged) {
    warning(paste0(
      "Model did not converge when using variable ",
      mediators[k],
      " as the response"
    ))
  }
  if (any(is.na(fit$coefficients))) {
    na_coefs <- names(stats::coef(fit))[is.na(stats::coef(fit))]
    warning(paste0(
      "The following coefficients were NA: ",
      paste(na_coefs, collapse = ", "),
      "when using variable ",
      mediators[k],
      " as the response"
    ))
  }

  # Set exposure
  a <- 0
  dat2 <- set_exposure(data = dat2, column_name = "X", exp_val = a)

  # Update mediators
  l = first:(k - 1)
  dat2[,
    paste0("M", l) := mget(med_outcome_all(l = l, first = first, a = a, K = K))
  ]

  # Perform counterfactual prediction and random draw
  dat2 <- cf_predict(
    fit = fit,
    data = dat2,
    var_name = med_outcome_all(l = k, first = first, a = a, K = K),
    n = n,
    family = fam_type[[k]]$family
  )

  dat2
}


# estimation IIE, IDE, and return results
compute_assign = function(
  dat2,
  fit,
  a,
  K,
  first,
  type,
  lnzero,
  p_ctr,
  results,
  effect_measure
) {
  if (type == "trt") {
    l = 1:K
    dat2[, paste0("M", l) := mget(med_outcome_name(a = a, l = l, K = K))]
  } else if (type == "all") {
    if (first > 1) {
      l = 1:(first - 1)
      dat2[, paste0("M", l) := mget(med_outcome_name(a = a, l = l, K = K))]
    }

    # all mediators of interest
    k = first:K
    dat2[,
      paste0("M", k) := mget(med_outcome_all(
        l = k,
        first = first,
        a = 0,
        K = K
      ))
    ]
  }
  y1 <- predict(fit, newdata = dat2, type = "response")
  avg_pred = mean(y1)
  if (type == "trt") {
    if (length(lnzero) > 1) {
      results[[paste0("p_trt_", a)]] <- avg_pred
      if (effect_measure %in% c("RD", "Diff")) {
        results[[paste0("TCE_", a, " (p_trt_", a, " - p_ctr)")]] <-
          avg_pred - p_ctr
      } else if (effect_measure == "RR") {
        results[[paste0("TCE_", a, " log(p_trt_", a, " / p_ctr)")]] <-
          log(avg_pred / p_ctr)
      }
    } else {
      results$p_trt <- avg_pred
      if (effect_measure %in% c("RD", "Diff")) {
        results[["TCE (p_trt - p_ctr)"]] <- avg_pred - p_ctr
      } else if (effect_measure == "RR") {
        results[["TCE log(p_trt / p_ctr)"]] <- log(avg_pred / p_ctr)
      }
    }
  } else if (type == "all") {
    if (length(lnzero) > 1) {
      results[[paste0("p_", type, "_", a)]] <- avg_pred
      # IIE
      if (effect_measure %in% c("RD", "Diff")) {
        results[[paste0(
          "IIE_",
          type,
          "_",
          a,
          " (p_trt_",
          a,
          " - p_all_",
          a,
          ")"
        )]] =
          results[[paste0("p_trt_", a)]] - avg_pred
      } else if (effect_measure == "RR") {
        results[[paste0(
          "IIE_",
          type,
          "_",
          a,
          " log(p_trt_",
          a,
          " / p_all_",
          a,
          ")"
        )]] =
          log(results[[paste0("p_trt_", a)]] / avg_pred)
      }

      # IDE
      if (effect_measure %in% c("RD", "Diff")) {
        results[[paste0("IDE_", type, "_", a, " (p_all_", a, " - p_ctr)")]] =
          avg_pred - p_ctr
      } else if (effect_measure == "RR") {
        results[[paste0("IDE_", type, "_", a, " log(p_all_", a, " / p_ctr)")]] =
          log(avg_pred / p_ctr)
      }
    } else {
      results[[paste0("p_", type)]] <- avg_pred
      # IIE
      if (effect_measure %in% c("RD", "Diff")) {
        results[[paste0(
          "IIE_",
          type,
          " (p_trt - p_",
          type,
          ")"
        )]] <- results$p_trt - avg_pred
      } else if (effect_measure == "RR") {
        results[[paste0(
          "IIE_",
          type,
          " log(p_trt / p_",
          type,
          ")"
        )]] <- log(results$p_trt / avg_pred)
      }
      # IDE
      if (effect_measure %in% c("RD", "Diff")) {
        results[[paste0("IDE_", type, " (p_", type, " - p_ctr)")]] <-
          avg_pred - p_ctr
      } else if (effect_measure == "RR") {
        results[[paste0("IDE_", type, " log(p_", type, " / p_ctr)")]] <-
          log(avg_pred / p_ctr)
      }
    }
  }
  results
}


compute_assign_loop = function(
  dat2,
  fit,
  a,
  K,
  first,
  type,
  lnzero,
  results,
  p_ctr,
  effect_measure
) {
  # get index for loop
  if (type == "shift_k") {
    index = first:K
    suffix = NULL
  } else if (type == "shift_k_order") {
    index = first:(K - 1)
    suffix = "_prime"
  }
  for (MM in index) {
    # prepare data
    if (type == "shift_k") {
      if (first > 1) {
        l = 1:(first - 1)
        dat2[, paste0("M", l) := mget(med_outcome_name(a = a, l = l, K = K))]
      }

      dat2[,
        paste0("M", MM) := get(paste0("m", MM, "_", 0, "_", strrep("m", K)))
      ]
      if (length(first:K) > 1) {
        k = setdiff(first:K, MM)
        dat2[,
          paste0("M", k) := mget(med_joint_other(
            k = k,
            a = a,
            MM = MM,
            K = K,
            ordering = FALSE
          ))
        ]
      }
    } else if (type == "shift_k_order") {
      if (MM != 1) {
        l = 1:(MM - 1)
        dat2[, paste0("M", l) := mget(med_outcome_name(a = a, l = l, K = K))]
      }
      dat2[,
        paste0("M", MM) := get(paste0("m", MM, "_", 0, "_", strrep("m", K)))
      ]
      if ((MM + 1) <= K) {
        k = (MM + 1):K
        dat2[,
          paste0("M", k) := mget(med_joint_other(k = k, a = a, MM = MM, K = K))
        ]
      }
    }
    y1 <- predict(fit, newdata = dat2, type = "response")
    avg_pred = mean(y1)

    if (length(lnzero) > 1) {
      results[[paste0("p_", MM - (first - 1), "_", a, suffix)]] <- avg_pred
      # IIE
      if (effect_measure %in% c("RD", "Diff")) {
        results[[paste0(
          "IIE_",
          MM - (first - 1),
          "_",
          a,
          suffix,
          " (p_trt_",
          a,
          " - p_",
          MM - (first - 1),
          "_",
          a,
          suffix,
          ")"
        )]] =
          results[[paste0("p_trt_", a)]] - avg_pred
      } else if (effect_measure == "RR") {
        results[[paste0(
          "IIE_",
          MM - (first - 1),
          "_",
          a,
          suffix,
          " log(p_trt_",
          a,
          " / p_",
          MM - (first - 1),
          "_",
          a,
          suffix,
          ")"
        )]] =
          log(results[[paste0("p_trt_", a)]] / avg_pred)
      }
      # IDE
      if (effect_measure %in% c("RD", "Diff")) {
        results[[paste0(
          "IDE_",
          MM - (first - 1),
          "_",
          a,
          suffix,
          " (p_",
          MM - (first - 1),
          "_",
          a,
          suffix,
          " - p_ctr)"
        )]] =
          avg_pred - p_ctr
      } else if (effect_measure == "RR") {
        results[[paste0(
          "IDE_",
          MM - (first - 1),
          "_",
          a,
          suffix,
          " log(p_",
          MM - (first - 1),
          "_",
          a,
          suffix,
          " / p_ctr)"
        )]] =
          log(avg_pred / p_ctr)
      }
    } else {
      results[[paste0("p_", MM - (first - 1), suffix)]] <- avg_pred
      # IIE
      if (effect_measure %in% c("RD", "Diff")) {
        results[[paste0(
          "IIE_",
          MM - (first - 1),
          suffix,
          " (p_trt - p_",
          MM - (first - 1),
          suffix,
          ")"
        )]] <-
          results$p_trt - avg_pred
      } else if (effect_measure == "RR") {
        results[[paste0(
          "IIE_",
          MM - (first - 1),
          suffix,
          " log(p_trt / p_",
          MM - (first - 1),
          suffix,
          ")"
        )]] <-
          log(results$p_trt / avg_pred)
      }
      # IDE
      if (effect_measure %in% c("RD", "Diff")) {
        results[[paste0(
          "IDE_",
          MM - (first - 1),
          suffix,
          " (p_",
          MM - (first - 1),
          suffix,
          " - p_ctr)"
        )]] <-
          avg_pred - p_ctr
      } else if (effect_measure == "RR") {
        results[[paste0(
          "IDE_",
          MM - (first - 1),
          suffix,
          " log(p_",
          MM - (first - 1),
          suffix,
          " / p_ctr)"
        )]] <-
          log(avg_pred / p_ctr)
      }
    }
  }

  results
}
