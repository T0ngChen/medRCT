
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
#' @param interactions_XC  A \code{character} string specifying the exposure-confounder or
#'  confounder-confounder interaction terms to include in the regression models for
#'  confounder adjustment.
#' @param exposure_level A numeric vector specifying the levels of the exposure
#'  (e.g., \code{c(0, 1)}) for which counterfactual predictions are performed.
#' @param n An integer specifying the number of observations for `dat2`.
joint_dist <- function(k, K, data, dat2, fam_type,
                       interactions_XC, exposure_level, n) {
  # Fit the model
  fit <- glm(
    as.formula(gen_formula(k = k,
                           interactions_XC = interactions_XC,
                           include_all = TRUE)),
    data = data,
    family = fam_type[[k]]
  )

  # Check convergence and coefficients
  flag <- (!fit$converged) | any(is.na(fit$coefficients))

  # Loop through exposure levels
  for (a in exposure_level) {
    # Set exposure
    dat2 <- set_exposure(data = dat2, column_name = "X", exp_val = a)

    # Handle mediators M1 to M(k-1) if k > 1
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
#' @param interactions_XC A \code{character} string specifying the exposure-confounder or
#'  confounder-confounder interaction terms to include in the regression models for
#'  confounder adjustment.
#' @param n An integer specifying the number of observations for `dat2`.
marg_dist <- function(k, first, K, data, dat2, fam_type, interactions_XC, n) {
  # Fit the model
  fit <- glm(
    as.formula(gen_formula(k = k,
                           interactions_XC = interactions_XC,
                           marginal = TRUE)),
    data = data,
    family = fam_type[[k]]
  )

  # Initialize a flag for convergence issues
  flag <- (!fit$converged) | any(is.na(fit$coefficients))

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




#' Estimation of Joint Distribution of Mediators under Exposed
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
#' @param interactions_XC A \code{character} string specifying the exposure-confounder or
#'  confounder-confounder interaction terms to include in the regression models for
#'  confounder adjustment.
#' @param lnzero A numeric vector specifying the non-zero levels of the exposure.
#' @param n An integer specifying the number of observations for `dat2`.
#' @param index An integer vector specifying the indices of all mediators, excluding the mediator
#'  specified by `MM`.
joint_X_nonzero <- function(MM, k, first, K, data, dat2,
                            fam_type, interactions_XC, lnzero, n, index) {
  # Check for intermediate confounders
  if (first == 1) {
    if (MM == 1 && k == index[1]) {
      fit <- glm(as.formula(gen_formula(k = k, MM = MM,
                                        first = first, K = K,
                                        interactions_XC = interactions_XC)),
                 data = data, family = fam_type[[k]])
    } else if (MM != 1 && k == index[1]) {
      return(dat2)
    } else {
      fit <- glm(as.formula(gen_formula(k = k, MM = MM,
                                        first = first, K = K,
                                        interactions_XC = interactions_XC)),
                 data = data, family = fam_type[[k]])
    }
  } else {
    fit <- glm(as.formula(gen_formula(k = k, MM = MM,
                                      first = first, K = K,
                                      interactions_XC = interactions_XC)),
               data = data, family = fam_type[[k]])
  }

  # Check model convergence and coefficients
  if ((!fit$converged) | any(is.na(fit$coefficients))) {
    flag <- TRUE
  }

  # loop over exposure levels
  for (a in lnzero) {
    dat2 <- set_exposure(data = dat2, column_name = "X", exp_val = a)

    # Update mediators
    if ((first == 1 && k != index[1]) || first != 1) {
      l <- setdiff(1:(k - 1), MM)
      dat2[, paste0("M", l) := mget(med_outcome_name(a = a, l = l, K = K))]
    }

    # Perform counterfactual prediction
    var_name <- paste0("m", k, "_", a, "_", paste0(c(
      rep(paste0(a), min(k - 1, MM - 1)),
      "m",
      rep(paste0(a), max(k - 1 - MM, 0)),
      rep("m", K - 1 - min(k - 1, MM - 1) - max(k - 1 - MM, 0))
    ), collapse = ""))

    dat2 <- cf_predict(fit = fit, data = dat2,
                       var_name = var_name, n = n,
                       family = fam_type[[k]]$family)
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
#' @param interactions_XC A \code{character} string specifying the exposure-confounder or
#'  confounder-confounder interaction terms to include in the regression models for
#'  confounder adjustment.
#' @param lnzero A numeric vector specifying the non-zero levels of the exposure.
#' @param n An integer specifying the number of observations for `dat2`.
con_exposed <- function(MM, k, K, data, dat2, fam_type, interactions_XC,
                        lnzero, n) {
  # Fit the model for the mediator k
  fit <- glm(
    as.formula(gen_formula(k = k,
                           interactions_XC = interactions_XC,
                           include_all = TRUE)),
    data = data,
    family = fam_type[[k]]
  )

  # Check for convergence or NA coefficients
  flag <- (!fit$converged) | any(is.na(fit$coefficients))

  # Iterate over exposure levels
  for (a in lnzero) {
    # Set exposure
    dat2 <- set_exposure(data = dat2, column_name = "X", exp_val = a)

    # Handle mediators before MM
    if (MM != 1) {
      l = 1:(MM - 1)
      dat2[, paste0("M", l) := mget(med_outcome_name(a = a,
                                                     l = l,
                                                     K = K))]
    }

    # Handle mediator MM
    dat2[, paste0("M", MM) := get(paste0("m", MM, "_", 0, "_", strrep("m", K)))]

    # Handle mediators between MM and k
    if (k > (MM + 1)) {
      for (l in (MM + 1):(k - 1)) {
        dat2[, paste0("M", l) := get(paste0("m", l, "_", a, "_", paste0(c(
          rep(paste0(a), MM - 1),
          0,
          rep(paste0(a), max(l - 1 - MM, 0)),
          rep("m", K - MM - max(l - 1 - MM, 0))
        ), collapse = "")))]
      }
    }

    # Perform counterfactual prediction
    var_name <- paste0("m", k, "_", a, "_", paste0(c(
      rep(paste0(a), MM - 1),
      0,
      rep(paste0(a), max(k - 1 - MM, 0)),
      rep("m", K - MM - max(k - 1 - MM, 0))
    ), collapse = ""))

    dat2 <- cf_predict(
      fit = fit,
      data = dat2,
      var_name = var_name,
      n = n,
      family = fam_type[[k]]$family
    )
  }
  dat2
}




#' Estimation of Joint Distribution of Mediators under Unexposed
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
#' @param interactions_XC A \code{character} string specifying the exposure-confounder or
#'  confounder-confounder interaction terms to include in the regression models for
#'  confounder adjustment.
#' @param n An integer specifying the number of observations for `dat2`.
joint_unexposed <- function(k, first, K, data, dat2, fam_type, interactions_XC, n) {
  # Fit the model for the mediator k
  fit <- glm(
    as.formula(gen_formula(k = k, interactions_XC = interactions_XC,
                           include_all = TRUE, first = first)),
    data = data,
    family = fam_type[[k]]
  )

  # Check for convergence or NA coefficients
  flag <- (!fit$converged) | any(is.na(fit$coefficients))

  # Set exposure
  a <- 0
  dat2 <- set_exposure(data = dat2, column_name = "X", exp_val = a)

  # Update mediators
  l = first:(k - 1)
  dat2[, paste0("M", l) := mget(med_outcome_all(l = l,
                                                first = first,
                                                a = a,
                                                K = K))]

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



