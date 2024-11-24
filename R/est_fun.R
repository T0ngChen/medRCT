#' Causal Mediation Analysis for Estimating Interventional Effects
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
#'    effects on its causal descendent mediators.
#' }
#' @param mcsim An \code{integer} specifying the number of Monte Carlo simulations to perform.
#'
#' @importFrom stats as.formula binomial glm predict rbinom rnorm df.residual
#' @importFrom data.table as.data.table ":="
medRCT.fun <- function(dat,
                       ind = 1:nrow(dat),
                       first = first,
                       K = K,
                       fam_type = fam_type,
                       interactions_XC = interactions_XC,
                       intervention_type = intervention_type,
                       mcsim) {
  # Take bootstrap sample
  data <- dat[ind, ]

  # Set flag to capture bootstrap samples to reject
  flag <- FALSE

  # Replicate dataset for simulations
  dat2 <- data.table::as.data.table(data)

  dat2[, 1:(2 + K) := lapply(.SD,
                             function(x) NA_integer_), .SDcols = 1:(2 + K)]

  dat2 <- zoo::coredata(dat2)[rep(seq(nrow(dat2)), mcsim), ]
  n <- nrow(dat2)

  # identify the exposure levels
  data$X = as.factor(data$X)
  exposure_level = sort(unique(as.numeric(dat$X)))
  lnzero = exposure_level[exposure_level!=0]

  # ESTIMATE DISTRIBUTIONS
  # Joint of M1 to MK under X=0 and X!=0 ...

  for (k in 1:K) {
    dat2 = joint_dist(k = k, K = K, data = data, dat2 = dat2,
                      fam_type = fam_type, interactions_XC = interactions_XC,
                      exposure_level = exposure_level, n = n)
  }

  # Estimating the target quantities
  # Marginals under X=0
  for (k in first:K) {
    dat2 <- marg_dist(
      k = k, first = first, K = K, data = data, dat2 = dat2,
      fam_type = fam_type, interactions_XC = interactions_XC, n = n
    )
  }


  # For p_first,..., p_K
  # Joint of others under X!=0
  if (any(intervention_type %in% c("all", "shift_k"))) {
    for (MM in first:K) {
      index = setdiff(first:K, MM)
      for (k in index) {
        dat2 <- joint_X_nonzero(
          MM = MM, k = k, first = first, K = K, data = data,
          dat2 = dat2, fam_type = fam_type, interactions_XC = interactions_XC,
          lnzero = lnzero, n = n, index = index
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
          MM = MM, k = k, K = K, data = data, dat2 = dat2,
          fam_type = fam_type, interactions_XC = interactions_XC,
          lnzero = lnzero, n = n
        )
      }
    }
  }



  # For p_all
  # Joint of main ones under X=0
  if (any(intervention_type %in% c("all", "shift_all"))) {
    for (k in (first + 1):K) {
      dat2 <- joint_unexposed(
        k = k, first = first, K = K, data = data, dat2 = dat2,
        fam_type = fam_type, interactions_XC = interactions_XC, n = n
      )
    }
  }

  # outcome
  outcome_type = family_type(data, "Y")
  fit <- glm(as.formula(paste0(
    "Y~(X+", paste0(paste0("M", 1:K), collapse = "+"), ")^2+",
    interactions_XC)),
    data = data,
    family = outcome_type[[1]])

  if ((!fit$converged) | any(is.na(fit$coefficients)))
    flag <- TRUE


  # ESTIMATE OUTCOME EXPECTATION IN EACH ARM & ESTIMATE EFFECTS
  # p_ctr

  a <- 0

  dat2 = set_exposure(data = dat2, column_name = "X", exp_val = a)

  l = 1:K
  dat2[, paste0("M", l) := mget(med_outcome_name(a = a,
                                                 l = l,
                                                 K = K))]

  y0 <- predict(fit, newdata = dat2, type = "response")

  p_ctr <- mean(y0)


  # p_trt
  for(a in lnzero){

    dat2 = set_exposure(data = dat2, column_name = "X", exp_val = a)

    dat2[, paste0("M", l) := mget(med_outcome_name(a = a,
                                                   l = l,
                                                   K = K))]

    y1 <- predict(fit, newdata = dat2, type = "response")

    if(length(lnzero) > 1){
      assign(paste0("p_trt_", a), mean(y1))
      assign(paste0("TCE_",a), get(paste0("p_trt_", a)) - p_ctr)
    } else {
      p_trt <- mean(y1)
      TCE <- p_trt - p_ctr
    }
  }


  # p_all
  if (any(intervention_type %in% c("all", "shift_all"))) {
    for(a in lnzero){

      dat2 = set_exposure(data = dat2, column_name = "X", exp_val = a)

      if (first > 1) {
        l = 1:(first - 1)
        dat2[, paste0("M", l) := mget(med_outcome_name(a = a,
                                                       l = l,
                                                       K = K))]
      }

      # all mediators of interest
      k = first:K
      dat2[, paste0("M", k) := mget(med_outcome_all(l=k,
                                                    first = first,
                                                    a = 0,
                                                    K=K))]

      y1 <- predict(fit, newdata = dat2, type = "response")

      if(length(lnzero) > 1){
        assign(paste0("p_all_",a), mean(y1))
        # IIE
        assign(paste0("IIE_all_",a), get(paste0("p_trt_", a))-get(paste0("p_all_", a)))
      } else {
        p_all <- mean(y1)
        # IIE
        IIE_all <- p_trt - p_all
      }
    }
  }


  # p_first....p_K
  if (any(intervention_type %in% c("all", "shift_k"))) {
    for(a in lnzero){

      dat2 = set_exposure(data = dat2, column_name = "X", exp_val = a)

      if (first > 1) {
        l = 1:(first - 1)
        dat2[, paste0("M", l) :=  mget(med_outcome_name(a = a,
                                                        l = l,
                                                        K = K))]
      }

      for (MM in first:K) {
        dat2[, paste0("M", MM) := get(paste0("m", MM, "_", 0, "_",
                                             strrep("m", K)))]
        k = setdiff(first:K, MM)
        dat2[, paste0("M", k) := mget(med_joint_other(k = k, a = a, MM = MM, K = K,
                                                      ordering = FALSE))]

        y0 <- predict(fit, newdata = dat2, type = "response")

        if(length(lnzero) > 1){
          assign(paste0("p_", MM, "_", a), mean(y0))
          # IIE
          assign(paste0("IIE_", MM, "_", a), get(paste0("p_trt_", a)) - get(paste0("p_", MM, "_", a)))
        } else {
          assign(paste0("p_", MM), mean(y0))
          # IIE
          assign(paste0("IIE_", MM), p_trt - get(paste0("p_", MM)))
        }
      }
    }
  }


  # p_first_prime....p_Kminus1_prime
  if (any(intervention_type %in% c("all", "shift_k_order"))) {
    for(a in lnzero){

      dat2 = set_exposure(data = dat2, column_name = "X", exp_val = a)

      for (MM in first:(K - 1)) {
        if (MM != 1) {
          l = 1:(MM - 1)

          dat2[, paste0("M", l) := mget(med_outcome_name(a = a,
                                                         l = l,
                                                         K = K))]
        }

        dat2[, paste0("M", MM) := get(paste0("m", MM, "_", 0, "_",
                                             strrep("m", K)))]

        if ((MM + 1) <= K) {
          k = (MM + 1):K
          dat2[, paste0("M", k) := mget(med_joint_other(k = k, a = a, MM = MM, K = K))]
        }

        y0 <- predict(fit, newdata = dat2, type = "response")

        if(length(lnzero) > 1){
          assign(paste0("p_", MM, "_prime", "_", a), mean(y0))
          # IIE
          assign(paste0("IIE_", MM,"_prime_", a), get(paste0("p_trt_", a)) - get(paste0("p_", MM, "_prime_", a)))
        } else {
          assign(paste0("p_", MM, "_prime"), mean(y0))
          # IIE
          assign(paste0("IIE_", MM, "_prime"), p_trt - get(paste0("p_", MM, "_prime")))
        }
      }
    }
  }




  # Collect and return results
  res <- vector()
  res_names <- vector()


  # save results for IIE_k
  if (any(intervention_type %in% c("all", "shift_k"))) {
    if(length(lnzero) > 1){
      res <- c(res, unlist(mget(paste0("IIE_",
                                       rep(first:K, length(lnzero)),
                                       "_",
                                       rep(lnzero, each = length(first:K))))))
      res_names <- c(res_names,
                     paste0(
                       "IIE_",
                       rep(first:K - (first - 1), length(lnzero)), "_",
                       rep(lnzero, each = length(first:K)),
                       " (p_trt_",
                       rep(lnzero, each = length(first:K)),
                       " - p_",
                       rep(first:K - (first - 1), length(lnzero)), "_",
                       rep(lnzero, each = length(first:K)),
                       ")"
                     ))
    } else {
      res <- c(res, unlist(mget(paste0("IIE_", first:K))))
      res_names <- c(res_names,
                     paste0(
                       "IIE_",
                       first:K - (first - 1),
                       " (p_trt - p_",
                       first:K - (first - 1),
                       ")"
                     ))
    }
  }

  # save results for IIE_k_prime
  if (any(intervention_type %in% c("all", "shift_k_order"))) {
    if(length(lnzero) > 1){
      res <- c(res, unlist(mget(paste0("IIE_",
                                       rep(first:(K-1), length(lnzero)),
                                       "_prime_",
                                       rep(lnzero, each = length(first:(K-1)))))))
      res_names <- c(res_names,
                     paste0(
                       "IIE_",
                       rep(first:(K - 1) - (first - 1), length(lnzero)), "_",
                       rep(lnzero, each = length(first:(K - 1))), "_prime",
                       " (p_trt_",
                       rep(lnzero, each = length(first:(K - 1))),
                       " - p_",
                       rep(first:(K - 1) - (first - 1), length(lnzero)), "_",
                       rep(lnzero, each = length(first:(K - 1))),
                       "_prime)"
                     ))
    } else {

      res <- c(res, unlist(mget(paste0("IIE_", first:(K - 1), "_prime"))))
      res_names <- c(res_names,
                     paste0(
                       "IIE_",
                       first:(K - 1) - (first - 1),
                       "_prime",
                       " (p_trt - p_",
                       first:(K - 1) - (first - 1),
                       "_prime)"
                     ))
    }
  }

  # save results for IIE_k_all

  if (any(intervention_type %in% c("all", "shift_all"))) {
    if(length(lnzero) > 1){
      res <- c(res, unlist(mget(paste0("IIE_all_", lnzero))))
      res_names <- c(res_names, paste0("IIE_all_", lnzero,
                                       " (p_trt_", lnzero, " - p_all_", lnzero, ")"))

    } else {
      res <- c(res, IIE_all)
      res_names <- c(res_names, "IIE_all (p_trt - p_all)")
    }
  }

  # TCE
  # p_trt & p_ctr
  if(length(lnzero) > 1){
    res <- c(res, unlist(mget(paste0("TCE_", lnzero))),
             unlist(mget(paste0("p_trt_", lnzero))), p_ctr)
    res_names <- c(res_names, paste0("TCE_", lnzero, " (p_trt_", lnzero, " - p_ctr)"),
                   paste0("p_trt_", lnzero), "p_ctr")
  } else {
    res <- c(res, TCE, p_trt, p_ctr)
    res_names <- c(res_names, "TCE (p_trt - p_ctr)", "p_trt", "p_ctr")
  }

  # p_k
  if (any(intervention_type %in% c("all", "shift_k"))) {
    if(length(lnzero) > 1){
      res <- c(res, unlist(mget(paste0("p_",
                                       rep(first:K, length(lnzero)),
                                       "_",
                                       rep(lnzero, each = length(first:K))))))
      res_names <- c(res_names,
                     paste0("p_",
                            rep(first:K - (first - 1), length(lnzero)),
                            "_",
                            rep(lnzero, each = length(first:K))))
    } else {
      res <- c(res, unlist(mget(paste0("p_", first:K))))
      res_names <- c(res_names, paste0("p_", first:K - (first - 1)))
    }
  }

  # p_k_prime
  if (any(intervention_type %in% c("all", "shift_k_order"))) {
    if(length(lnzero) > 1){
      res <- c(res, unlist(mget(paste0("p_",
                                       rep(first:(K-1), length(lnzero)),
                                       "_prime_",
                                       rep(lnzero, each = length(first:(K-1)))))))

      res_names <- c(res_names,
                     paste0("p_",
                            rep(first:(K-1) - (first - 1), length(lnzero)),
                            "_prime_",
                            rep(lnzero, each = length(first:(K-1)))))

    } else {
      res <- c(res, unlist(mget(paste0("p_", first:(K - 1), "_prime"))))
      res_names <- c(res_names, paste0("p_", first:(K - 1) - (first - 1), "_prime"))
    }
  }

  # p_all
  if (any(intervention_type %in% c("all", "shift_all"))) {
    if(length(lnzero) > 1){
      res <- c(res, unlist(mget(paste0("p_all_", lnzero))))
      res_names <- c(res_names, paste0("p_all_", lnzero))
    } else {
      res <- c(res, p_all)
      res_names <- c(res_names, "p_all")
    }
  }
  names(res) = res_names

  if (!flag)
    return(res)
  else
    return(rep(NA, length(res)))
}
