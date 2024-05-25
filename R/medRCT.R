utils::globalVariables(".SD")


#' Causal mediation analysis for estimating the interventional effect
#'
#' 'medRCT' is used to estimating the interventional effects that are defined by a mapping to a target trial. It can handle multiple
#' mediators, including some not of primary interest but that are intermediate confounders
#'
#' @param dat a data.frame with the data for analysis
#' @param exposure a character string representing the name of the exposure variable in the data. The exposure variable must be binary, with \code{1} indicating exposed or treated, and \code{0} indicating unexposed or control.
#' @param outcome a character string representing the name of the outcome variable in the data. The outcome variable must be binary.
#' @param mediators a character vector representing the names of the mediators in the data. The mediators must be binary variables.
#' @param intermediate_confs a character vector representing the names of the intermediate confounders in the data. The intermediate
#' confounders must be binary variables.
#' @param confounders a character vector representing the names of the confounders in the data, which must be of the required class
#' (e.g. factor if appropriate)
#' @param interactions_XC a character string specifying the exposure-confounder or confounder-confounder interaction terms
#' to include in all regression models in the procedure. Default is 'all', which includes all two-way exposure-confounder interactions
#' and excludes confounder-confounder interactions.
#' @param intervention_type a character string indicating the type of the interventional effect to be estimated.
#' Can be 'all', 'shift_all', 'shift_k', 'shift_k_order'. Default is 'all', under which all effects will be estimated.
#' ‘shift all’ is to estimate the interventional effect under shifting the joint distribution of all mediatios.
#' ‘shift k’ is to estimate the interventional effect under shifting the distribution of mediator k, independently of other mediators.
#' ‘shift_k_order’ is to estimate the interventional effect under shifting the distribution of mediator k, independently of the its
#' antecedent mediators, while allowing for the flow on effect of mediator k on its descendent mediators.
#' @param mcsim the number of Monte Carlo simulations to conduct
#' @param boostrap logical. If \code{TRUE}, bootstrap will be conducted.
#' @param boot_args a \code{list} of bootstrapping arguments. \code{R} is the number of bootstrap replicates.
#' \code{stype} indicates what the second argument of \code{statistics} in the \code{boot} function represents
#' @param ... other arguments passed to the \code{boot} function in the \code{boot} package.
#'
#' @export
medRCT <- function(dat, exposure, outcome, mediators, intermediate_confs, confounders,
                   interactions_XC = "all",
                   intervention_type = c("all", "shift_all", "shift_k", "shift_k_order"), mcsim,
                   boostrap = TRUE,
                   boot_args = list(R = 100, stype = "i"), ...) {
  intervention_type = sapply(intervention_type, function(arg) match.arg(arg, choices = c("all", "shift_all", "shift_k", "shift_k_order")))
  ci.type = "perc"
  mediators = c(intermediate_confs, mediators)
  first = length(intermediate_confs) + 1
  K <- length(mediators)
  dat <- as.data.frame(dat)
  # Rename all variables & prepare dataset
  dat$X <- dat[, exposure]
  dat$Y <- dat[, outcome]
  for (k in 1:K)
    dat[, paste("M", k, sep = "")] <- dat[, mediators[k]]
  dat <- dat[, c("X", paste("M", 1:K, sep = ""), "Y", confounders)]

  # count the No. of missing values
  no.miss = nrow(dat) - sum(stats::complete.cases(dat))
  message(paste0("Conducting complete case analysis, ", no.miss, " observations deleted\n"))
  dat <- dat[stats::complete.cases(dat),]
  # Prepare confounder terms for formulae
  # (defaults to all exposure-confounder interactions if not provided)
  if (interactions_XC == "all"){
    interactions_XC <- paste(paste(rep("X", length(confounders)), confounders, sep ="*"),
                       collapse = "+")
  } else {
    interactions_XC <- gsub(exposure, "X", interactions_XC)
  }

  if (boostrap == FALSE) {
    boot_args$R = 1
  }

  boot.out <- boot::boot(
    data = dat,
    statistic = medRCT.fun,
    first = first,
    K = K,
    mcsim = mcsim,
    interactions_XC = interactions_XC,
    intervention_type = intervention_type,
    stype = boot_args$stype,
    R = boot_args$R,
    ...
  )

  if (boostrap == TRUE) {
    est = boot.out$t0
    se = apply(boot.out$t, 2, stats::sd)
    pval <- 2 * (1 - stats::pnorm(q = abs(est / se)))
    cilow = ciupp = numeric()
    for (i in 1:length(boot.out$t0)) {
      bt <- boot::boot.ci(boot.out, index = i, type = ci.type)
      cilow <- c(cilow, bt$percent[4])
      ciupp <- c(ciupp, bt$percent[5])
    }
    out = list(
      est = est,
      se = se,
      pval = pval,
      cilow = cilow,
      ciupp = ciupp,
      sample.size = nrow(dat),
      mcsim = mcsim,
      boostrap = boostrap
    )
  } else {
    out = list(
      est = boot.out$t0,
      boostrap = boostrap
    )
  }

  class(out) <- "medRCT"
  out
}


#' Causal mediation analysis for estimating the interventional effect mapped to a target trial
#'
#' @param dat A data.frame with the data for analysis
#' @param ind A vector of indices that define the sample from dat on which to conduct the analysis. Defaults to all rows of the data.
#' This facilitates the use of this function within the boot() function from the boot package
#' @param first index of first mediator of interest after combining intermediate confounders with mediators
#' @param K the number of mediators
#' @param interactions_XC a character string specifying the exposure-confounder or confounder-confounder interaction terms
#' to include in all regression models in the procedure. Defaults to include all two-way exposure-confounder interactions
#' and no confounder-confounder interactions.
#' @param intervention_type a character string indicating the type of the interventional effect to be estimated.
#' Can be 'all', 'shift_all', 'shift_k', 'shift_k_order'. Default is 'all'.
#' @param mcsim the number of Monte Carlo simulations to conduct
#'
#' @importFrom stats as.formula binomial glm predict rbinom
#' @importFrom data.table as.data.table ":="
medRCT.fun <- function(dat,
                       ind = 1:nrow(dat),
                       first = first,
                       K = K,
                       interactions_XC = interactions_XC,
                       intervention_type = intervention_type,
                       mcsim) {
  # Take boostrap sample
  data <- dat[ind, ]

  # Set flag to capure bootstrap samples to reject
  flag <- FALSE

  # Replicate dataset for simulations
  dat2 <- data.table::as.data.table(data)
  dat2[, 1:(2 + K) := lapply(.SD, function(x) NA_integer_), .SDcols = 1:(2 + K)]

  dat2 <- zoo::coredata(dat2)[rep(seq(nrow(dat2)), mcsim), ]
  n <- nrow(dat2)

  # ESTIMATE DISTRIBUTIONS
  # Joint of M1 to MK under X=0 and X=1

  for (k in 1:K) {
    if (k == 1)
      fit <- glm(as.formula(paste("M", k, "~X+", interactions_XC, sep = "")),
                 data = data, family = binomial)
    else
      fit <- glm(as.formula(paste("M", k, "~(X+",
                                  paste(paste("M", 1:(k - 1), sep = ""), collapse = "+"), ")^2+",
                                  interactions_XC, sep = "")),
                 data = data, family = binomial)
    if ((!fit$converged) | any(is.na(fit$coefficients)))
      flag <- TRUE

    for (a in c(0, 1)) {
      dat2$X <- a

      if (k != 1) {
        for (l in 1:(k - 1))
          dat2[, paste("M", l, sep = "") := get(
            paste("m", l, "_", a, "_",
                  paste(c(rep(paste(a), (l - 1)), rep("m", K - (l - 1))), collapse = ""), sep = ""))]
      }

      dat2[, paste("m", k, "_", a, "_",
                   paste(c(rep(paste(a), (k - 1)), rep("m", K - (k - 1))), collapse = ""),
                   sep = "") := rbinom(n, 1, predict(fit, newdata = dat2, type = "response"))]
    }
  }

  # Estimating the target quantities
  # Marginals under X=0
  for (k in first:K) {
    fit <- glm(as.formula(paste("M", k, "~X+", interactions_XC, sep = "")), data =
                 data, family = binomial)

    if ((!fit$converged) | any(is.na(fit$coefficients)))
      flag <- TRUE

    a <- 0
    dat2[, 'X' := a]
    dat2[, paste("m", k, "_", a, "_", paste(rep("m", K), collapse = ""),
                 sep = "") := rbinom(n, 1, predict(fit, newdata = dat2, type = "response"))]
  }


  # For p_first,..., p_K
  # Joint of others under X=1
  if (any(intervention_type %in% c("all", "shift_k"))) {
    for (MM in first:K) {
      for (k in setdiff(first:K, MM)) {
        fit <- glm(as.formula(paste("M", k, "~(X+",
                                    paste(paste("M", setdiff(1:(k - 1), MM), sep = ""),
                                          collapse = "+"),")^2+", interactions_XC, sep = "")),
                   data = data, family = binomial)

        if ((!fit$converged) | any(is.na(fit$coefficients)))
          flag <- TRUE

        a <- 1
        dat2[, 'X' := a]

        if (k != 1) {
          for (l in setdiff(1:(k - 1), MM))
            dat2[, paste("M", l, sep = "") := get(
              paste("m", l, "_", a, "_", paste(c(rep(paste(a), (l - 1)), rep("m", K - (l - 1))),
                                               collapse = ""), sep = ""))]
        }
        dat2[, paste("m", k, "_", a, "_",
                     paste(c(rep(paste(a), min(k - 1, MM - 1)), "m", rep(paste(a), max(k - 1 - MM, 0)),
                             rep("m", K - 1 - min(k - 1, MM - 1) - max(k - 1 - MM, 0))), collapse = ""),
                     sep = "") := rbinom(n, 1, predict(fit, newdata = dat2, type = "response"))]
      }
    }
  }


  # For p_first_prime,...., p_K_prime
  # Conditionals under X=1
  if (any(intervention_type %in% c("all", "shift_k_order"))) {
    for (MM in first:(K - 1)) {
      for (k in (MM + 1):K) {
        fit <- glm(as.formula(paste("M", k, "~(X+",
                                    paste(paste("M", 1:(k - 1), sep = ""), collapse = "+"),
                                    ")^2+", interactions_XC, sep = "")),
                   data = data, family = binomial)

        if ((!fit$converged) | any(is.na(fit$coefficients)))
          flag <- TRUE

        a <- 1
        dat2[, 'X' := a]

        if (MM != 1) {
          for (l in 1:(MM - 1))
            dat2[, paste("M", l, sep = "") := get(
              paste("m", l, "_", a, "_", paste(c(rep(paste(a), (l - 1)), rep("m", K - (l - 1))),
                                               collapse = ""), sep = ""))]
        }
        dat2[, paste("M", MM, sep = "") := get(paste("m", MM, "_", 0, "_",
                                                      paste(rep("m", K), collapse = ""), sep = ""))]
        if (k > (MM + 1)) {
          for (l in (MM + 1):(k - 1))
            dat2[, paste("M", l, sep = "") := get(
              paste("m", l, "_", a, "_", paste(c(rep(paste(a), MM - 1), 0, rep(paste(a), max(l - 1 - MM, 0)),
                                                 rep("m", K - MM - max(l - 1 - MM, 0))), collapse = ""),
                    sep = ""))]
        }

        dat2[, paste("m", k, "_", a, "_", paste(c(rep(paste(a), MM - 1), 0, rep(paste(a), max(k - 1 - MM, 0)),
                                                  rep("m", K - MM - max(k - 1 - MM, 0))), collapse = ""),
                      sep = "") := rbinom(n, 1, predict(fit, newdata = dat2, type = "response"))]
      }
    }
  }



  # For p_all
  # Joint of main ones under X=0
  if (any(intervention_type %in% c("all", "shift_all"))) {
    for (k in (first + 1):K) {
      fit <- glm(as.formula(paste("M", k, "~(X+",
                                  paste(paste("M", first:(k - 1), sep = ""), collapse = "+"),
                                  ")^2+", interactions_XC, sep = "")),
                 data = data, family = binomial)

      if ((!fit$converged) | any(is.na(fit$coefficients)))
        flag <- TRUE

      a <- 0
      dat2[, 'X' := a]

      for (l in first:(k - 1)) {
        dat2[, paste("M", l, sep = "") := get(
          paste("m", l, "_", a, "_", paste(c(rep("m", first - 1), rep(paste(a), (l - first)),
                                             rep("m", K - l + 1)), collapse = ""), sep = ""))]
      }

      dat2[, paste("m", k, "_", a, "_",
                   paste(c(rep("m", first - 1), rep(paste(a), k - first), rep("m", K + 1 - k)),
                         collapse = ""), sep = "") := rbinom(n, 1, predict(fit, newdata = dat2, type = "response"))]
    }
  }


  # outcome
  # Y

  fit <- glm(as.formula(
    paste("Y~(X+", paste(paste("M", 1:K, sep = ""), collapse = "+"), ")^2+", interactions_XC, sep ="")),
    data = data, family = binomial)
  if ((!fit$converged) | any(is.na(fit$coefficients)))
    flag <- TRUE


  # ESTIMATE OUTCOME EXPECTATION IN EACH ARM & ESTIMATE EFFECTS
  # p_ctr

  a <- 0
  dat2[, 'X' := a]
  for (k in 1:K) {
    dat2[, paste("M", k, sep = "") := get(
      paste("m", k, "_", a, "_", paste(c(rep(paste(a), (k - 1)), rep("m", K - (k - 1))),
                                       collapse = ""), sep = ""))]
  }

  y0 <- predict(fit, newdata = dat2, type = "response")

  p_ctr <- mean(y0)


  # p_trt

  a <- 1
  dat2[, 'X' := a]
  for (k in 1:K) {
    dat2[, paste("M", k, sep = "") := get(
      paste("m", k, "_", a, "_", paste(c(rep(paste(a), (k - 1)), rep("m", K - (k - 1))),
                                       collapse = ""), sep = ""))]
  }

  y1 <- predict(fit, newdata = dat2, type = "response")

  p_trt <- mean(y1)


  # p_all
  if (any(intervention_type %in% c("all", "shift_all"))) {
    a <- 1
    dat2[, 'X' := a]
    for (k in 1:(first - 1)) {
      dat2[, paste("M", k, sep = "") := get(
        paste("m", k, "_", a, "_", paste(c(rep(paste(a), (k - 1)), rep("m", K - (k - 1))),
                                         collapse = ""), sep = ""))]
    }

    a <- 0
    for (k in first:K) {
      dat2[, paste("M", k, sep = "") := get(
        paste("m", k, "_", a, "_", paste(c(rep("m", first - 1), rep(paste(a), k - first), rep("m", K + 1 - k)),
                                         collapse = ""), sep = ""))]
    }

    y1 <- predict(fit, newdata = dat2, type = "response")

    p_all <- mean(y1)
    # Interventional effects
    IIE_all <- p_trt - p_all
  }


  # p_first....p_K
  if (any(intervention_type %in% c("all", "shift_k"))) {
    a <- 1
    dat2[, 'X' := a]

    for (k in 1:(first - 1)) {
      dat2[, paste("M", k, sep = "") :=  get(
        paste("m", k, "_", a, "_", paste(c(rep(paste(a), (k - 1)), rep("m", K - (k - 1))),
                                         collapse = ""), sep = ""))]
    }

    for (MM in first:K) {
      dat2[, paste("M", MM, sep = "") := get(
        paste("m", MM, "_", 0, "_", paste(paste(rep("m", K), sep = ""), collapse = ""), sep = ""))]

      for (k in setdiff(first:K, MM)) {
        dat2[, paste("M", k, sep = "") := get(
          paste("m", k, "_", a, "_", paste(c(rep(paste(a), min(k - 1, MM - 1)), "m", rep(paste(a), max(k - 1 - MM, 0)),
                                             rep("m", K - 1 - min(k - 1, MM - 1) - max(k - 1 - MM, 0))), collapse = ""), sep = ""))]
      }

      y0 <- predict(fit, newdata = dat2, type = "response")

      assign(paste("p_", MM, sep = ""), mean(y0))
    }
    # Interventional effects
    for (k in first:K)
      assign(paste("IIE_", k, sep = ""), p_trt - get(paste("p_", k, sep ="")))
  }


  # p_first_prime....p_Kminus1_prime
  if (any(intervention_type %in% c("all", "shift_k_order"))) {
    a <- 1
    dat2$X <- a

    for (MM in first:(K - 1)) {
      if (MM != 1) {
        for (k in 1:(MM - 1)) {
          dat2[, paste("M", k, sep = "") := get(
            paste("m", k, "_", a, "_", paste(c(rep(paste(a), (k - 1)), rep("m", K - (k - 1))),
                                             collapse = ""), sep = ""))]
        }
      }

      dat2[, paste("M", MM, sep = "") := get(
        paste("m", MM, "_", 0, "_", paste(paste(rep("m", K), sep = ""), collapse = ""), sep = ""))]

      if ((MM + 1) <= K) {
        for (k in (MM + 1):K) {
          dat2[, paste("M", k, sep = "") := get(
            paste("m", k, "_", a, "_", paste(c(rep(paste(a), MM - 1), 0, rep(paste(a), max(k - 1 - MM, 0)),
                                               rep("m", K - MM - max(k - 1 - MM, 0))), collapse = ""), sep = ""))]
        }
      }

      y0 <- predict(fit, newdata = dat2, type = "response")

      assign(paste("p_", MM, "_prime", sep = ""), mean(y0))
    }
    # Interventional effects
    for (k in first:(K - 1))
      assign(paste("IIE_", k, "_prime", sep = ""), p_trt - get(paste("p_", k, "_prime", sep ="")))
  }




  # Collect and return results
  res <- vector()
  res_names <- vector()

  if (any(intervention_type %in% c("all", "shift_k"))) {
    res <- c(res, unlist(lapply(first:K, function(k) get(paste("IIE_", k, sep = "")))))
    res_names <- c(res_names, paste0("IIE_", first:K - (first - 1),
                                     " (p_trt - p_", first:K - (first - 1), ")"))
  }

  if (any(intervention_type %in% c("all", "shift_k_order"))) {
    res <- c(res, unlist(lapply(first:(K - 1), function(k) get(paste("IIE_", k, "_prime", sep = "")))))
    res_names <- c(res_names, paste0("IIE_", first:(K - 1) - (first - 1), "_prime",
                                     " (p_trt - p_", first:(K - 1) - (first - 1), "_prime)"))
  }
  if (any(intervention_type %in% c("all", "shift_all"))) {
    res <- c(res, IIE_all)
    res_names <- c(res_names, "IIE_all (p_trt - p_all)")
  }
  res <- c(res, p_trt, p_ctr)
  res_names <- c(res_names, "p_trt", "p_ctr")
  if (any(intervention_type %in% c("all", "shift_k"))) {
    res <- c(res, unlist(lapply(first:K, function(k) get(paste("p_", k, sep = "")))))
    res_names <- c(res_names, paste("p_", first:K - (first - 1), sep = ""))
  }

  if (any(intervention_type %in% c("all", "shift_k_order"))) {
    res <- c(res, unlist(lapply(first:(K - 1), function(k) get(paste("p_", k, "_prime", sep = "")))))
    res_names <- c(res_names, paste("p_", first:(K - 1) - (first - 1), "_prime", sep = ""))
  }
  if (any(intervention_type %in% c("all", "shift_all"))) {
    res <- c(res, p_all)
    res_names <- c(res_names, "p_all")
  }
  names(res) = res_names

  if (!flag)
    return(res)
  else
    return(rep(NA, length(res)))
}



