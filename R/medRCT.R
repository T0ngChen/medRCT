


#' Causal mediation analysis for estimating the interventional effect
#'
#' 'medRCT' is used to estimating the interventional effects that emulate a target trial with any number of multiple
#' mediators, including some not of primary interest but that are intermediate confounders
#'
#' @param dat a data.frame with the data for analysis
#' @param exposure a character string with name of the exposure in the data. The exposure variable can only be binary.
#' @param outcome a character string with name of the outcome in the data. The outcome variable con only be binary.
#' @param mediators a character vector with the names of the mediators in data. The mediators can only be binary variables.
#' @param int.confounders a character vector with the names of the intermediate confounders in data. The intermediate
#' confounders can only be binary variables.
#' @param confounders a character vector with names of the confounders in data, which must be of the required class
#' (e.g. factor if appropriate)
#' @param Xconfsint a character string specifying the exposure-confounder or confounder-confounder interaction terms
#' to include in all regression models in the procedure. Defaults to include all two-way exposure-confounder interactions
#' and no confounder-confounder interactions.
#' @param effect_type a character string indicating the type of the interventional effect to be estimated.
#' Can be 'all', 'E_all', 'E_k', 'E_prime'. Default is 'all'.
#' @param mcsim the number of Monte Carlo simulations to conduct
#' @param boot_args a \code{list} of bootstrapping arguments. \code{R} is the number of bootstrap replicates.
#' \code{stype} indicates what the second argument of \code{statistics} in the \code{boot} function represents
#' @param ... other arguments passed to the \code{boot} package.
#'
#' @export
medRCT <- function(dat, exposure, outcome, mediators, int.confounders, confounders,
                   Xconfsint = paste(paste(rep("X", length(confounders)), confounders, sep ="*"), collapse = "+"),
                   effect_type = c("all", "E_all", "E_k", "E_prime"), mcsim,
                   boot_args = list(R = 100, stype = "i"), ...) {
  effect_type = match.arg(effect_type)
  ci.type = "perc"
  mediators = c(int.confounders, mediators)
  first = length(int.confounders) + 1
  K <- length(mediators)
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
  if (is.null(Xconfsint))
    Xconfsint <- paste(paste(rep("X", length(confounders)), confounders, sep ="*"),
                       collapse = "+")

  boot.out <- boot::boot(
    data = dat,
    statistic = medRCT.fun,
    first = first,
    K = K,
    mcsim = mcsim,
    Xconfsint = Xconfsint,
    effect_type = effect_type,
    stype = boot_args$stype,
    R = boot_args$R,
    ...
  )

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
    mcsim = mcsim
  )

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
#' @param Xconfsint a character string specifying the exposure-confounder or confounder-confounder interaction terms
#' to include in all regression models in the procedure. Defaults to include all two-way exposure-confounder interactions
#' and no confounder-confounder interactions.
#' @param effect_type a character string indicating the type of the interventional effect to be estimated.
#' Can be 'all', 'E_all', 'E_k', 'E_prime'. Default is 'all'.
#' @param mcsim the number of Monte Carlo simulations to conduct
#'
#' @importFrom stats as.formula binomial glm predict rbinom
medRCT.fun <- function(dat,
                       ind = 1:nrow(dat),
                       first = first,
                       K = K,
                       Xconfsint = Xconfsint,
                       effect_type = effect_type,
                       mcsim) {
  # Take boostrap sample
  data <- dat[ind, ]

  # Set flag to capure bootstrap samples to reject
  flag <- FALSE

  # Replicate dataset for simulations
  dat2 <- data
  dat2[, 1:(2 + K)] <- NA_integer_
  dat2 <- zoo::coredata(dat2)[rep(seq(nrow(dat2)), mcsim), ]
  n <- nrow(dat2)

  # ESTIMATE DISTRIBUTIONS
  # Joint of M1 to MK under X=0 and X=1

  for (k in 1:K) {
    if (k == 1)
      fit <- glm(as.formula(paste("M", k, "~X+", Xconfsint, sep = "")),
                 data = data, family = binomial)
    else
      fit <- glm(as.formula(paste("M", k, "~(X+",
                                  paste(paste("M", 1:(k - 1), sep = ""), collapse = "+"), ")^2+",
                                  Xconfsint, sep = "")),
                 data = data, family = binomial)
    if ((!fit$converged) | any(is.na(fit$coefficients)))
      flag <- TRUE

    for (a in c(0, 1)) {
      dat2$X <- a

      if (k != 1) {
        for (l in 1:(k - 1))
          dat2[, paste("M", l, sep = "")] <- get(
            paste("m", l, "_", a, "_",
                  paste(c(rep(paste(a), (l - 1)), rep("m", K - (l - 1))), collapse = ""), sep = ""))
      }

      assign(paste("m", k, "_", a, "_",
                   paste(c(rep(paste(a), (k - 1)), rep("m", K - (k - 1))), collapse = ""),
                   sep = ""),
             rbinom(n, 1, predict(fit, newdata = dat2, type = "response")))
    }
  }

  # Estimating the target quantities
  # Marginals under X=0
  for (k in first:K) {
    fit <- glm(as.formula(paste("M", k, "~X+", Xconfsint, sep = "")), data =
                 data, family = binomial)

    if ((!fit$converged) | any(is.na(fit$coefficients)))
      flag <- TRUE

    a <- 0
    dat2$X <- a
    assign(paste("m", k, "_", a, "_", paste(rep("m", K), collapse = ""), sep = ""),
           rbinom(n, 1, predict(fit, newdata = dat2, type = "response")))
  }


  # For p_first,..., p_K
  # Joint of others under X=1
  if (effect_type == "all" | effect_type == "E_k") {
    for (MM in first:K) {
      for (k in setdiff(first:K, MM)) {
        fit <- glm(as.formula(paste("M", k, "~(X+",
                                    paste(paste("M", setdiff(1:(k - 1), MM), sep = ""),
                                          collapse = "+"),")^2+", Xconfsint, sep = "")),
                   data = data, family = binomial)

        if ((!fit$converged) | any(is.na(fit$coefficients)))
          flag <- TRUE

        a <- 1
        dat2$X <- a

        if (k != 1) {
          for (l in setdiff(1:(k - 1), MM))
            dat2[, paste("M", l, sep = "")] <- get(
              paste("m", l, "_", a, "_", paste(c(rep(paste(a), (l - 1)), rep("m", K - (l - 1))),
                                               collapse = ""), sep = ""))
        }
        assign(paste("m", k, "_", a, "_",
                     paste(c(rep(paste(a), min(k - 1, MM - 1)), "m", rep(paste(a), max(k - 1 - MM, 0)),
                             rep("m", K - 1 - min(k - 1, MM - 1) - max(k - 1 - MM, 0))), collapse = ""), sep = ""),
               rbinom(n, 1, predict(fit, newdata = dat2, type = "response")))
      }
    }
  }


  # For p_first_prime,...., p_K_prime
  # Conditionals under X=1
  if (effect_type == "all" | effect_type == "E_prime") {
    for (MM in first:(K - 1)) {
      for (k in (MM + 1):K) {
        fit <- glm(as.formula(paste("M", k, "~(X+",
                                    paste(paste("M", 1:(k - 1), sep = ""), collapse = "+"),
                                    ")^2+", Xconfsint, sep = "")),
                   data = data, family = binomial)

        if ((!fit$converged) | any(is.na(fit$coefficients)))
          flag <- TRUE

        a <- 1
        dat2$X <- a

        if (MM != 1) {
          for (l in 1:(MM - 1))
            dat2[, paste("M", l, sep = "")] <- get(
              paste("m", l, "_", a, "_", paste(c(rep(paste(a), (l - 1)), rep("m", K - (l - 1))),
                                               collapse = ""), sep = ""))
        }
        dat2[, paste("M", MM, sep = "")] <- get(paste("m", MM, "_", 0, "_",
                                                      paste(rep("m", K), collapse = ""), sep = ""))
        if (k > (MM + 1)) {
          for (l in (MM + 1):(k - 1))
            dat2[, paste("M", l, sep = "")] <- get(
              paste("m", l, "_", a, "_", paste(c(rep(paste(a), MM - 1), 0, rep(paste(a), max(l - 1 - MM, 0)),
                                                 rep("m", K - MM - max(l - 1 - MM, 0))), collapse = ""),
                    sep = ""))
        }

        assign(paste("m", k, "_", a, "_", paste(c(rep(paste(a), MM - 1), 0, rep(paste(a), max(k - 1 - MM, 0)),
                                                  rep("m", K - MM - max(k - 1 - MM, 0))), collapse = ""), sep = ""),
               rbinom(n, 1, predict(fit, newdata = dat2, type = "response")))
      }
    }
  }



  # For p_all
  # Joint of main ones under X=0
  if (effect_type == "all" | effect_type == "E_all") {
    for (k in (first + 1):K) {
      fit <- glm(as.formula(paste("M", k, "~(X+",
                                  paste(paste("M", first:(k - 1), sep = ""), collapse = "+"),
                                  ")^2+", Xconfsint, sep = "")),
                 data = data, family = binomial)

      if ((!fit$converged) | any(is.na(fit$coefficients)))
        flag <- TRUE

      a <- 0
      dat2$X <- a

      for (l in first:(k - 1)) {
        dat2[, paste("M", l, sep = "")] <- get(
          paste("m", l, "_", a, "_", paste(c(rep("m", first - 1), rep(paste(a), (l - first)),
                                             rep("m", K - l + 1)), collapse = ""), sep = ""))
      }

      assign(paste("m", k, "_", a, "_",
                   paste(c(rep("m", first - 1), rep(paste(a), k - first), rep("m", K + 1 - k)),
                         collapse = ""), sep = ""),
             rbinom(n, 1, predict(fit, newdata = dat2, type = "response")))
    }
  }


  # outcome
  # Y

  fit <- glm(as.formula(
    paste("Y~(X+", paste(paste("M", 1:K, sep = ""), collapse = "+"), ")^2+", Xconfsint, sep ="")),
    data = data, family = binomial)
  if ((!fit$converged) | any(is.na(fit$coefficients)))
    flag <- TRUE


  # ESTIMATE OUTCOME EXPECTATION IN EACH ARM & ESTIMATE EFFECTS
  # p_ctr

  a <- 0
  dat2$X <- a
  for (k in 1:K) {
    dat2[, paste("M", k, sep = "")] <- get(
      paste("m", k, "_", a, "_", paste(c(rep(paste(a), (k - 1)), rep("m", K - (k - 1))),
                                       collapse = ""), sep = ""))
  }

  y0 <- predict(fit, newdata = dat2, type = "response")

  p_ctr <- mean(y0)


  # p_trt

  a <- 1
  dat2$X <- a
  for (k in 1:K) {
    dat2[, paste("M", k, sep = "")] <- get(
      paste("m", k, "_", a, "_", paste(c(rep(paste(a), (k - 1)), rep("m", K - (k - 1))),
                                       collapse = ""), sep = ""))
  }

  y1 <- predict(fit, newdata = dat2, type = "response")

  p_trt <- mean(y1)


  # p_all
  if (effect_type == "all" | effect_type == "E_all") {
    a <- 1
    dat2$X <- a
    for (k in 1:(first - 1)) {
      dat2[, paste("M", k, sep = "")] <- get(
        paste("m", k, "_", a, "_", paste(c(rep(paste(a), (k - 1)), rep("m", K - (k - 1))),
                                         collapse = ""), sep = ""))
    }

    a <- 0
    for (k in first:K) {
      dat2[, paste("M", k, sep = "")] <- get(
        paste("m", k, "_", a, "_", paste(c(rep("m", first - 1), rep(paste(a), k - first), rep("m", K + 1 - k)),
                                         collapse = ""), sep = ""))
    }

    y1 <- predict(fit, newdata = dat2, type = "response")

    p_all <- mean(y1)
    # Interventional effects
    IIE_all <- p_trt - p_all
  }


  # p_first....p_K
  if (effect_type == "all" | effect_type == "E_k") {
    a <- 1
    dat2$X <- a

    for (k in 1:(first - 1)) {
      dat2[, paste("M", k, sep = "")] <- get(
        paste("m", k, "_", a, "_", paste(c(rep(paste(a), (k - 1)), rep("m", K - (k - 1))),
                                         collapse = ""), sep = ""))
    }

    for (MM in first:K) {
      dat2[, paste("M", MM, sep = "")] <- get(
        paste("m", MM, "_", 0, "_", paste(paste(rep("m", K), sep = ""), collapse = ""), sep = ""))

      for (k in setdiff(first:K, MM)) {
        dat2[, paste("M", k, sep = "")] <- get(
          paste("m", k, "_", a, "_", paste(c(rep(paste(a), min(k - 1, MM - 1)), "m", rep(paste(a), max(k - 1 - MM, 0)),
                                             rep("m", K - 1 - min(k - 1, MM - 1) - max(k - 1 - MM, 0))), collapse = ""), sep = ""))
      }

      y0 <- predict(fit, newdata = dat2, type = "response")

      assign(paste("p_", MM, sep = ""), mean(y0))
    }
    # Interventional effects
    for (k in first:K)
      assign(paste("IIE_", k, sep = ""), p_trt - get(paste("p_", k, sep ="")))
  }


  # p_first_prime....p_Kminus1_prime
  if (effect_type == "all" | effect_type == "E_prime") {
    a <- 1
    dat2$X <- a

    for (MM in first:(K - 1)) {
      if (MM != 1) {
        for (k in 1:(MM - 1)) {
          dat2[, paste("M", k, sep = "")] <- get(
            paste("m", k, "_", a, "_", paste(c(rep(paste(a), (k - 1)), rep("m", K - (k - 1))),
                                             collapse = ""), sep = ""))
        }
      }

      dat2[, paste("M", MM, sep = "")] <- get(
        paste("m", MM, "_", 0, "_", paste(paste(rep("m", K), sep = ""), collapse = ""), sep = ""))

      if ((MM + 1) <= K) {
        for (k in (MM + 1):K) {
          dat2[, paste("M", k, sep = "")] <- get(
            paste("m", k, "_", a, "_", paste(c(rep(paste(a), MM - 1), 0, rep(paste(a), max(k - 1 - MM, 0)),
                                               rep("m", K - MM - max(k - 1 - MM, 0))), collapse = ""), sep = ""))
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

  if (effect_type == "all" | effect_type == "E_k") {
    res <- c(res, unlist(lapply(first:K, function(k) get(paste("IIE_", k, sep = "")))))
    res_names <- c(res_names, paste0("IIE_", first:K - (first - 1),
                                     " (p_trt - p_", first:K - (first - 1), ")"))
  }

  if (effect_type == "all" | effect_type == "E_prime") {
    res <- c(res, unlist(lapply(first:(K - 1), function(k) get(paste("IIE_", k, "_prime", sep = "")))))
    res_names <- c(res_names, paste0("IIE_", first:(K - 1) - (first - 1), "_prime",
                                     " (p_trt - p_", first:(K - 1) - (first - 1), "_prime)"))
  }
  if (effect_type == "all" | effect_type == "E_all") {
    res <- c(res, IIE_all)
    res_names <- c(res_names, "IIE_all (p_trt - p_all)")
  }
  res <- c(res, p_trt, p_ctr)
  res_names <- c(res_names, "p_trt", "p_ctr")
  if (effect_type == "all" | effect_type == "E_k") {
    res <- c(res, unlist(lapply(first:K, function(k) get(paste("p_", k, sep = "")))))
    res_names <- c(res_names, paste("p_", first:K - (first - 1), sep = ""))
  }

  if (effect_type == "all" | effect_type == "E_prime") {
    res <- c(res, unlist(lapply(first:(K - 1), function(k) get(paste("p_", k, "_prime", sep = "")))))
    res_names <- c(res_names, paste("p_", first:(K - 1) - (first - 1), "_prime", sep = ""))
  }
  if (effect_type == "all" | effect_type == "E_all") {
    res <- c(res, p_all)
    res_names <- c(res_names, "p_all")
  }
  names(res) = res_names

  if (!flag)
    return(res)
  else
    return(rep(NA, length(res)))
}



