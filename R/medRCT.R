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
#' @param bootstrap logical. If \code{TRUE}, bootstrap will be conducted.
#' @param boot_args a \code{list} of bootstrapping arguments. \code{R} is the number of bootstrap replicates.
#' \code{stype} indicates what the second argument of \code{statistics} in the \code{boot} function represents
#' @param ... other arguments passed to the \code{boot} function in the \code{boot} package.
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
                   boot_args = list(R = 100, stype = "i"),
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

  # set boot ci type for now
  ci.type = "perc"

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
      bt <- boot::boot.ci(boot.out, index = i, type = ci.type)
      cilow <- c(cilow, bt$percent[4])
      ciupp <- c(ciupp, bt$percent[5])
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


#' Causal mediation analysis for estimating the interventional effect mapped to a target trial
#'
#' @param dat A data.frame with the data for analysis
#' @param ind A vector of indices that define the sample from dat on which to conduct the analysis. Defaults to all rows of the data.
#' This facilitates the use of this function within the boot() function from the boot package
#' @param first index of first mediator of interest after combining intermediate confounders with mediators
#' @param K the number of mediators
#' @param fam_type family type for mediators
#' @param interactions_XC a character string specifying the exposure-confounder or confounder-confounder interaction terms
#' to include in all regression models in the procedure. Defaults to include all two-way exposure-confounder interactions
#' and no confounder-confounder interactions.
#' @param intervention_type a character string indicating the type of the interventional effect to be estimated.
#' Can be 'all', 'shift_all', 'shift_k', 'shift_k_order'. Default is 'all'.
#' @param mcsim the number of Monte Carlo simulations to conduct
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

  # Set flag to capure bootstrap samples to reject
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
  # Joint of M1 to MK under X=0 and X=1 ...

  for (k in 1:K) {

    formula_str <- if (k == 1) {
      paste0("M", k, "~ X +", interactions_XC)
    } else {
      paste0("M", k,
             "~ (X +",
             paste0(paste0("M", 1:(k - 1)), collapse = "+"),
             ")^2 +",
             interactions_XC)
    }

    fit <- glm(as.formula(formula_str), data = data, family = fam_type[[k]])

    if ((!fit$converged) | any(is.na(fit$coefficients)))
      flag <- TRUE

    for(a in exposure_level) {

      dat2 = set_exposure(data = dat2, column_name = "X", exp_val = a)

      if (k != 1) {
        l = 1:(k - 1)
        dat2[, paste0("M", l) := mget(med_outcome_name(a = a,
                                                       l = l,
                                                       K = K))]
      }


      if(fam_type[[k]]$family == "binomial"){
        dat2[, paste0("m", k, "_", a, "_",
                     paste0(c(rep(paste0(a), (k - 1)),
                             rep("m", K - (k - 1))),
                           collapse = "")) :=
               stats::rbinom(n, 1, predict(fit, newdata = dat2, type = "response"))]
      } else if (fam_type[[k]]$family == "gaussian") {
        dat2[, paste0("m", k, "_", a, "_",
                     paste0(c(rep(paste0(a), (k - 1)),
                             rep("m", K - (k - 1))),
                           collapse = "")) :=
               stats::rnorm(n, mean = predict(fit,newdata=dat2,type="response"),
                     sd = sqrt(sum(fit$residuals^2)/stats::df.residual(fit)))]
      }
    }
  }

  # Estimating the target quantities
  # Marginals under X=0
  for (k in first:K) {
    fit <- glm(as.formula(paste0("M", k, "~ X +", interactions_XC)),
               data = data, family = fam_type[[k]])

    if ((!fit$converged) | any(is.na(fit$coefficients)))
      flag <- TRUE

    a <- 0

    # covert X to the correct class
    dat2 = set_exposure(data = dat2, column_name = "X", exp_val = a)


    if(fam_type[[k]]$family == "binomial"){
      dat2[, paste0("m", k, "_", a, "_", paste0(rep("m", K), collapse = "")) :=
             stats::rbinom(n, 1, predict(fit, newdata = dat2, type = "response"))]
    } else if (fam_type[[k]]$family == "gaussian") {
      dat2[, paste0("m", k, "_", a, "_", paste0(rep("m", K), collapse = "")) :=
             stats::rnorm(n, mean = predict(fit,newdata=dat2,type="response"),
                   sd = sqrt(sum(fit$residuals^2)/stats::df.residual(fit)))]
    }
  }


  # For p_first,..., p_K
  # Joint of others under X!=0
  if (any(intervention_type %in% c("all", "shift_k"))) {
    for (MM in first:K) {
      for (k in setdiff(first:K, MM)) {
        # without intermediate confounders
        if (first == 1) {
          if (MM == 1 & k == setdiff(first:K, MM)[1]) {
            fit <- glm(as.formula(paste0("M", k, "~X+", interactions_XC)),
                       data = data,
                       family = fam_type[[k]])
          } else if (!(MM != 1 & k == setdiff(first:K, MM)[1])) {
            fit <- glm(as.formula(
              paste0("M", k, "~(X+", paste0(paste0("M", setdiff(1:(k - 1), MM)), collapse = "+"),
                ")^2+", interactions_XC)),
              data = data,
              family = fam_type[[k]])
          }
          if ((!fit$converged) | any(is.na(fit$coefficients)))
            flag <- TRUE

          for(a in lnzero){

            dat2 = set_exposure(data = dat2, column_name = "X", exp_val = a)

            if (k != setdiff(first:K, MM)[1]) {
              l = setdiff(1:(k - 1), MM)
              dat2[, paste0("M", l) := mget(med_outcome_name(a = a,
                                                             l = l,
                                                             K = K))]
            }

            if(fam_type[[k]]$family == "binomial"){
              dat2[, paste0("m", k, "_", a, "_", paste0(c(
                rep(paste0(a), min(k - 1, MM - 1)),
                "m",
                rep(paste0(a), max(k - 1 - MM, 0)),
                rep("m", K - 1 - min(k - 1, MM - 1) - max(k - 1 - MM, 0))
              ), collapse = "")) :=
                stats::rbinom(n, 1, predict(fit, newdata = dat2, type = "response"))]
            } else if (fam_type[[k]]$family == "gaussian") {
              dat2[, paste0("m", k, "_", a, "_", paste0(c(
                rep(paste0(a), min(k - 1, MM - 1)),
                "m",
                rep(paste0(a), max(k - 1 - MM, 0)),
                rep("m", K - 1 - min(k - 1, MM - 1) - max(k - 1 - MM, 0))
              ), collapse = "")) :=
                stats::rnorm(n, mean = predict(fit,newdata=dat2,type="response"),
                      sd = sqrt(sum(fit$residuals^2)/stats::df.residual(fit)))]
            }
          }
        } else {
          # with intermediate confounders
          fit <- glm(as.formula(
            paste0("M", k, "~(X+",
              paste0(paste0("M", setdiff(1:(k - 1), MM)), collapse = "+"),
              ")^2+", interactions_XC)),
            data = data,
            family = fam_type[[k]])

          if ((!fit$converged) | any(is.na(fit$coefficients)))
            flag <- TRUE

          for(a in lnzero){

            dat2 = set_exposure(data = dat2, column_name = "X", exp_val = a)

            if(k!=1){
              l = setdiff(1:(k - 1), MM)
              dat2[, paste0("M", l) := mget(med_outcome_name(a = a,
                                                             l = l,
                                                             K = K))]
            }

            if(fam_type[[k]]$family == "binomial"){
              dat2[, paste0("m", k, "_", a, "_", paste0(c(
                rep(paste0(a), min(k - 1, MM - 1)),
                "m",
                rep(paste0(a), max(k - 1 - MM, 0)),
                rep("m", K - 1 - min(k - 1, MM - 1) - max(k - 1 - MM, 0))
              ), collapse = "")) :=
                stats::rbinom(n, 1, predict(fit, newdata = dat2, type = "response"))]
            } else if (fam_type[[k]]$family == "gaussian") {
              dat2[, paste0("m", k, "_", a, "_", paste0(c(
                rep(paste0(a), min(k - 1, MM - 1)),
                "m",
                rep(paste0(a), max(k - 1 - MM, 0)),
                rep("m", K - 1 - min(k - 1, MM - 1) - max(k - 1 - MM, 0))
              ), collapse = "")) :=
                stats::rnorm(n, mean = predict(fit,newdata=dat2,type="response"),
                      sd = sqrt(sum(fit$residuals^2)/stats::df.residual(fit)))]
            }
          }
        }
      }
    }
  }

  # For p_first_prime,...., p_K_prime
  # Conditionals under X=1
  if (any(intervention_type %in% c("all", "shift_k_order"))) {
    for (MM in first:(K - 1)) {
      for (k in (MM + 1):K) {
        fit <- glm(as.formula(paste0(
          "M", k, "~(X+", paste0(paste0("M", 1:(k - 1)), collapse = "+"),
          ")^2+",
          interactions_XC)),
          data = data,
          family = fam_type[[k]])

        if ((!fit$converged) | any(is.na(fit$coefficients)))
          flag <- TRUE

        for(a in lnzero){

          dat2 = set_exposure(data = dat2, column_name = "X", exp_val = a)

          if (MM != 1) {
            l = 1:(MM - 1)
            dat2[, paste0("M", l) := mget(med_outcome_name(a = a,
                                                           l = l,
                                                           K = K))]
          }

          dat2[, paste0("M", MM) := get(paste0("m", MM, "_", 0, "_",
                                               paste0(rep("m", K), collapse = "")))]

          if (k > (MM + 1)) {
            for (l in (MM + 1):(k - 1))
              dat2[, paste0("M", l) := get(paste0("m", l, "_", a, "_", paste0(c(
                rep(paste0(a), MM - 1),
                0,
                rep(paste0(a), max(l - 1 - MM, 0)),
                rep("m", K - MM - max(l - 1 - MM, 0))
              ), collapse = "")))]
          }

          if(fam_type[[k]]$family == "binomial"){
            dat2[, paste0("m", k, "_", a, "_", paste0(c(
              rep(paste0(a), MM - 1),
              0,
              rep(paste0(a), max(k - 1 - MM, 0)),
              rep("m", K - MM - max(k - 1 - MM, 0))), collapse = "")) :=
                stats::rbinom(n, 1, predict(fit, newdata = dat2, type = "response"))]
          } else if (fam_type[[k]]$family == "gaussian") {
            dat2[, paste0("m", k, "_", a, "_", paste0(c(
              rep(paste0(a), MM - 1),
              0,
              rep(paste0(a), max(k - 1 - MM, 0)),
              rep("m", K - MM - max(k - 1 - MM, 0))), collapse = "")) :=
                stats::rnorm(n, mean = predict(fit,newdata=dat2,type="response"),
                      sd = sqrt(sum(fit$residuals^2)/stats::df.residual(fit)))]
          }
        }
      }
    }
  }



  # For p_all
  # Joint of main ones under X=0
  if (any(intervention_type %in% c("all", "shift_all"))) {
    for (k in (first + 1):K) {
      fit <- glm(as.formula(
        paste0("M", k, "~(X+",
          paste0(paste0("M", first:(k - 1)), collapse = "+"),
          ")^2+",
          interactions_XC)),
        data = data,
        family = fam_type[[k]])

      if ((!fit$converged) | any(is.na(fit$coefficients)))
        flag <- TRUE

      a <- 0

      dat2 = set_exposure(data = dat2, column_name = "X", exp_val = a)

      for (l in first:(k - 1)) {
        dat2[, paste0("M", l) := get(paste0("m", l, "_", a, "_", paste0(c(
          rep("m", first - 1), rep(paste0(a), (l - first)), rep("m", K - l + 1)
        ), collapse = "")))]
      }

      if(fam_type[[k]]$family == "binomial"){
        dat2[, paste0("m", k, "_", a, "_", paste0(c(
          rep("m", first - 1), rep(paste0(a), k - first), rep("m", K + 1 - k)
        ), collapse = "")) := stats::rbinom(n, 1, predict(fit, newdata = dat2, type = "response"))]
      } else if (fam_type[[k]]$family == "gaussian") {
        dat2[, paste0("m", k, "_", a, "_", paste0(c(
          rep("m", first - 1), rep(paste0(a), k - first), rep("m", K + 1 - k)
        ), collapse = "")) := stats::rnorm(n, mean = predict(fit,newdata=dat2,type="response"),
                                    sd = sqrt(sum(fit$residuals^2)/stats::df.residual(fit)))]
      }
    }
  }


  # outcome
  # Y
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
    } else {
      p_trt <- mean(y1)
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
      for (k in first:K) {
        dat2[, paste0("M", k) := get(paste0("m", k, "_", 0, "_", paste0(c(
          rep("m", first - 1), rep(paste0(0), k - first), rep("m", K + 1 - k)
        ), collapse = "")))]
      }

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
        dat2[, paste0("M", MM) := get(paste0("m", MM, "_", 0, "_", paste0(paste0(
          rep("m", K)), collapse = "")))]

        for (k in setdiff(first:K, MM)) {
          dat2[, paste0("M", k) := get(paste0("m", k, "_", a, "_", paste0(c(
            rep(paste0(a), min(k - 1, MM - 1)),
            "m",
            rep(paste0(a), max(k - 1 - MM, 0)),
            rep("m", K - 1 - min(k - 1, MM - 1) - max(k - 1 - MM, 0))
          ), collapse = "")))]
        }

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
                                             paste0(paste0(rep("m", K)), collapse = "")))]

        if ((MM + 1) <= K) {
          for (k in (MM + 1):K) {
            dat2[, paste0("M", k) := get(paste0("m", k, "_", a, "_", paste0(c(
              rep(paste0(a), MM - 1),
              0,
              rep(paste0(a), max(k - 1 - MM, 0)),
              rep("m", K - MM - max(k - 1 - MM, 0))
            ), collapse = "")))]
          }
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
      res <- c(res, unlist(sapply(lnzero, function(a) {
        sapply(first:K, function(k) {
          get(paste0("IIE_", k, "_", a))
        })
      })))

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
      res <- c(res, sapply(first:K, function(k)
        get(paste0("IIE_", k))))

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
      res <- c(res, unlist(sapply(lnzero, function(a) {
        sapply(first:(K-1), function(k) {
          get(paste0("IIE_", k, "_prime_", a))
        })
      })))

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

      res <- c(res, sapply(first:(K - 1), function(k)
        get(paste0("IIE_", k, "_prime"))))

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
      res <- c(res, sapply(lnzero, function(a) {
        get(paste0("IIE_all_", a))
      }))

      res_names <- c(res_names, paste0("IIE_all_", lnzero,
                     " (p_trt_", lnzero, " - p_all_", lnzero, ")"))

    } else {
      res <- c(res, IIE_all)
      res_names <- c(res_names, "IIE_all (p_trt - p_all)")
      }
  }

  # p_trt & p_ctr
  if(length(lnzero) > 1){
    res <- c(res, sapply(lnzero, function(a) {
      get(paste0("p_trt_", a))
    }), p_ctr)
    res_names <- c(res_names, paste0("p_trt_", lnzero), "p_ctr")
  } else {
    res <- c(res, p_trt, p_ctr)
    res_names <- c(res_names, "p_trt", "p_ctr")
  }

  # p_k
  if (any(intervention_type %in% c("all", "shift_k"))) {
    if(length(lnzero) > 1){
      res <- c(res, unlist(sapply(lnzero, function(a) {
        sapply(first:K, function(k) {
          get(paste0("p_", k, "_", a))
        })
      })))

      res_names <- c(res_names,
                     paste0("p_",
                            rep(first:K - (first - 1), length(lnzero)),
                            "_",
                            rep(lnzero, each = length(first:K))))
    } else {
      res <- c(res, sapply(first:K, function(k)
        get(paste0("p_", k))))
      res_names <- c(res_names, paste0("p_", first:K - (first - 1)))
    }
  }

  # p_k_prime
  if (any(intervention_type %in% c("all", "shift_k_order"))) {
    if(length(lnzero) > 1){
      res <- c(res, unlist(sapply(lnzero, function(a) {
        sapply(first:(K-1), function(k) {
          get(paste0("p_", k, "_prime_", a))
        })
      })))

      res_names <- c(res_names,
                     paste0("p_",
                            rep(first:(K-1) - (first - 1), length(lnzero)),
                            "_prime_",
                            rep(lnzero, each = length(first:(K-1)))))

    } else {
      res <- c(res, sapply(first:(K - 1), function(k)
        get(paste0("p_", k, "_prime"))))
      res_names <- c(res_names, paste0("p_", first:(K - 1) - (first - 1), "_prime"))
    }
  }

  # p_all
  if (any(intervention_type %in% c("all", "shift_all"))) {
    if(length(lnzero) > 1){
      res <- c(res, sapply(lnzero, function(a) {
        get(paste0("p_all_", a))
      }))

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



