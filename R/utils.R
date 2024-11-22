


#' Summarise Results from Mediation Analysis
#'
#' Provides a summary of the results obtained from a mediation analysis conducted
#' using the \code{medRCT} function.
#'
#' @param object An object of class \code{medRCT}.
#' @param ... Additional arguments to customise the summary output.
#' @method summary medRCT
#' @export
summary.medRCT <- function(object, ...){
  if (object$bootstrap == TRUE) {
    index <- grep("^p", names(object$est))[1]
    out1 = cbind(object$est, object$se, object$cilow, object$ciupp, object$pval)[1:index-1,]
    out1 = matrix(out1, ncol = 5)
    row.names(out1) = names(object$est)[1:index-1]
    colnames(out1) <- c("Estimate", "Std. Error", "CI Lower",
                        "CI Upper", "p-value")
    cat("\n")
    cat("Estimated interventional effect: \n\n")
    stats::printCoefmat(out1, signif.stars = F, tst.ind=NULL, digits = 3)
    out2 = cbind(object$est, object$se, object$cilow, object$ciupp, object$pval)[index:length(object$est),]
    colnames(out2) <- c("Estimate", "Std. Error", "CI Lower",
                        "CI Upper", "p-value")
    cat("\n")
    cat("Estimated expected outcome in each trial arm:\n\n")
    stats::printCoefmat(out2, signif.stars = F, tst.ind=NULL, digits = 3)
    cat("\n")
    cat("Sample Size:", object$sample.size,"\n")
    cat("\n")
    cat("Simulations:", object$mcsim,"\n\n")

    # Return as a list for easy extraction of coefficients
    invisible(list(
      IIE = out1,
      expected_outcome = out2,
      sample_size = object$sample.size,
      n.sim = object$mcsim
    ))

  } else {
    object$est
  }
}





#' Determine Appropriate Family Types for GLM
#'
#' Identifies the appropriate family type (\code{"binomial()"} or \code{"gaussian()"}) for a set of variables.
#'
#' @param data A \code{data.frame} containing the variables to be analyzed.
#' @param variable_names A \code{character} vector specifying the names of the variables to evaluate.
#' Each variable name should correspond to a column in the \code{data.frame}.
#' @param unique_threshold An \code{integer} value specifying the minimum number of unique values required for a variable
#' to be classified as continuous. Defaults to \code{10}.
#'
#' @return A \code{list} where each element corresponds to a variable in \code{variable_names}, and the value indicates
#' the family type: either \code{"binomial()"} for binary variables or \code{"gaussian()"} for continuous variables.
#'
#' @importFrom stats binomial gaussian
#' @export
family_type <- function(data, variable_names, unique_threshold = 10) {
  result <- lapply(variable_names, function(var) {
    unique_vals <- unique(data[[var]])
    if (length(unique_vals) == 2) {
      return(stats::binomial())
    } else if (length(unique_vals) > unique_threshold) {
      return(stats::gaussian())
    } else {
      stop("Error: The variable must be either continuous or binary.")
    }
  })
  return(result)
}



#' Set Exposure Column to a Specific Value and Convert to Factor
#'
#' @param data A `data.table` object.
#' @param column_name A character string specifying the name of the column to modify.
#' @param exp_val The exposure value to assign to the column.
#'
#' @importFrom collapse qF
set_exposure = function(data, column_name, exp_val) {
  # Convert column to numeric
  data[, (column_name) := as.numeric(get(column_name))]
  # Assign the new value
  data[, (column_name) := exp_val]
  # Convert column to factor
  data[, (column_name) := collapse::qF(get(column_name))]
  data
}


#' Counterfactual Prediction and Random Draw
#'
#' @param fit A fitted glm model object
#' @param data A `data.table` containing the data to which the counterfactual predictions
#'             will be applied.
#' @param var_name A character string specifying the name of the variable to store
#'                 the random draws.
#' @param n An integer specifying the number of observations in the dataset.
#' @param family A character string specifying the distribution family to use for
#'               generating random draws. Must be either `"binomial"` or `"gaussian"`.
cf_predict = function(fit, data, var_name, n, family){
  # get predictions
  predictions <- predict(fit, newdata = data, type = "response")

  # counterfactual prediction and random draw
  data[, (var_name) := switch(
    family,
    "binomial" = stats::rbinom(n, 1, predictions),
    "gaussian" = stats::rnorm(
      n,
      mean = predictions,
      sd = sqrt(sum(fit$residuals^2) / stats::df.residual(fit))
    )
  )]

  data
}


med_outcome_name = function(l, a, K) {
  paste0("m", l, "_", a, "_", paste0(
    strrep(a, (l - 1)), strrep("m", K - (l - 1))
  ))
}

med_outcome_all = function(l, first, a, K){
  paste0("m", l, "_", a, "_", paste0(
    strrep("m", first - 1), strrep(a, (l - first)), strrep("m", K - l + 1)))
}









