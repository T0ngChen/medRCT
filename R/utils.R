


#' Summarizing Output from Mediation Analysis
#'
#' @param object a \code{medRCT} object
#' @param ... other arguments
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
    return(list(
      IIE = out1,
      expected_outcome = out2,
      sample_size = object$sample.size,
      n.sim = object$mcsim
    ))

  } else {
    object$est
  }
}




#' Determine Family Type for GLM
#'
#' @param data A data frame containing the variables to be analyzed.
#' @param variable_names A character vector of variable names to check. Each variable name should correspond
#' to a column in the data frame.
#' @param unique_threshold An integer value specifying the minimum number of unique values for a variable to
#' be considered continuous. Default is 10.
#'
#' @return A list where each element corresponds to a family type for variables from `variable_names`
#' and the value is either "binomial()" or "gaussian()", indicating the type of the family.
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



set_exposure = function(data, column_name, exp_val) {
  # Convert column to numeric
  data[, (column_name) := as.numeric(get(column_name))]
  # Assign the new value
  data[, (column_name) := exp_val]
  # Convert column to factor
  data[, (column_name) := as.factor(get(column_name))]
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





