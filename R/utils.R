


#' Summarizing Output from Mediation Analysis
#'
#' @param object a \code{medRCT} object
#' @param ... other arguments
#' @method summary medRCT
#' @export
summary.medRCT <- function(object, ...){
  index <- grep("^p", names(object$est))[1]

  out1 = cbind(object$est, object$se, object$cilow, object$ciupp, object$pval)[1:index-1,]
  colnames(out1) <- c("Estimate", "Std. Error", "CI Lower",
                     "CI Upper", "p-value")
  cat("\n")
  cat("Estimated interventional effect: \n\n")
  stats::printCoefmat(out1, digits = 3, signif.stars = F)
  out2 = cbind(object$est, object$se, object$cilow, object$ciupp, object$pval)[index:length(object$est),]
  colnames(out2) <- c("Estimate", "Std. Error", "CI Lower",
                      "CI Upper", "p-value")
  cat("\n")
  cat("Estimated expected outcome in each trial arm:\n\n")
  stats::printCoefmat(out2, digits = 3, signif.stars = F)
  cat("\n")
  cat("Sample Size:", object$sample.size,"\n")
  cat("\n")
  cat("Simulations:", object$mcsim,"\n\n")
}

