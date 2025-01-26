#' Simulated Data: the Longitudinal Study of Australian Children (LSAC)
#'
#' This dataset contains simulated data based on the Longitudinal Study of Australian Children (LSAC).
#' The data comprises 5,107 rows and 11 columns.
#'
#' @usage LSACdata
#'
#' @format A \code{data.frame} with 5,107 rows, containing the following variables:
#' \describe{
#'   \item{sep}{A binary exposure variable}
#'   \item{fam_stress}{A binary intermediate confounder}
#'   \item{parent_mh}{A binary mediator}
#'   \item{preschool_att}{A binary mediator}
#'   \item{child_mh}{A binary outcome variable}
#'   \item{child_SDQscore}{A continuous outcome variable}
#'   \item{child_sex}{A binary baseline confounder}
#'   \item{child_atsi}{A binary baseline confounder}
#'   \item{mat_cob}{A binary baseline confounder}
#'   \item{mat_engl}{A binary baseline confounder}
#'   \item{mat_age}{A binary baseline confounder}
#' }
"LSACdata"


