.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste0(
    "medRCT", ": ", utils::packageDescription("medRCT")$Title, "\n",
    "IMPORTANT: when estimating the effect type 'shift_k_order', the order of mediators specified in the 'mediators' argument is important."
  ))
}
