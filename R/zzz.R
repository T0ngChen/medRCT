.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste0(
    gsub("\n", " ", paste0("medRCT: ", utils::packageDescription("medRCT", fields = "Title"))), "\n",
    "Note: When setting intervention_type = 'shift_k_order', the order of the mediators as specified in the 'mediators' argument is important."
  ))
}


