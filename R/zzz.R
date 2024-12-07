.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste0(
    gsub("\n", " ", paste0("medRCT: ", utils::packageDescription("medRCT", fields = "Title"))), "\n",
    "When estimating the effect type 'shift_k_order', the order of mediators specified in the 'mediators' argument is important."
  ))
}


