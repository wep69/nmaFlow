#' nmaFlow package
#'
#' Integrated workflows for Network Meta-Analysis.
#'
#' @keywords internal
"_PACKAGE"

# Global variables for ggplot2 aesthetics (avoids R CMD check NOTEs)
utils::globalVariables(c(
  "x", "y", "xend", "yend", "degree", "studies", "treatment", "value",
  ".nma_ymin", ".nma_ymax"
))

.onLoad <- function(libname, pkgname) {
  op <- options()
  defaults <- list(
    nmaFlow.digits = 3,
    nmaFlow.seed = 260819,
    nmaFlow.auto_repair = FALSE
  )
  to_set <- defaults[!names(defaults) %in% names(op)]
  if (length(to_set)) options(to_set)
  invisible()
}
