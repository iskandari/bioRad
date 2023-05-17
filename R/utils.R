#' A wrapper for spTransform
#'
#' This function calls spTransform, setting the sp_evolution_status option to 2L.
#'
#' @param ... Parameters to be passed to spTransform.
#'
#' @return The result of calling spTransform with the specified parameters and
#'   sp_evolution_status set to 2L.
#'
#' @export
spTransform2 <- function(...) {
  withr::with_options(list(sp_evolution_status = 2L), sp::spTransform(...))
}