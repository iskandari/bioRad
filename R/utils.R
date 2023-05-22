#' A wrapper for spTransform
#'
#' This function calls spTransform, setting the sp_evolution_status option to 2L.
#'
#' @param ... Parameters to be passed to spTransform.
#'
#' @return The result of calling spTransform with the specified parameters and
#'   sp_evolution_status set to 2L.
#'
#' @keywords internal
spTransform2 <- function(x, CRSobj, ...) {
       as(sf::st_transform(sf::st_as_sf(x), sf::st_crs(CRSobj)), "Spatial")
}
