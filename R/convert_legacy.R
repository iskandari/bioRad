#' Convert legacy bioRad objects
#'
#' Convert legacy bioRad objects (`vp`, `vpts`, `ppi`) and make them compatible with
#' the current bioRad version. Conversion includes renaming `HGHT` to `height`
#' and converting legacy `ppi$data` to `terra::SpatRaster`. Spatial conversion
#' preserves layer names, values, grid and CRS; the original object is unchanged.
#' Legacy `sp`/`raster` inputs warn and require the corresponding optional packages.
#'
#' @param x A `vp`, `vpts` or `ppi` object.
#' @returns An updated object of the same class as the input.
#' @family conversion functions
#' @export
#' @examples
#' # Convert a vp object
#' vp <- convert_legacy(example_vp)
#'
#' # Convert a vpts object
#' vpts <- convert_legacy(example_vpts)
convert_legacy <- function(x) {
  UseMethod("convert_legacy", x)
}

#' @rdname convert_legacy
#' @export
convert_legacy.vp <- function(x) {
  assertthat::assert_that(inherits(x, "vp"))
  names(x$data) <- sub("HGHT", "height", names(x$data))
  x
}

#' @rdname convert_legacy
#' @export
convert_legacy.vpts <- function(x) {
  assertthat::assert_that(inherits(x, "vpts"))
  names(x) <- sub("heights", "height", names(x))
  names(x) <- sub("dates", "datetime", names(x))
  x
}

#' @rdname convert_legacy
#' @export
convert_legacy.ppi <- function(x) {
  x$data <- .spatial_raster(x$data)
  x
}
