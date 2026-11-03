# Compatibility is confined to these adapters; new calculations use sf/terra.
.spatial_raster <- function(x) {
  if (inherits(x, c("Raster", "SpatialGridDataFrame", "SpatialPixelsDataFrame"))) {
    lifecycle::deprecate_warn("1.0.0", I("Support for legacy spatial objects"),
      details = "Use terra::SpatRaster inputs and convert_legacy() for saved ppi objects.",
      id = "legacy-spatial")
    rlang::check_installed("sp", reason = "to convert legacy spatial objects")
    if (inherits(x, "Raster")) rlang::check_installed("raster", reason = "to convert legacy Raster objects")
    x <- terra::rast(x)
  }
  if (!inherits(x, "SpatRaster")) stop("Expected a terra::SpatRaster")
  x
}

.spatial_crs <- function(x) {
  if (inherits(x, "CRS")) {
    lifecycle::deprecate_warn("1.0.0", I("Support for sp::CRS objects"),
      details = "Supply an sf CRS, EPSG code or CRS string instead of sp::CRS.",
      id = "legacy-crs")
  }
  tryCatch(sf::st_crs(x)$wkt, error = function(e) stop("invalid crs: ", as.character(x)))
}

.raster_points <- function(x, crs = terra::crs(x)) {
  points <- sf::st_as_sf(as.data.frame(terra::xyFromCell(x, seq_len(terra::ncell(x)))),
    coords = c("x", "y"), crs = terra::crs(x))
  sf::st_transform(points, crs)
}

.raster_bbox <- function(x) {
  b <- sf::st_bbox(sf::st_transform(sf::st_as_sfc(sf::st_bbox(
    c(xmin = terra::xmin(x), ymin = terra::ymin(x), xmax = terra::xmax(x), ymax = terra::ymax(x)),
    crs = sf::st_crs(terra::crs(x)))), 4326))
  matrix(b, nrow = 2, dimnames = list(c("lon", "lat"), c("min", "max")))
}
