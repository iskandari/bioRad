#' Convert a polar scan into a simple feature object
#'
#' Georeferences the center of  pixels for a scan into a
#' [`sf::sf`] object.
#'
#' @inheritParams beam_height
#' @param scan a scan (sweep) of class scan
#' @param lat Geodetic latitude of the radar in degrees. If missing taken from `scan`.
#' @param lon Geodetic longitude of the radar in degrees. If missing taken from `scan`.
#' @returns An `sf` object with point geometries.
#' @family raster functions
#' @export
#' @details Beam altitude accounts for the curvature of the earth, using [beam_height].
#' Distance from the radar over the earth's surface is calculated using [beam_distance].
#' @examples
#' # load example scan:
#' data(example_scan)
#'
#' # convert to an sf object:
#' scan_to_sf(example_scan)
scan_to_sf <- function(scan, lat, lon, k = 4 / 3, re = 6378, rp = 6357) {
  assertthat::assert_that(is.scan(scan))
  assertthat::assert_that(assertthat::is.number(k))
  assertthat::assert_that(assertthat::is.number(re))
  assertthat::assert_that(assertthat::is.number(rp))
  if (is.null(scan$geo$lat) && missing(lat)) stop("radar latitude cannot be found in scan, specify using 'lat' argument")
  if (is.null(scan$geo$lon) && missing(lon)) stop("radar longitude cannot be found in scan, specify using 'lon' argument")
  if (missing(lat)) lat <- scan$geo$lat
  if (missing(lon)) lon <- scan$geo$lon
  assertthat::assert_that(assertthat::is.number(lat))
  assertthat::assert_that(assertthat::is.number(lon))

  crs <- paste("+proj=aeqd +lat_0=", lat,
    " +lon_0=", lon,
    " +units=m",
    sep = ""
  )

  rscale <- scan$geo$rscale
  ascale <- scan$geo$ascale
  elev <- scan$geo$elangle
  rstart <- ifelse(is.null(scan$geo$rstart), 0, scan$geo$rstart)
  astart <- ifelse(is.null(scan$geo$astart), 0, scan$geo$astart)


  data <- data.frame(
    azim = astart + c(t(matrix(rep(seq(0, dim(scan)[3] - 1) * ascale + ascale / 2, dim(scan)[2]), nrow = dim(scan)[3]))),
    range = rstart + rep(seq(1, dim(scan)[2]) * rscale, dim(scan)[3]) - rscale / 2
  )
  data$distance <- beam_distance(range = data$range, elev = elev, k = k, lat = lat, re = re, rp = rp)
  data$height <- scan$geo$height + beam_height(data$range, elev, k = k, lat = lat, re = re, rp = rp)
  data <- cbind(data, as.data.frame(sapply(scan$params, c)))
  coords <- data.frame(
    x = data$distance * cos(pi / 2 - data$azim * pi / 180),
    y = data$distance * sin(pi / 2 - data$azim * pi / 180)
  )
  sf::st_as_sf(cbind(data, coords), coords = c("x", "y"), crs = crs)
}

#' Convert a polar scan into a spatial object
#'
#' `r lifecycle::badge("deprecated")`
#'
#' Retained for compatibility during the spatial migration. Use the corresponding
#' `sf` or `terra` conversion for new code. The legacy package must be installed.
#'
#' @inheritParams scan_to_sf
#' @returns A `SpatialPointsDataFrame`.
#' @family deprecated spatial functions
#' @export
scan_to_spatial <- function(scan, lat, lon, k = 4 / 3, re = 6378, rp = 6357) {
  lifecycle::deprecate_warn("1.0.0", "scan_to_spatial()", "scan_to_sf()")
  call <- match.call()
  call[[1]] <- quote(scan_to_sf)
  rlang::check_installed("sp", reason = "to return a legacy Spatial object")
  sf::as_Spatial(eval(call, parent.frame()))
}

#' Convert a polar scan into a SpatRaster
#'
#' Converts an object of class `scan` into a [`terra::SpatRaster`].
#'
#' @inheritParams scan_to_sf
#' @param nx number of raster pixels in the x (longitude) dimension
#' @param ny number of raster pixels in the y (latitude) dimension
#' @param xlim x (longitude) range
#' @param ylim y (latitude) range
#' @param param scan parameters to include. If `NA` include all scan parameters. Reducing the number
#' of scan parameters speeds up evaluation.
#' @param crs A CRS string, EPSG code or `sf::crs` object. Legacy `sp::CRS`
#' objects are deprecated.
#' When 'NA' (default), an azimuthal equidistant projection with origin at the radar location is used.
#' To use a WSG84 (lat,lon) projection, use crs="+proj=longlat +datum=WGS84"
#' @param res numeric vector of length 1 or 2 to set the resolution of the raster (see [res][terra::res]).
#' If this argument is used, arguments `nx` and `ny` are ignored. Unit is identical to `xlim` and `ylim`.
#' @param raster (optional) `terra::SpatRaster` (legacy `RasterLayer` inputs are deprecated) with a CRS. When specified
#' this raster topology is used for the output, and nx, ny, res arguments are ignored.
#' @returns A `SpatRaster`.
#' @details Samples the polar scan at output cell centres, using the range gate
#' and azimuth bin containing each centre. Cells outside the scan are `NA`.
#' @family raster functions
#' @export
#' @examples
#' \donttest{
#' # default projects full extent on 100x100 pixel raster:
#' scan_to_spatraster(example_scan)
#'
#' # crop the scan and project at a resolution of 0.1 degree:
#' scan_to_spatraster(example_scan, ylim = c(55, 57), xlim = c(12, 13), res = .1)
#'
#' # using a template raster (a terra SpatRaster can be passed directly):
#' template_raster <- terra::rast(terra::ext(12, 13, 56, 58), crs = "epsg:4326")
#' scan_to_spatraster(example_scan, raster = template_raster)
#' }
scan_to_spatraster <- function(scan, nx = 100, ny = 100, xlim, ylim, res = NA, param, raster = NA, lat, lon, crs = NA, k = 4 / 3, re = 6378, rp = 6357) {
  if (!is.scan(scan)) stop("'scan' should be an object of class scan")
  if (get_elevation_angles(scan) == 90) stop("georeferencing of 90 degree birdbath scan not supported")
  if (!identical(raster, NA)) raster <- .spatial_raster(raster)
  if (!assertthat::is.number(nx) && missing(res)) stop("'nx' should be an integer")
  if (!assertthat::is.number(ny) && missing(res)) stop("'ny' should be an integer")
  if (!missing(xlim)) {
    if (length(xlim) != 2 & !is.numeric(xlim)) stop("'xlim' should be an integer vector of length two")
    if (is.na(xlim[1]) | is.na(xlim[2]) | xlim[1] > xlim[2]) stop("'xlim' should be a vector with two numeric values for upper and lower bound")
  }
  if (!missing(ylim)) {
    if (length(ylim) != 2 & !is.numeric(ylim)) stop("'ylim' should be an integer vector of length two")
    if (is.na(ylim[1]) | is.na(ylim[2]) | ylim[1] > ylim[2]) stop("'ylim' should be a vector with two numeric values for upper and lower bound")
  }
  if (!missing(res) && !all(is.na(res))) {
    assertthat::assert_that(is.numeric(res))
    assertthat::assert_that(length(res) <= 2)
  }
  if (!missing(param)) {
    if (FALSE %in% (param %in% c(names(scan$params), "azim", "range", "distance"))) stop("'param' contains scan parameter not found in scan")
    if (!(FALSE %in% (param %in% c("azim", "range", "distance")))) stop("'param' should contain the name of one or more scan parameters contained in 'scan'")

    param_to_use <- param
  }
  else {
    param_to_use <- names(scan$params)
  }
  if (!assertthat::are_equal(raster, NA)) {
    assertthat::assert_that(inherits(raster, "SpatRaster"))
  }

  if (is.null(scan$geo$lat) && missing(lat)) stop("radar latitude cannot be found in scan, specify using 'lat' argument")
  if (is.null(scan$geo$lon) && missing(lon)) stop("radar longitude cannot be found in scan, specify using 'lon' argument")

  if (missing(lat)) lat <- scan$geo$lat
  if (missing(lon)) lon <- scan$geo$lon

  assertthat::assert_that(assertthat::is.number(lat))
  assertthat::assert_that(assertthat::is.number(lon))
  localCrs <- paste("+proj=aeqd +lat_0=", lat,
    " +lon_0=", lon,
    " +units=m",
    sep = ""
  )
  if (missing(crs) || (length(crs) == 1 && is.na(crs))) {
    crs <- localCrs
  }
  else {
    # normalize modern and legacy CRS inputs
    crs <- .spatial_crs(crs)
  }
  if (!assertthat::are_equal(raster, NA)) {
    crs <- terra::crs(raster)
  }
  assertthat::assert_that(assertthat::is.number(k))
  assertthat::assert_that(assertthat::is.number(re))
  assertthat::assert_that(assertthat::is.number(rp))

  rscale <- scan$geo$rscale
  ascale <- scan$geo$ascale
  rstart <- ifelse(is.null(scan$geo$rstart), 0, scan$geo$rstart)
  astart <- ifelse(is.null(scan$geo$astart), 0, scan$geo$astart)


  nrang <- dim(scan)[2]
  nazim <- dim(scan)[3]

  # extent not fully specified, determine it
  if (identical(raster, NA) && (missing(xlim) || missing(ylim))) {
    # georeference the data
    points <- scan_to_sf(scan, k = k, lat = lat, lon = lon, re = re, rp = rp)
    # keep only selected scan parameters
    if (!missing(param)) points <- points[param]
    # transform spatialpoints to coordinate system of the raster
    points <- sf::st_transform(points, sf::st_crs(crs))
    # get extent of the available data
    points_extent <- sf::st_bbox(points)
    # prepare a raster matching the data extent (or user-specified extent)
    if (missing(xlim)) xlim <- points_extent[c("xmin", "xmax")]
    if (missing(ylim)) ylim <- points_extent[c("ymin", "ymax")]
  }
  if (!assertthat::are_equal(raster, NA)) {
    r <- terra::rast(raster, nlyrs = 1)
  } else {
    if (missing(res) || all(is.na(res))) {
      r <- terra::rast(ncols = nx, nrows = ny, ext = terra::ext(c(xlim, ylim)), crs = crs)
    }
    else {
      r <- terra::rast(ext = terra::ext(c(xlim, ylim)), crs = crs, resolution = res)
    }
  }

  crds <- sf::st_coordinates(.raster_points(r, localCrs))

  # convert raster coordinates to polar indices
  polar_coords <- cartesian_to_polar(crds, elev = scan$geo$elangle, k = k, lat = lat, re = re, rp = rp)
  index <- polar_to_index(polar_coords, rangebin = rscale, azimbin = ascale, rangestart = rstart, azimstart = astart)
  # set indices outside the scan's matrix to NA
  index$row[index$row > nrang] <- NA
  index$col[index$col > nazim] <- NA
  # rstart can result in locations outside of the radar scope close to the radar
  index$row[index$row < 1] <- NA
  stopifnot(all(index$col >= 1))

  # convert 2D index to 1D index
  index <- (index$col - 1) * nrang + index$row

  # generate brick for each scan parameter
  param_to_add <- setdiff(param_to_use, c("distance", "range", "azim"))
  output <- terra::rast(r, nlyrs = length(param_to_add))
  names(output) <- param_to_add
  terra::values(output) <- vapply(param_to_add, function(name) as.vector(scan$params[[name]][index]), numeric(length(index)))

  if ("distance" %in% param_to_use) output$distance <- beam_distance(polar_coords$range, elev = scan$geo$elangle, k = k, lat = lat, re = re, rp = rp)
  if ("range" %in% param_to_use) output$range <- polar_coords$range
  if ("azim" %in% param_to_use) output$azim <- polar_coords$azim
  output
}

#' Convert a polar scan into a raster
#'
#' `r lifecycle::badge("deprecated")`
#'
#' Retained for compatibility during the spatial migration. Use the corresponding
#' `sf` or `terra` conversion for new code. The legacy package must be installed.
#'
#' @inheritParams scan_to_spatraster
#' @returns A `RasterBrick`.
#' @family deprecated spatial functions
#' @export
scan_to_raster <- function(scan, nx = 100, ny = 100, xlim, ylim, res = NA, param, raster = NA, lat, lon, crs = NA, k = 4 / 3, re = 6378, rp = 6357) {
  lifecycle::deprecate_warn("1.0.0", "scan_to_raster()", "scan_to_spatraster()")
  call <- match.call()
  call[[1]] <- quote(scan_to_spatraster)
  rlang::check_installed("raster", reason = "to return a legacy RasterBrick")
  raster::brick(eval(call, parent.frame()))
}

# Internal point sampler used by the deprecated scan_to_spdf() adapter.
# allows for faster projections of multiple scans to the same grid in
# integrate_to_ppi function, see issue #293
.scan_to_spdf <- function(scan, spdf, param, lat, lon, k = 4 / 3, re = 6378, rp = 6357) {
  if (!is.scan(scan)) stop("'scan' should be an object of class scan")
  if (get_elevation_angles(scan) == 90) stop("georeferencing of 90 degree birdbath scan not supported")
  if (!inherits(spdf, "sf")) {
    stop("spdf should be of class sf")
  }
  if (!missing(param)) {
    if (FALSE %in% (param %in% c(names(scan$params), "azim", "range", "distance"))) stop("'param' contains scan parameter not found in scan")
    if (!(FALSE %in% (param %in% c("azim", "range", "distance")))) stop("'param' should contain the name of one or more scan parameters contained in 'scan'")

    param_to_use <- param
  }
  else {
    param_to_use <- names(scan$params)
  }

  if (is.null(scan$geo$lat) && missing(lat)) stop("radar latitude cannot be found in scan, specify using 'lat' argument")
  if (is.null(scan$geo$lon) && missing(lon)) stop("radar longitude cannot be found in scan, specify using 'lon' argument")

  if (missing(lat)) lat <- scan$geo$lat
  if (missing(lon)) lon <- scan$geo$lon

  assertthat::assert_that(assertthat::is.number(lat))
  assertthat::assert_that(assertthat::is.number(lon))
  localCrs <- paste("+proj=aeqd +lat_0=", lat,
    " +lon_0=", lon,
    " +units=m",
    sep = ""
  )
  assertthat::assert_that(assertthat::is.number(k))
  assertthat::assert_that(assertthat::is.number(re))
  assertthat::assert_that(assertthat::is.number(rp))
  stopifnot(sf::st_crs(localCrs) == sf::st_crs(spdf))

  rscale <- scan$geo$rscale
  ascale <- scan$geo$ascale
  rstart <- ifelse(is.null(scan$geo$rstart), 0, scan$geo$rstart)
  astart <- ifelse(is.null(scan$geo$astart), 0, scan$geo$astart)

  nrang <- dim(scan)[2]
  nazim <- dim(scan)[3]

  crds <- sf::st_coordinates(spdf)
  # convert raster coordinates to polar indices
  polar_coords <- cartesian_to_polar(crds, elev = scan$geo$elangle, k = k, lat = lat, re = re, rp = rp)
  index <- polar_to_index(polar_coords, rangebin = rscale, azimbin = ascale, rangestart = rstart, azimstart = astart)
  # set indices outside the scan's matrix to NA
  index$row[index$row > nrang] <- NA
  index$col[index$col > nazim] <- NA
  # rstart can result in locations outside of the radar scope close to the radar
  index$row[index$row < 1] <- NA
  stopifnot(all(index$col >= 1))
  # convert 2D index to 1D index
  index <- (index$col - 1) * nrang + index$row

  # generate brick for each scan parameter
  param_to_add <- setdiff(param_to_use, c("distance", "range", "azim"))
  output <- spdf[, F]
  # fill the rasterbrick
  for (name in param_to_add) {
    output[[name]] <- scan$params[[name]][index]
  }
  if ("distance" %in% param_to_use) output$distance <- beam_distance(polar_coords$range, elev = scan$geo$elangle, k = k, lat = lat, re = re, rp = rp)
  if ("range" %in% param_to_use) output$range <- polar_coords$range
  if ("azim" %in% param_to_use) output$azim <- polar_coords$azim
  output
}

scan_to_spdf <- function(scan, spdf, param, lat, lon, k = 4 / 3, re = 6378, rp = 6357) {
  lifecycle::deprecate_warn("1.0.0", "scan_to_spdf()", "scan_to_sf()")
  call <- match.call()
  call[[1]] <- quote(.scan_to_spdf)
  rlang::check_installed("sp", reason = "to return a legacy Spatial object")
  call$spdf <- sf::st_as_sf(spdf)
  sf::as_Spatial(eval(call, parent.frame()))
}
