#' Project a scan (`scan`) or parameter (`param`) to a plan
#' position indicator (`ppi`)
#'
#' Make a plan position indicator (ppi)
#'
#' @inheritParams beam_height
#' @param x An object of class `param` or `scan`.
#' @param grid_size Cartesian grid size in m.
#' @param range_max Maximum range in m.
#' @param ylim The range of latitudes to include.
#' @param xlim The range of longitudes to include.
#' @param project Whether to vertically project onto earth's surface.
#' @param raster (optional) `terra::SpatRaster` (legacy `RasterLayer` inputs are deprecated) with a CRS. When specified
#' this raster topology is used for the output, and `grid_size`, `range_max`, `xlim`, `ylim` are ignored.
#' @returns An object of class '[ppi][summary.ppi]'.
#' @family ppi calculation functions
#' @export
#' @details The returned PPI is in Azimuthal Equidistant Projection.
#' @examples
#' # load a polar scan example object:
#' data(example_scan)
#' example_scan
#'
#' # plot the scan:
#' plot(example_scan)
#'
#' # make PPIs for all scan parameters in the scan:
#' ppi <- project_as_ppi(example_scan)
#'
#' # print summary info for the ppi:
#' ppi
#'
#' # plot the ppi:
#' plot(ppi)
#'
#' # extract the DBZH scan parameter of the volume to a new
#' # object 'param':
#' param <- get_param(example_scan, "VRADH")
#'
#' # make a ppi for the new 'param' object:
#' ppi <- project_as_ppi(param)
#'
#' # print summary info for this ppi:
#' ppi
#'
#' # plot the ppi:
#' plot(ppi)
project_as_ppi <- function(x, grid_size = 500, range_max = 50000,
                           project = TRUE, ylim = NULL, xlim = NULL, raster = NA, k = 4 / 3, re = 6378, rp = 6357) {
  UseMethod("project_as_ppi", x)
}


#' @describeIn project_as_ppi Project as `ppi` for a single scan parameter.
#' @export
project_as_ppi.param <- function(x, grid_size = 500, range_max = 50000,
                                 project = TRUE, ylim = NULL, xlim = NULL, raster = NA, k = 4 / 3, re = 6378, rp = 6357) {
  stopifnot(inherits(x, "param"))
  if (!identical(raster, NA)) grid_size <- .spatial_raster(raster)
  data <- sample_polar(x, grid_size, range_max, project, ylim, xlim, k = k, re = re, rp = rp)
  geo <- attributes(x)$geo
  geo$bbox <- attr(data, "bboxlatlon")
  if (is.null(geo$bbox)) geo$bbox <- .raster_bbox(data)
  geo$merged <- FALSE
  structure(list(radar = attr(x, "radar"), datetime = attr(x, "datetime"),
    data = data, geo = geo), class = "ppi")
}

#' @describeIn project_as_ppi Project multiple `ppi`'s for all scan
#' parameters in a scan
#' @export
project_as_ppi.scan <- function(x, grid_size = 500, range_max = 50000,
                                project = TRUE, ylim = NULL, xlim = NULL, raster = NA, k = 4 / 3, re = 6378, rp = 6357) {
  stopifnot(inherits(x, "scan"))

  if (!identical(raster, NA)) grid_size <- .spatial_raster(raster)
  layers <- lapply(x$params, sample_polar, grid_size = grid_size, range_max = range_max,
    project = project, ylim = ylim, xlim = xlim, k = k, re = re, rp = rp)
  data <- do.call(c, unname(layers))
  geo <- x$geo
  geo$bbox <- attr(layers[[1]], "bboxlatlon")
  if (is.null(geo$bbox)) geo$bbox <- .raster_bbox(data)
  geo$merged <- FALSE
  structure(list(radar = x$radar, datetime = x$datetime, data = data, geo = geo), class = "ppi")
}


sample_polar <- function(param, grid_size, range_max, project, ylim, xlim, k = 4 / 3, re = 6378, rp = 6357) {
  proj4string <- (paste("+proj=aeqd +lat_0=", attributes(param)$geo$lat,
    " +lon_0=", attributes(param)$geo$lon,
    " +units=m",
    sep = ""
  ))
  if (inherits(grid_size, c("Raster", "SpatRaster"))) {
    output <- terra::rast(.spatial_raster(grid_size), nlyrs = 1)
    gridTopoCrds <- sf::st_coordinates(.raster_points(output, proj4string))
  } else {
    bboxlatlon <- matrix(sf::st_bbox(proj_to_wgs(c(-range_max, range_max),
      c(-range_max, range_max), proj4string)), nrow = 2,
      dimnames = list(c("lon", "lat"), c("min", "max")))
    if (!is.null(ylim)) bboxlatlon["lat", ] <- ylim
    if (!is.null(xlim)) bboxlatlon["lon", ] <- xlim
    bbox <- sf::st_coordinates(wgs_to_proj(bboxlatlon["lon", ], bboxlatlon["lat", ], proj4string))
    # Without geographic limits, retain the exact requested projected bounds.
    # A projection round trip can otherwise add a cell through ceiling().
    if (is.null(xlim) && is.null(ylim)) {
      bbox <- cbind(c(-range_max, range_max), c(-range_max, range_max))
    }
    offset <- apply(bbox, 2, min)
    cells <- apply(bbox, 2, function(z) diff(range(z))) / grid_size
    nearest <- round(cells)
    close <- abs(cells - nearest) < 1e-10
    cells[close] <- nearest[close]
    dims <- ceiling(cells)
    output <- terra::rast(ncols = dims[1], nrows = dims[2],
      xmin = offset[1] - grid_size / 2, xmax = offset[1] + (dims[1] - .5) * grid_size,
      ymin = offset[2] - grid_size / 2, ymax = offset[2] + (dims[2] - .5) * grid_size,
      crs = proj4string)
    gridTopoCrds <- terra::xyFromCell(output, seq_len(terra::ncell(output)))
  }
  # if projecting, account for elevation angle
  if (project) {
    elev <- attributes(param)$geo$elangle
  } else {
    elev <- 0
  }

  # get scan parameter indices, and extract data
  index <- polar_to_index(
    cartesian_to_polar(gridTopoCrds, elev, k = k, lat = attributes(param)$geo$lat, re = re, rp = rp),
    rangebin = attributes(param)$geo$rscale,
    azimbin = attributes(param)$geo$ascale,
    rangestart = max(c(0, attributes(param)$geo$rstart), na.rm = TRUE),
    azimstart = max(c(0, attributes(param)$geo$astart), na.rm = TRUE)
  )

  # set indices outside the scan's matrix to NA
  nrang <- dim(param)[1]
  nazim <- dim(param)[2]
  index$row[index$row > nrang] <- NA
  index$col[index$col > nazim] <- NA
  # rstart can result in locations outside of the radar scope close to the radar
  index$row[index$row < 1] <- NA
  stopifnot(all(index$col >= 1))
  # convert 2D index to 1D index
  index <- (index$col - 1) * nrang + index$row
  data <- as.data.frame(param[index])

  colnames(data) <- attributes(param)$param

  terra::values(output) <- data
  names(output) <- colnames(data)
  if (exists("bboxlatlon", inherits = FALSE)) attr(output, "bboxlatlon") <- bboxlatlon

  output
}

cartesian_to_polar <- function(coords, elev = 0, k = 4 / 3, lat = 35, re = 6378, rp = 6357) {
  range <- beam_range(sqrt(coords[, 1]^2 + coords[, 2]^2), elev, k = k, lat = lat, re = re, rp = rp)
  azim <- (0.5 * pi - atan2(coords[, 2], coords[, 1])) %% (2 * pi)
  data.frame(range = range, azim = azim * 180 / pi)
}

safe_subset <- function(data, indexx, indexy) {
  datadim <- dim(data)
  if (indexx < 1 || indexx > datadim[1] || indexy < 1 || indexy > datadim[2]) {
    out <- NA
  } else {
    out <- data[indexx, indexy]
  }
  out
}

polar_to_index <- function(coords_polar, rangebin = 1, azimbin = 1, rangestart = 0, azimstart = 0) {
  row <- floor(1 + (coords_polar$range - rangestart) / c(rangebin))
  azimuth_bin <- ((coords_polar$azim - azimstart + 360) %% 360) / azimbin
  # Give points on bin boundaries the same bin across PROJ/platform round-off.
  nearest <- round(azimuth_bin)
  close <- abs(azimuth_bin - nearest) < 1e-10
  azimuth_bin[close] <- nearest[close]
  col <- 1 + floor(azimuth_bin %% (360 / azimbin))
  data.frame(row = row, col = col)
}
