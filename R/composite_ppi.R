#' Create a composite of multiple plan position indicators (`ppi`)
#'
#' Combines multiple plan position indicators (`ppi`) into a single
#' `ppi`. Can be used to make a composite of `ppi`'s from multiple
#' radars.
#'
#' @inheritParams integrate_to_ppi
#' @param x A list of `ppi` objects.
#' @param param Character (vector). One or more parameter name(s) to composite.
#'   To composite all available scan parameters use `all` (default).
#' @param method Character (vector). Compositing method(s), either `mean`,
#'   `min`, `max` or `idw`. To apply different methods for each of the
#'   parameters, provide a vector with the same length as `param`.
#' @param idw_max_distance Numeric. Maximum distance from the radar to consider
#'   in inverse distance weighting. Measurements beyond this distance will have
#'   a weighting factor of zero.
#' @param idp Numeric. Inverse distance weighting power.
#' @param coverage Logical. When `TRUE`, adds an additional `coverage` parameter
#'   to the `ppi` indicating the number of `ppi`s covering a single composite
#'   pixel.
#' @returns A `ppi` object.
#' @family ppi calculation functions
#' @export
#' @details
#' The function can combine multiple `ppi`s of different scan elevations of the
#' same radar or `ppi`s of different radars. The coordinates of the returned
#' `ppi` object are in the WGS84 datum, unless a different `crs` is provided. If
#' only `res` is provided, but no `crs` is set, `res` is in meters and the
#' origin of the composite `ppi` is set to the mean(lat, lon) location.
#'
#' The `method` parameter determines how values of different `ppi`s at the same
#' geographic location are combined:
#' * `mean`: Compute the average value.
#' * `max`: Compute the maximum value. If `ppi`s are of the same radar and the
#' same polar volume, this computes a max product, showing the maximum detected
#' signal at that geographic location.
#' * `min`: Compute the minimum value.
#' * `idw`: This option is useful primarily when compositing `ppi`s of multiple
#' radars. Performs an inverse distance weighting, where values are weighted
#' according to 1/(distance from the radar)^`idp`.
#'
#' Argument `method` determines how values of different ppi's at the same
#' geographic location are combined.
#'
#' * `mean`: Compute the average value
#' * `max`: Compute the maximum value. If ppi's are of the same radar and the
#'   same polar volume, this computes a max product, showing the maximum
#'   detected signal at that geographic location.
#' * `min`: Compute the minimum value
#' * `idw`: This option is useful primarily when compositing ppi's of multiple
#'   radars. Performs an inverse distance weighting, where values are weighted
#'   according to 1/(distance from the radar)^`idp`
#'
#' The coordinates system of the returned `ppi` is a WGS84
#' (lat, lon) datum, unless a different `crs` is provided. If only
#' `res` is provided, but no `crs` is set, `res` is in
#' meter units and the origin of the composite `ppi` is set to the
#' mean (lat, lon) location.
#'
#' This function is a prototype and under active development
#' @examples
#' # Locate and read the polar volume example file
#' \donttest{
#' pvolfile <- system.file("extdata", "volume.h5", package = "bioRad")
#' pvol <- read_pvolfile(pvolfile)
#'
#' # Calculate a ppi for each elevation scan, 1000m grid resolution
#' ppis <- lapply(pvol$scans, project_as_ppi, grid_size=1000)
#'
#' # Overlay the ppis, calculating the maximum value observed
#' # across the available scans at each geographic location
#' composite <- composite_ppi(ppis, method = "max", res=1000)
#'
#' # Plot the calculated max product on the basemap
#' if (all(sapply(c("ggspatial", "prettymapr", "rosm"), requireNamespace, quietly = TRUE))) {
#' map(composite)
#' }
#' }
composite_ppi <-
  function(x,
           param = "all",
           nx = 100,
           ny = 100,
           xlim,
           ylim,
           res,
           crs,
           raster = NA,
           method = "max",
           idp = 2,
           idw_max_distance = NA,
           coverage = FALSE) {
  if (FALSE %in% sapply(x, is.ppi)) {
    stop("'x' should be an object of class ppi, or a list of objects of class ppi")
  }
  if (!assertthat::is.count(nx) && missing(res)) stop("'nx' should be an integer")
  if (!assertthat::is.count(ny) && missing(res)) stop("'ny' should be an integer")
  if (!missing(xlim)) {
    if (length(xlim) != 2 ||
        !is.numeric(xlim)) {
      stop("'xlim' should be a numeric vector of length two")
    }
    if (is.na(xlim[1]) |
        is.na(xlim[2]) |
        xlim[1] > xlim[2]) {
      stop(paste("'xlim' should be a vector with two numeric values",
                 "for lower and upper bound respectively"))
    }
  }
  if (!missing(ylim)) {
    if (length(ylim) != 2 ||
        !is.numeric(ylim)) {
      stop("'ylim' should be a numeric vector of length two")
    }
    if (is.na(ylim[1]) |
        is.na(ylim[2]) |
        ylim[1] > ylim[2]) {
      stop(
        paste(
          "'ylim' should be a vector with two numeric values for",
          "lower and upper bound respectively"
        )
      )
    }
  }
  if (!missing(res)) {
    assertthat::assert_that(is.numeric(res))
    assertthat::assert_that(length(res) <= 2)
    t_res <- res
  } else {
    t_res <- NULL
  }
  # normalize modern and legacy CRS inputs
  if (!missing(crs)) {
    t_crs <- .spatial_crs(crs)
  } else {
    t_crs <- NULL
  }

  if (!all(method %in% c("max", "min", "mean", "idw"))) stop("'method' should be one or multiple of 'max', 'mean', 'min' or 'idw'")
  if (length(method) != length(param) & length(method) != 1) stop("'method' should be of length 1 or length(param)")
  assertthat::assert_that(assertthat::is.flag(coverage))

  if (length(param) == 1 && param == "all") {
    param <- names(x[[1]]$data)
  }
  ppis <- lapply(x, `[.ppi`, i = param)

  lons <- sapply(ppis, function(x) x$geo$bbox["lon", ])
  lats <- sapply(ppis, function(x) x$geo$bbox["lat", ])
  if (!missing(xlim)) {
    lons <- xlim
  }
  if (!missing(ylim)) {
    lats <- ylim
  }
  lons.radar <- sapply(ppis, function(x) x$geo$lon)
  lats.radar <- sapply(ppis, function(x) x$geo$lat)
  elangles <- sapply(ppis, function(x) x$geo$elangle)
  bbox <- matrix(c(min(lons), min(lats), max(lons), max(lats)),
    nrow = 2,
    ncol = 2, dimnames = dimnames(ppis[[1]]$geo$bbox)
  )

  if (!identical(raster, NA)) {
    r <- terra::rast(.spatial_raster(raster), nlyrs = 1)
  } else {
    r <- terra::rast(ncols = nx, nrows = ny,
      ext = terra::ext(min(lons), max(lons), min(lats), max(lats)), crs = "EPSG:4326")
    if (is.null(t_crs) && !is.null(t_res)) {
      t_crs <- paste0("+proj=aeqd +units=m +ellps=WGS84 +lat_0=", mean(lats), " +lon_0=", mean(lons))
    }
    if (!is.null(t_crs)) r <- terra::project(r, t_crs)
    if (!is.null(t_res)) terra::res(r) <- t_res
  }
  assertthat::assert_that(terra::ncol(r) > 1, msg = "output raster has only a single pixel in x dimension, increase raster resolution")
  assertthat::assert_that(terra::nrow(r) > 1, msg = "output raster has only a single pixel in y dimension, increase raster resolution")
  points <- .raster_points(r)
  # Sample source cells at target centres, preserving the original nearest-cell semantics.
  projs <- lapply(ppis, function(x) {
    source <- .spatial_raster(x$data)
    pts <- terra::vect(sf::st_transform(points, terra::crs(source)))
    data <- terra::extract(source, pts, ID = FALSE)
    if (coverage) {
      xy <- terra::crds(pts)
      data$coverage <- as.numeric(!is.na(terra::cellFromXY(source, xy)))
    }
    data
  })
  result <- list()
  for (p in param) {
    values <- do.call(cbind, lapply(projs, `[[`, p))
    param_method <- if (length(method) == 1) method else method[match(p, param)]
    if (param_method == "idw") {
      # Distance in metres, including when the target CRS is geographic.
      weights <- vapply(ppis, function(x) {
        radar <- sf::st_as_sf(data.frame(lon = x$geo$lon[1], lat = x$geo$lat[1]),
          coords = c("lon", "lat"), crs = 4326)
        d <- as.numeric(sf::st_distance(sf::st_transform(points, 4326), radar))
        if (!is.na(idw_max_distance)) d[d > idw_max_distance] <- NA_real_
        1 / d^idp
      }, numeric(nrow(points)))
      weights[is.na(values)] <- NA_real_
      # At a radar centre use only the coincident radar(s).
      zero <- rowSums(is.infinite(weights)) > 0
      weights[zero, ] <- ifelse(is.infinite(weights[zero, , drop = FALSE]), 1, 0)
      result[[p]] <- rowSums(values * weights, na.rm = TRUE) / rowSums(weights, na.rm = TRUE)
    } else {
      result[[p]] <- switch(param_method,
        mean = rowMeans(values, na.rm = TRUE),
        max = do.call(pmax, c(as.data.frame(values), na.rm = TRUE)),
        min = do.call(pmin, c(as.data.frame(values), na.rm = TRUE)))
    }
    result[[p]][rowSums(!is.na(values)) == 0] <- NA_real_
  }
  if (coverage) result$coverage <- rowSums(do.call(cbind, lapply(projs, `[[`, "coverage")))
  output <- terra::rast(r, nlyrs = length(result))
  terra::values(output) <- as.data.frame(result)
  names(output) <- names(result)
  structure(list(data = output, geo = list(lat = lats.radar, lon = lons.radar,
    elangle = elangles, bbox = .raster_bbox(output), merged = TRUE)), class = "ppi")
}
