# Load the polar volume that ships with bioRad as example data
pvolfile <- system.file("extdata", "volume.h5", package = "bioRad")
pvol <- read_pvolfile(pvolfile)

# Build a terra SpatRaster matching the topology of the lowest scan and fill it
# with synthetic data. scan_to_spatraster() returns a terra object in the
# radar's local azimuthal-equidistant projection.
make_synthetic_spatraster <- function(scan) {
  r <- scan_to_spatraster(scan, nx = 50, ny = 50, param = "DBZH")
  # synthetic, fully populated values so extraction within the extent is non-NA
  terra::values(r) <- seq_len(terra::ncell(r))
  names(r) <- "synthetic"
  r
}

lowest_scan <- get_scan(pvol, 0.5)
# Modern workflows run without the optional legacy packages. Compatibility
# with RasterLayer inputs is tested below.
synthetic_spatraster <- make_synthetic_spatraster(lowest_scan)
synthetic_raster <- synthetic_spatraster

test_that("add_param() requires its arguments", {
  # `param` has no default and must be supplied
  expect_error(
    add_param(lowest_scan, synthetic_raster),
    regexp = "missing subscript",
    fixed = TRUE
  )
  # the generic only dispatches for `scan` and `pvol` objects
  expect_error(
    add_param("not_a_scan_or_pvol", synthetic_spatraster, "synthetic"),
    regexp = "no applicable method for 'add_param'"
  )
})

test_that("add_param.scan() errors when x is not a scan", {
  expect_error(
    add_param.scan("not_a_scan", synthetic_raster, "synthetic"),
    regexp = 'inherits(x, "scan") is not TRUE',
    fixed = TRUE
  )
})

test_that("add_param() adds a parameter to a scan (SpatRaster input)", {
  scan_out <- add_param(lowest_scan, synthetic_raster, "synthetic")

  # the returned object is still a scan
  expect_s3_class(scan_out, "scan")
  # the new parameter has been added under the requested name
  expect_true("synthetic" %in% names(scan_out$params))
  # the added parameter is a proper scan parameter
  expect_s3_class(scan_out$params[["synthetic"]], "param")
  # it shares the dimensions of the existing scan parameters
  expect_equal(
    dim(scan_out$params[["synthetic"]]),
    dim(scan_out$params[[1]])
  )
  # the `param` attribute is set to the supplied name
  expect_equal(attributes(scan_out$params[["synthetic"]])$param, "synthetic")
  # at least some range gates fall within the raster and receive a value
  expect_false(all(is.na(scan_out$params[["synthetic"]])))
})

test_that("add_param() honours the param argument name (SpatRaster input)", {
  scan_out <- add_param(lowest_scan, synthetic_spatraster, "HGHT")
  expect_s3_class(scan_out, "scan")
  expect_true("HGHT" %in% names(scan_out$params))
  expect_equal(attributes(scan_out$params[["HGHT"]])$param, "HGHT")
})

test_that("add_param() adds a parameter to every scan of a pvol (SpatRaster input)", {
  expect_s4_class(synthetic_raster, "SpatRaster")

  pvol_out <- add_param(pvol, synthetic_raster, "synthetic")

  # the returned object is still a pvol
  expect_s3_class(pvol_out, "pvol")
  # every scan in the volume now carries the new parameter
  expect_true(all(vapply(
    pvol_out$scans,
    function(s) "synthetic" %in% names(s$params),
    logical(1)
  )))
  expect_equal(
    attributes(get_scan(pvol_out, 0.5)$params[["synthetic"]])$param,
    "synthetic"
  )
})

test_that("add_param() accepts a terra SpatRaster for a pvol", {
  expect_s4_class(synthetic_spatraster, "SpatRaster")

  pvol_out <- add_param(pvol, synthetic_spatraster, "synthetic")

  expect_s3_class(pvol_out, "pvol")
  expect_true(all(vapply(
    pvol_out$scans,
    function(s) "synthetic" %in% names(s$params),
    logical(1)
  )))
})

test_that("add_param() extracts geographic raster values in the correct CRS", {
  dem <- terra::rast(nrows = 10, ncols = 10, xmin = 12, xmax = 14,
    ymin = 55, ymax = 58, crs = "EPSG:4326")
  terra::values(dem) <- 123
  scan <- add_param(example_scan, dem, "HGHT")
  points <- sf::st_transform(scan_to_sf(example_scan), 4326)
  xy <- sf::st_coordinates(points)
  inside <- xy[,1] > 12 & xy[,1] < 14 & xy[,2] > 55 & xy[,2] < 58
  expect_true(any(inside))
  expect_equal(as.vector(scan$params$HGHT)[inside], rep(123, sum(inside)))
})

test_that("add_param() warns for legacy raster arguments", {
  skip_if_not_installed("raster")
  template <- terra::rast(nrows = 3, ncols = 3, xmin = 12.8, xmax = 12.9,
    ymin = 56.3, ymax = 56.4, crs = "EPSG:4326")
  terra::values(template) <- 100
  lifecycle::expect_deprecated(add_param(example_scan, raster::raster(template), "HGHT"))
})

test_that("add_param() preserves polar-volume extraction with legacy rasters", {
  skip_if_not_installed("raster")
  pvol <- read_pvolfile(system.file("extdata", "volume.h5", package = "bioRad"))
  pvol$scans <- pvol$scans[1:2]
  dem <- scan_to_spatraster(pvol$scans[[1]], nx = 10, ny = 10, param = "DBZH")
  terra::values(dem) <- seq_len(terra::ncell(dem))
  modern <- add_param(pvol, dem, "HGHT")
  lifecycle::expect_deprecated(legacy <- add_param(pvol, raster::raster(dem), "HGHT"))
  expect_equal(lapply(legacy$scans, function(x) x$params$HGHT),
    lapply(modern$scans, function(x) x$params$HGHT))
})
