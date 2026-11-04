test_that("legacy scan conversion functions are deprecated", {
  skip_if_not_installed("sp")
  skip_if_not_installed("raster")
  points <- sf::as_Spatial(scan_to_sf(example_scan))

  lifecycle::expect_deprecated(scan_to_spatial(example_scan))
  lifecycle::expect_deprecated(scan_to_raster(example_scan, nx = 10, ny = 10))
  lifecycle::expect_deprecated(scan_to_spdf(example_scan, points, param = "DBZH"))
})

test_that("scan_to_sf() returns error on incorrect parameters", {
  expect_error(scan_to_sf(scan = "a"),
               regexp = "is.scan(x = scan) is not TRUE",
               fixed = TRUE)
  expect_error(scan_to_sf(example_scan, k = "a"),
               regexp = "k is not a number (a length one numeric vector).",
               fixed = TRUE)
  expect_error(scan_to_sf(example_scan, re = "a"),
               regexp = "re is not a number (a length one numeric vector).",
               fixed = TRUE)
  expect_error(scan_to_sf(example_scan, rp = "a"),
               regexp = "rp is not a number (a length one numeric vector).",
               fixed = TRUE)
  expect_error(scan_to_sf(example_scan, lat = "a"),
               regexp = "lat is not a number (a length one numeric vector).",
               fixed = TRUE)
  expect_error(scan_to_sf(example_scan, lon = "a"),
               regexp = "lon is not a number (a length one numeric vector).",
               fixed = TRUE)
  missing_lat_scan <- example_scan
  missing_lon_scan <- example_scan
  missing_lat_scan$geo$lat <- NULL
  missing_lon_scan$geo$lon <- NULL
  expect_error(scan_to_sf(missing_lat_scan),
               regexp = "radar latitude cannot be found in scan, specify using 'lat' argument",
               fixed = TRUE)
  expect_error(scan_to_sf(missing_lon_scan),
               regexp = "radar longitude cannot be found in scan, specify using 'lon' argument",
               fixed = TRUE)
})

test_that("scan_to_spatraster() returns error on incorrect parameters", {
  expect_error(scan_to_spatraster("a"),
               regexp = "'scan' should be an object of class scan",
               fixed = TRUE)
  square_birdbath <- example_scan
  square_birdbath$attributes$where$elangle <- 90
  expect_error(scan_to_spatraster(square_birdbath),
               regexp = "georeferencing of 90 degree birdbath scan not supported",
               fixed = TRUE)
  expect_error(scan_to_spatraster(example_scan, nx = "a"),
               regexp = "'nx' should be an integer",
               fixed = TRUE)
  expect_error(scan_to_spatraster(example_scan, ny = "a"),
               regexp = "'ny' should be an integer",
               fixed = TRUE)
  expect_error(scan_to_spatraster(example_scan, xlim = "a"),
               regexp = "'xlim' should be an integer vector of length two",
               fixed = TRUE)
  expect_error(scan_to_spatraster(example_scan, xlim = 2),
               regexp = "'xlim' should be a vector with two numeric values for upper and lower bound",
               fixed = TRUE)
  expect_error(scan_to_spatraster(example_scan, xlim = c("a",2)),
               regexp = "'xlim' should be a vector with two numeric values for upper and lower bound",
               fixed = TRUE)
  expect_error(scan_to_spatraster(example_scan, xlim = c(pi, 2)),
               regexp = "'xlim' should be a vector with two numeric values for upper and lower bound",
               fixed = TRUE)
  expect_error(scan_to_spatraster(example_scan, ylim = "a"),
               regexp = "'ylim' should be an integer vector of length two",
               fixed = TRUE)
  expect_error(scan_to_spatraster(example_scan, ylim = 2),
               regexp = "'ylim' should be a vector with two numeric values for upper and lower bound",
               fixed = TRUE)
  expect_error(scan_to_spatraster(example_scan, ylim = c("a",2)),
               regexp = "'ylim' should be a vector with two numeric values for upper and lower bound",
               fixed = TRUE)
  expect_error(scan_to_spatraster(example_scan, ylim = c(pi, 2)),
               regexp = "'ylim' should be a vector with two numeric values for upper and lower bound",
               fixed = TRUE)
  missing_lat_scan <- example_scan
  missing_lon_scan <- example_scan
  missing_lat_scan$geo$lat <- NULL
  missing_lon_scan$geo$lon <- NULL
  expect_error(scan_to_spatraster(missing_lat_scan),
               regexp = "radar latitude cannot be found in scan, specify using 'lat' argument",
               fixed = TRUE)
  expect_error(scan_to_spatraster(missing_lon_scan),
               regexp = "radar longitude cannot be found in scan, specify using 'lon' argument",
               fixed = TRUE)
})

test_that("scan_to_spatraster() returns error on wrongly formed param argument",{
  expect_error(scan_to_spatraster(example_scan, param = "unexisting_param"),
               regexp = "'param' contains scan parameter not found in scan",
               fixed = TRUE)
  expect_error(scan_to_spatraster(example_scan, param = "azim"),
               regexp = "'param' should contain the name of one or more scan parameters contained in 'scan'",
               fixed = TRUE)
})

test_that("scan_to_spatraster() raster argument produces expected raster output", {
  data(example_vpts)
  expect_s4_class(b <- scan_to_spatraster(example_scan, ylim = c(55, 57), xlim = c(12, 13), res = .1), "SpatRaster")
  from_template <- scan_to_spatraster(example_scan, raster = b)
  expect_true(terra::compareGeom(b, from_template, stopOnError = FALSE))
  expect_equal(terra::values(b), terra::values(from_template))
})

test_that("scan_to_spatraster() accepts raster templates", {
  skip_if_not_installed("raster")
  template <- raster::raster(scan_to_spatraster(example_scan, ylim = c(55, 57), xlim = c(12, 13), res = .1))
  # a terra SpatRaster template must produce the same output as the equivalent RasterLayer
  from_terra <- scan_to_spatraster(example_scan, raster = terra::rast(template))
  lifecycle::expect_deprecated(from_raster <- scan_to_spatraster(example_scan, raster = template))
  expect_true(terra::compareGeom(from_terra, from_raster, stopOnError = FALSE))
  expect_equal(terra::values(from_terra), terra::values(from_raster))
})

test_that("new scan conversions preserve legacy output", {
  skip_if_not_installed("sp")
  skip_if_not_installed("raster")
  modern_points <- scan_to_sf(example_scan)
  legacy_points <- suppressWarnings(scan_to_spatial(example_scan))

  expect_equal(sf::st_drop_geometry(modern_points), legacy_points@data)
  expect_equal(
    unname(sf::st_coordinates(modern_points)),
    unname(sp::coordinates(legacy_points))
  )
  expect_true(sf::st_crs(modern_points) == sf::st_crs(legacy_points))

  modern_raster <- scan_to_spatraster(example_scan, nx = 20, ny = 20)
  legacy_raster <- suppressWarnings(scan_to_raster(example_scan, nx = 20, ny = 20))
  expect_equal(names(modern_raster), names(legacy_raster))
  expect_equal(terra::values(modern_raster), raster::values(legacy_raster))
  expect_equal(unname(as.vector(terra::ext(modern_raster))), unname(as.vector(raster::extent(legacy_raster))))
  expect_equal(terra::res(modern_raster), raster::res(legacy_raster))
  expect_true(sf::st_crs(modern_raster) == sf::st_crs(legacy_raster))
})

test_that(".scan_to_spdf() returns error error on incorrect parameters",{
  expect_error(bioRad:::.scan_to_spdf("a"),
               regexp = "'scan' should be an object of class scan",
               fixed = TRUE)
  square_birdbath <- example_scan
  square_birdbath$attributes$where$elangle <- 90
  expect_error(bioRad:::.scan_to_spdf(square_birdbath),
               regexp = "georeferencing of 90 degree birdbath scan not supported",
               fixed = TRUE)
})

test_that(".scan_to_spdf() returns error on wrongly formed param argument",{
  points <- scan_to_sf(example_scan)
  expect_error(bioRad:::.scan_to_spdf(example_scan, param = "unexisting_param",
                            spdf = points),
               regexp = "'param' contains scan parameter not found in scan",
               fixed = TRUE)
  expect_error(bioRad:::.scan_to_spdf(example_scan, param = "azim",
                            spdf = points),
               regexp = "'param' should contain the name of one or more scan parameters contained in 'scan'",
               fixed = TRUE)
  expect_error(bioRad:::.scan_to_spdf(example_scan),
               regexp = 'argument "spdf" is missing, with no default',
               fixed = TRUE)
  expect_error(
    bioRad:::.scan_to_spdf(example_scan,
                 points,
                 lat = "a"),
    regexp = "lat is not a number (a length one numeric vector).",
    fixed = TRUE
  )
  expect_error(
    bioRad:::.scan_to_spdf(example_scan,
                 points,
                 lon = "a"),
    regexp = "lon is not a number (a length one numeric vector).",
    fixed = TRUE
  )
  missing_lat_scan <- example_scan
  missing_lon_scan <- example_scan
  missing_lat_scan$geo$lat <- NULL
  missing_lon_scan$geo$lon <- NULL
  expect_error(
    bioRad:::.scan_to_spdf(missing_lat_scan, points),
    regexp = "radar latitude cannot be found in scan, specify using 'lat' argument",
    fixed = TRUE
  )
  expect_error(
    bioRad:::.scan_to_spdf(missing_lon_scan, points),
    regexp = "radar longitude cannot be found in scan, specify using 'lon' argument",
    fixed = TRUE
  )
  expect_error(
    bioRad:::.scan_to_spdf(example_scan, points, k = "a"),
    regexp = "k is not a number (a length one numeric vector).",
    fixed = TRUE
  )
  expect_error(
    bioRad:::.scan_to_spdf(example_scan, points, re = "a"),
    regexp = "re is not a number (a length one numeric vector).",
    fixed = TRUE
  )
  expect_error(
    bioRad:::.scan_to_spdf(example_scan, points, rp = "a"),
    regexp = "rp is not a number (a length one numeric vector).",
    fixed = TRUE
  )
  # TODO add coverage for crs mismatch between local crs and spdf crs
})

test_that("scan_to_spatraster() preserves raster templates", {
  template <- terra::rast(nrows = 4, ncols = 5, xmin = 12.8, xmax = 12.9,
    ymin = 56.3, ymax = 56.4, crs = "EPSG:4326")
  terra::values(template) <- rep(c(NA_real_, 100), 10)
  original <- terra::values(template)
  scan <- scan_to_spatraster(example_scan, raster = template)
  expect_true(terra::compareGeom(scan, template, lyrs = FALSE))
  expect_equal(terra::values(template), original)
  expect_equal(terra::ncell(scan), 20)
  expect_false(all(is.na(terra::values(scan))))
})

test_that("scan_to_spatraster() warns for legacy CRS arguments", {
  skip_if_not_installed("sp")
  lifecycle::expect_deprecated(
    scan_to_spatraster(example_scan, nx = 3, ny = 3, crs = sp::CRS("EPSG:4326"))
  )
})

test_that("scan_to_spatraster() preserves values and grid from upstream master", {
  # Generated with adokter/bioRad d3e03b48 using the bundled example data.
  reference <- readRDS(test_path("fixtures", "spatial-master-d3e03b48.rds"))[["scan"]]
  r <- scan_to_spatraster(example_scan, nx = 30, ny = 25)
  values <- terra::values(r)
  expected <- reference$values
  # The old projection round trip chose either adjacent azimuth bin
  # at exact bin boundaries, depending on the PROJ version/platform.
  # Compare historical values only where that choice is unambiguous.
  geo <- example_scan$geo
  local_crs <- paste0("+proj=aeqd +lat_0=", geo$lat,
    " +lon_0=", geo$lon, " +units=m")
  xy <- sf::st_coordinates(bioRad:::.raster_points(r, local_crs))
  azimuth <- (90 - atan2(xy[, 2], xy[, 1]) * 180 / pi) %% 360
  bins <- azimuth / example_scan$geo$ascale
  stable <- abs(bins - round(bins)) > 1e-8
  values <- values[stable, , drop = FALSE]
  expected <- expected[stable, , drop = FALSE]
  expect_equal(values, expected, tolerance = 1e-8)
  expect_equal(unname(as.vector(terra::ext(r))), reference$extent, tolerance = 1e-8)
  expect_equal(unname(terra::res(r)), reference$res, tolerance = 1e-8)
  expect_equal(names(r), reference$names)
})
