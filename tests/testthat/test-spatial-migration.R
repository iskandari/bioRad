test_that("legacy PPI conversion preserves data and supports downstream methods", {
  skip_if_not_installed("raster")
  skip_if_not_installed("sp")
  modern <- project_as_ppi(example_scan, grid_size = 2000, range_max = 10000)
  legacy <- modern
  legacy$data <- methods::as(raster::brick(modern$data), "SpatialGridDataFrame")
  lifecycle::expect_deprecated(converted <- convert_legacy(legacy))
  expect_true(terra::compareGeom(converted$data, modern$data))
  expect_equal(terra::values(converted$data), terra::values(modern$data))
  expect_equal(converted$geo, modern$geo)
  expect_s4_class(legacy$data, "SpatialGridDataFrame")
  lifecycle::expect_deprecated(subset <- legacy[1:2])
  expect_equal(names(subset$data), names(modern$data)[1:2])
  lifecycle::expect_deprecated(calculated <- calculate_param(legacy, doubled = DBZH * 2))
  expect_equal(terra::values(calculated$data$doubled), terra::values(modern$data$DBZH) * 2,
    ignore_attr = TRUE)
})

test_that("templates retain missing cells and do not change the caller's raster", {
  template <- terra::rast(nrows = 4, ncols = 5, xmin = 12.8, xmax = 12.9,
    ymin = 56.3, ymax = 56.4, crs = "EPSG:4326")
  terra::values(template) <- rep(c(NA_real_, 100), 10)
  original <- terra::values(template)
  ppi <- project_as_ppi(example_scan, raster = template)
  scan <- scan_to_spatraster(example_scan, raster = template)
  expect_true(terra::compareGeom(ppi$data, template, lyrs = FALSE))
  expect_true(terra::compareGeom(scan, template, lyrs = FALSE))
  expect_equal(terra::values(template), original)
  expect_equal(terra::ncell(ppi$data), 20)
  expect_false(all(is.na(terra::values(ppi$data))))
  param <- project_as_ppi(example_scan$params[[1]], raster = template)
  expect_equal(terra::values(param$data), terra::values(ppi$data[[1]]))
})

test_that("geographic DEM extraction samples points in the correct CRS", {
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

test_that("composites preserve missing values, coverage and weighted means", {
  r <- terra::rast(nrows = 3, ncols = 3, xmin = -1500, xmax = 1500,
    ymin = -1500, ymax = 1500, crs = "+proj=aeqd +lat_0=56 +lon_0=13 +datum=WGS84")
  terra::values(r) <- c(NA, rep(2, 8))
  names(r) <- "DBZH"
  a <- structure(list(data = r, geo = list(lon = 13, lat = 56, elangle = .5,
    bbox = bioRad:::.raster_bbox(r))), class = "ppi")
  b <- a
  b$data <- r * 2
  for (method in c("min", "max", "mean", "idw")) {
    out <- composite_ppi(list(a, b), raster = r, method = method, coverage = TRUE)
    expected <- switch(method, min = 2, max = 4, mean = 3, idw = 3)
    expect_equal(as.vector(terra::values(out$data$DBZH)), c(NA, rep(expected, 8)))
    expect_equal(as.vector(terra::values(out$data$coverage)), rep(2, 9))
  }
  one <- composite_ppi(list(a), raster = r, method = "mean")
  expect_equal(terra::values(one$data), terra::values(r))
  expect_equal(terra::values(a$data), terra::values(r))
})

test_that("legacy CRS and raster arguments have migration warnings", {
  skip_if_not_installed("sp")
  skip_if_not_installed("raster")
  template <- terra::rast(nrows=3, ncols=3, xmin=12.8, xmax=12.9,
    ymin=56.3, ymax=56.4, crs="EPSG:4326")
  terra::values(template) <- 100
  lifecycle::expect_deprecated(add_param(example_scan, raster::raster(template), "HGHT"))
  lifecycle::expect_deprecated(scan_to_spatraster(example_scan, nx=3, ny=3, crs=sp::CRS("EPSG:4326")))
})

test_that("spatial values and grids match upstream master before migration", {
  # Generated with adokter/bioRad d3e03b48 using the bundled example data.
  reference <- readRDS(test_path("fixtures", "spatial-master-d3e03b48.rds"))
  pvol <- read_pvolfile(system.file("extdata", "volume.h5", package = "bioRad"))
  current <- list(
    ppi = project_as_ppi(example_scan, grid_size = 2000, range_max = 20000)$data,
    scan = scan_to_spatraster(example_scan, nx = 30, ny = 25),
    integrated = integrate_to_ppi(pvol, example_vp, nx = 20, ny = 20)$data
  )
  for (name in names(current)) {
    r <- current[[name]]
    values <- terra::values(r)
    expected <- reference[[name]]$values
    {
      # The old projection round trip chose either adjacent azimuth bin
      # at exact bin boundaries, depending on the PROJ version/platform.
      # Compare historical values only where that choice is unambiguous.
      geo <- if (name == "integrated") pvol$geo else example_scan$geo
      local_crs <- paste0("+proj=aeqd +lat_0=", geo$lat,
        " +lon_0=", geo$lon, " +units=m")
      xy <- sf::st_coordinates(bioRad:::.raster_points(r, local_crs))
      azimuth <- (90 - atan2(xy[, 2], xy[, 1]) * 180 / pi) %% 360
      bins <- azimuth / example_scan$geo$ascale
      stable <- abs(bins - round(bins)) > 1e-8
      values <- values[stable, , drop = FALSE]
      expected <- expected[stable, , drop = FALSE]
    }
    expect_equal(values, expected, tolerance = 1e-8)
    expect_equal(unname(as.vector(terra::ext(r))), reference[[name]]$extent, tolerance = 1e-8)
    expect_equal(unname(terra::res(r)), reference[[name]]$res, tolerance = 1e-8)
    expect_equal(names(r), reference[[name]]$names)
  }
})

test_that("legacy elevation rasters preserve polar-volume extraction", {
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


test_that("projection round-off does not change requested PPI grids", {
  original <- project_as_ppi(example_scan, grid_size = 2000, range_max = 20000)
  transform <- bioRad:::wgs_to_proj
  local_mocked_bindings(wgs_to_proj = function(...) {
    points <- transform(...)
    sf::st_geometry(points) <- sf::st_geometry(points) + c(-1e-9, 1e-9)
    points
  }, .package = "bioRad")
  perturbed <- project_as_ppi(example_scan, grid_size = 2000, range_max = 20000)
  expect_equal(dim(perturbed$data), c(20L, 20L, 5L))
  expect_equal(as.vector(terra::ext(perturbed$data)),
    as.vector(terra::ext(original$data)))
  expect_equal(terra::values(perturbed$data), terra::values(original$data))
})


test_that("azimuth-bin boundaries are stable under projection round-off", {
  angles <- c(0, 45, 225, 359, 360)
  sample <- function(delta) bioRad:::polar_to_index(
    data.frame(range = 1000, azim = angles + delta), azimbin = 1)$col
  expect_equal(sample(0), c(1, 46, 226, 360, 1))
  expect_equal(sample(-1e-12), sample(0))
  expect_equal(sample(1e-12), sample(0))
})
