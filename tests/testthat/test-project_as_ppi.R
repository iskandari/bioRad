
remove_bboxlatlon <- function(x) { attributes(x)$bboxlatlon <- NULL; return(x) }

test_that("returns error on incorrect parameters", {
  # I don't think a user could naturally reach these errors:
  expect_error(
    project_as_ppi.param(example_scan),
    regexp = 'inherits(x, "param") is not TRUE',
    fixed = TRUE
  )
  expect_error(
    project_as_ppi.scan(get_param(example_scan, "DBZH")),
    regexp = 'inherits(x, "scan") is not TRUE',
    fixed = TRUE
  )
})

test_that("sample_polar works", {
  data("example_scan")
  a <- example_scan$params[[3]]
  expect_s4_class(b <- bioRad:::sample_polar(example_scan$params[[3]], 500, 10000, xlim = 12.9 + c(-1, 1), ylim = 56.4 + c(-1, 1), project = F), "SpatRaster")
  expect_equal(terra::values(b), terra::values(bioRad:::sample_polar(example_scan$params[[3]], b, project = FALSE)))
  expect_s4_class(bb <- bioRad:::sample_polar(example_scan$params[[3]], 5000, 10000, xlim = 12.9 + c(-1, 1), ylim = 56.4 + c(-1, 1), project = T), "SpatRaster")
  expect_equal(terra::values(bb), terra::values(bioRad:::sample_polar(example_scan$params[[3]], bb, project = TRUE)))
})

test_that("project_as_ppi works", {
  data("example_scan")
  expect_s3_class(b <- project_as_ppi(example_scan, 500, 10000, project = F), "ppi")
  expect_s3_class(bb <- project_as_ppi(example_scan, raster = b$data, project = F), "ppi")
  expect_equal(bb$radar, b$radar, ignore_attr = TRUE)
  expect_equal(bb$datetime, b$datetime, ignore_attr = TRUE)
  expect_true(terra::compareGeom(bb$data, b$data, stopOnError = FALSE))
  expect_equal(terra::values(bb$data), terra::values(b$data), ignore_attr = TRUE)
  expect_s3_class(b <- project_as_ppi(example_scan, 50, 1000, project = T), "ppi")
  expect_s3_class(bb <- project_as_ppi(example_scan, raster = b$data, project = T), "ppi")
  expect_equal(bb$radar, b$radar, ignore_attr = TRUE)
  expect_equal(bb$datetime, b$datetime, ignore_attr = TRUE)
  expect_true(terra::compareGeom(bb$data, b$data, stopOnError = FALSE))
  expect_equal(terra::values(bb$data), terra::values(b$data), ignore_attr = TRUE)
})

test_that("project_as_ppi preserves projected values and geometry", {
  legacy <- bioRad:::sample_polar(
    example_scan$params[[1]], 500, 10000, project = FALSE,
    ylim = NULL, xlim = NULL
  )
  ppi <- project_as_ppi(example_scan$params[[1]], 500, 10000, project = FALSE)

  expect_s4_class(ppi$data, "SpatRaster")
  expect_equal(terra::values(ppi$data)[, 1], terra::values(legacy)[, 1])
  expect_equal(
    unname(as.vector(terra::ext(ppi$data))),
    unname(as.vector(terra::ext(legacy)))
  )
  expect_true(sf::st_crs(ppi$data) == sf::st_crs(legacy))
})

test_that("project_as_ppi() preserves its geographic bounding box", {
  data("example_scan")
  range_max <- 10000
  ppi <- project_as_ppi(example_scan, 500, range_max)
  local_crs <- paste0(
    "+proj=aeqd +lat_0=", example_scan$geo$lat,
    " +lon_0=", example_scan$geo$lon, " +units=m"
  )
  corners <- sf::st_as_sf(data.frame(x=c(-range_max, range_max), y=c(-range_max, range_max)),
    coords=c("x", "y"), crs=local_crs)
  expected <- matrix(sf::st_bbox(sf::st_transform(corners, 4326)), nrow=2,
    dimnames=list(c("lon","lat"),c("min","max")))
  rownames(expected) <- c("lon", "lat")

  expect_equal(ppi$geo$bbox, expected, tolerance = 1e-7)

  limited <- project_as_ppi(
    example_scan, 500, range_max,
    xlim = c(12.8, 13), ylim = c(56.3, 56.5)
  )
  expect_equal(unname(limited$geo$bbox["lon", ]), c(12.8, 13))
  expect_equal(unname(limited$geo$bbox["lat", ]), c(56.3, 56.5))
})


test_that("project_as_ppi() accepts a terra SpatRaster as raster argument", {
  data("example_scan")
  b <- project_as_ppi(example_scan, 500, 10000, project = FALSE)
  skip_if_not_installed("raster")
  template <- raster::raster(b$data)
  # a terra SpatRaster template must produce the same ppi as the equivalent RasterLayer
  ppi_spat <- project_as_ppi(example_scan, raster = terra::rast(template), project = FALSE)
  lifecycle::expect_deprecated(ppi_rast <- project_as_ppi(example_scan, raster = template, project = FALSE))
  expect_s3_class(ppi_spat, "ppi")
  expect_true(terra::compareGeom(ppi_spat$data, ppi_rast$data, stopOnError = FALSE))
  expect_equal(terra::values(ppi_spat$data), terra::values(ppi_rast$data), ignore_attr = TRUE)
})

test_that("project_as_ppi works from different projection", {
  data("example_scan")
  expect_s3_class(b <- project_as_ppi(example_scan, 3000, 5000, project = F), "ppi")
  legacy_grid <- bioRad:::sample_polar(
    example_scan$params[[1]], 3000, 5000, project = F,
    ylim = NULL, xlim = NULL
  )
  skip_if_not_installed("raster")
  skip_if_not_installed("sp")
  legacy_grid <- methods::as(raster::brick(legacy_grid), "SpatialGridDataFrame")
  expect_s4_class(r <- raster::rasterFromXYZ(sp::SpatialPoints(sp::spTransform(as(legacy_grid, "SpatialPoints")[s <- c(1, 11), ], "+proj=longlat"), proj4string= sp::CRS("+proj=longlat")), crs = "+proj=longlat"), "RasterLayer")
  lifecycle::expect_deprecated(bb <- project_as_ppi(example_scan, raster = r, project = FALSE))
  expect_s3_class(bb, "ppi")
  expect_equal(terra::values(b$data)[s, ], terra::values(bb$data)[c(1, 4), ], ignore_attr = TRUE)
})

test_that("project_as_ppi() preserves raster templates", {
  template <- terra::rast(nrows = 4, ncols = 5, xmin = 12.8, xmax = 12.9,
    ymin = 56.3, ymax = 56.4, crs = "EPSG:4326")
  terra::values(template) <- rep(c(NA_real_, 100), 10)
  original <- terra::values(template)
  ppi <- project_as_ppi(example_scan, raster = template)
  expect_true(terra::compareGeom(ppi$data, template, lyrs = FALSE))
  expect_equal(terra::values(template), original)
  expect_equal(terra::ncell(ppi$data), 20)
  expect_false(all(is.na(terra::values(ppi$data))))
  param <- project_as_ppi(example_scan$params[[1]], raster = template)
  expect_equal(terra::values(param$data), terra::values(ppi$data[[1]]))
})

test_that("project_as_ppi() preserves values and grid from upstream master", {
  # Generated with adokter/bioRad d3e03b48 using the bundled example data.
  reference <- readRDS(test_path("fixtures", "spatial-master-d3e03b48.rds"))[["ppi"]]
  r <- project_as_ppi(example_scan, grid_size = 2000, range_max = 20000)$data
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

test_that("project_as_ppi() preserves grids under projection round-off", {
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

test_that("polar_to_index() preserves azimuth bins under projection round-off", {
  angles <- c(0, 45, 225, 359, 360)
  sample <- function(delta) bioRad:::polar_to_index(
    data.frame(range = 1000, azim = angles + delta), azimbin = 1)$col
  expect_equal(sample(0), c(1, 46, 226, 360, 1))
  expect_equal(sample(-1e-12), sample(0))
  expect_equal(sample(1e-12), sample(0))
})
