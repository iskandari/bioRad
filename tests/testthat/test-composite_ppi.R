test_that("composite_ppi() returns error on incorrect parameters", {
  pvol <- read_pvolfile(system.file("extdata", "volume.h5", package = "bioRad"))
  ppis <- lapply(pvol$scans, project_as_ppi)
  expect_error(
    composite_ppi(example_vp),
    regexp = "'x' should be an object of class ppi, or a list of objects of class ppi",
    fixed = TRUE
  )
  expect_error(
    composite_ppi(ppis[2:3], nx = "a"),
    regexp = "'nx' should be an integer",
    fixed = TRUE
  )
  expect_error(
    composite_ppi(ppis[2:3], ny = "a"),
    regexp = "'ny' should be an integer",
    fixed = TRUE
  )
  expect_error(
    composite_ppi(ppis[2:3], xlim = "a"),
    regexp = "'xlim' should be a numeric vector of length two",
    fixed = TRUE
  )
  expect_error(
    composite_ppi(ppis[2:3], xlim = c(1,9,4)),
    regexp = "'xlim' should be a numeric vector of length two",
    fixed = TRUE
  )
  expect_error(
    composite_ppi(ppis[2:3], xlim = c(9,NA)),
    regexp = paste("'xlim' should be a vector with two numeric values",
                     "for lower and upper bound respectively"),
    fixed = TRUE
  )
  expect_error(
    composite_ppi(ppis[2:3], ylim = "a"),
    regexp = "'ylim' should be a numeric vector of length two",
    fixed = TRUE
  )
  expect_error(
    composite_ppi(ppis[2:3], ylim = c(1,9,4)),
    regexp = "'ylim' should be a numeric vector of length two",
    fixed = TRUE
  )
  expect_error(
    composite_ppi(ppis[2:3], ylim = c(9,NA)),
    regexp = paste("'ylim' should be a vector with two numeric values",
                   "for lower and upper bound respectively"),
    fixed = TRUE
  )
  expect_error(
    composite_ppi(ppis[1:2], res = seq(3)),
    regexp = "length(res) not less than or equal to 2",
    fixed = TRUE
  )
  expect_error(
    composite_ppi(ppis[1:2], res = c("a",2)),
    regexp = "res is not a numeric or integer vector",
    fixed = TRUE
  )
  expect_error(
    composite_ppi(ppis[1:2], crs = "a"),
    regexp = "invalid crs: a",
    fixed = TRUE
  )
  expect_error(
    composite_ppi(ppis[1:2], method = "not_a_method"),
    regexp = "'method' should be one or multiple of 'max', 'mean', 'min' or 'idw'",
    fixed = TRUE
  )
  expect_error(
    composite_ppi(ppis[1:2], param = "all", method = c("max","min","mean")),
    regexp = "'method' should be of length 1 or length(param)",
    fixed = TRUE
  )
  expect_error(
    composite_ppi(ppis[1:2], param = "DBZH", method = c("max","min")),
    regexp = "'method' should be of length 1 or length(param)",
    fixed = TRUE
  )
  expect_error(
    composite_ppi(ppis[1:2],
                  param = c("DBZH","ZDR"),
                  method = c("max","min","mean")),
    regexp = "'method' should be of length 1 or length(param)",
    fixed = TRUE
  )
  expect_error(
    composite_ppi(ppis[1:2], coverage = 1),
    regexp = "coverage is not a flag (a length one logical vector).",
    fixed = TRUE
  )
  expect_error(
    composite_ppi(ppis[1:2], coverage = c(TRUE,TRUE)),
    regexp = "coverage is not a flag (a length one logical vector).",
    fixed = TRUE
  )
  expect_error(
    composite_ppi(ppis[1:2], coverage = NA),
    regexp = "missing value where TRUE/FALSE needed",
    fixed = TRUE
  )
})

test_that("composite_ppi() returns terra SpatRaster data", {
  pvol <- read_pvolfile(system.file("extdata", "volume.h5", package = "bioRad"))
  ppis <- lapply(pvol$scans[1:2], project_as_ppi)
  ppi <- composite_ppi(ppis, nx = 20, ny = 20, param = "DBZH")

  expect_s3_class(ppi, "ppi")
  expect_s4_class(ppi$data, "SpatRaster")
  expect_equal(names(ppi$data), "DBZH")
})

test_that("composite_ppi() preserves missing values, coverage and weighted means", {
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
