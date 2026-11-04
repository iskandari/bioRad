# Make legacy objects
legacy_vp <- example_vp
legacy_vpts <- example_vpts
names(legacy_vp$data) <- sub("height", "HGHT", names(legacy_vp$data))
names(legacy_vpts$data) <- sub("height", "heights", names(legacy_vpts$data))
names(legacy_vpts$data) <- sub("dates", "datetime", names(legacy_vpts$data))

test_that("convert_legacy() returns error on incorrect parameters", {
  expect_error(convert_legacy("not_a_object"))
})

test_that("convert_legacy() returns same object type", {
  expect_s3_class(convert_legacy(legacy_vp), "vp")
  expect_s3_class(convert_legacy(legacy_vpts), "vpts")
})

test_that("convert_legacy() renames certain data colnames", {
  expect_identical(convert_legacy(legacy_vp), example_vp)
  expect_identical(convert_legacy(legacy_vpts), example_vpts)
})

test_that("convert_legacy() preserves PPI data and geometry", {
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
})
