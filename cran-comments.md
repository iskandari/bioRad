## bioRad 1.0.0

* Migrate core spatial processing to sf and terra. New PPIs contain SpatRaster data.
* Retain deprecated legacy scan conversions and spatial input adapters for compatibility.
* Move sp and raster to Suggests; retain optional ggspatial/prettymapr/rosm mapping without legacy raster reprojection.
* Add a migration vignette and numerical regression tests against the previous implementation.
* Rebuild documentation and release metadata; update obsolete documentation URLs.

## Validation

* Local macOS Sequoia 15.7.5, arm64, R 4.5.0: R CMD check --as-cran completed
  with 0 errors, 0 warnings and 0 notes. This includes the extended examples,
  tests, vignette rebuilding, and PDF/HTML manuals with all suggested packages.
* The full local development test suite passes. Spatial regression tests
  compare values and grids against upstream master d3e03b48, and cover legacy
  adapters, geographic DEMs, missing cells, compositing and map rendering.
* The README, roxygen documentation, migration guide and pkgdown site were rebuilt.
  CITATION.cff was regenerated and validated, and codemeta.json regenerated.
* GitHub Actions and R-hub validation is tracked on the fork's
  feature/sf-terra-migration branch. Verify green runs for the final commit
  before submission.

This is release preparation; no CRAN submission has been made.
