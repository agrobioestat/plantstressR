# The 0.1.0 dashboard died because nothing tied it to the API it called: the
# 0.2.0 rewrite removed eleven functions and the app went on referring to all of
# them. These tests are that missing tie.

app_source <- function() {
  installed <- system.file("shiny", "app.R", package = "plantstressR")
  if (nzchar(installed) && file.exists(installed)) {
    return(installed)
  }
  # Running from the source tree, before the package is installed.
  local <- testthat::test_path("..", "..", "inst", "shiny", "app.R")
  if (file.exists(local)) local else NA_character_
}

test_that("the app is shipped with the package", {
  path <- app_source()
  skip_if(is.na(path), "The bundled app was not found in this layout.")
  expect_true(file.exists(path))
})

test_that("every plantstressR function the app calls is still exported", {
  path <- app_source()
  skip_if(is.na(path), "The bundled app was not found in this layout.")

  code <- readLines(path, warn = FALSE)
  calls <- regmatches(code, gregexpr("plantstressR::[A-Za-z_.][A-Za-z0-9_.]*", code))
  used <- unique(sub("^plantstressR::", "", unlist(calls)))

  expect_gt(length(used), 5L)
  exported <- getNamespaceExports("plantstressR")
  expect_setequal(setdiff(used, exported), character())
})

test_that("the app parses as valid R code", {
  path <- app_source()
  skip_if(is.na(path), "The bundled app was not found in this layout.")
  expect_silent(parse(path))
})

test_that("run_plantstress_app fails with a clear message without shiny", {
  skip_if(requireNamespace("shiny", quietly = TRUE), "shiny is installed")
  expect_error(run_plantstress_app(), class = "plantstressR_error")
})
