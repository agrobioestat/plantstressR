test_that("the dictionary orients well-known traits", {
  d <- trait_directions(c("Fv_Fm", "DI0_RC", "Ci", "Fo", "shoot_biomass"))

  expect_equal(unname(d[["Fv_Fm"]]), "higher_is_better")
  expect_equal(unname(d[["shoot_biomass"]]), "higher_is_better")
  expect_equal(unname(d[["DI0_RC"]]), "lower_is_better")
  expect_equal(unname(d[["Ci"]]), "lower_is_better")
  expect_equal(unname(d[["Fo"]]), "lower_is_better")
})

test_that("dictionary matching ignores separators and case", {
  expect_equal(
    unname(trait_directions("dio.rc")[["dio.rc"]]),
    "lower_is_better"
  )
})

test_that("explicit arguments override the dictionary", {
  d <- trait_directions(c("Ci", "WUE"), higher_is_better = "Ci", lower_is_better = "WUE")

  expect_equal(unname(d[["Ci"]]), "higher_is_better")
  expect_equal(unname(d[["WUE"]]), "lower_is_better")
})

test_that("use_dictionary = FALSE falls back to the default", {
  d <- trait_directions(c("Ci", "Fo"), use_dictionary = FALSE)

  expect_true(all(d == "higher_is_better"))
})

test_that("conflicting and unknown overrides are reported", {
  expect_error(
    trait_directions("A", lower_is_better = "A", higher_is_better = "A"),
    class = "plantstressR_error"
  )
  expect_warning(
    trait_directions("A", lower_is_better = "missing_trait"),
    class = "plantstressR_warning"
  )
})

test_that("numeric directions are accepted by calculate_sri()", {
  dat <- toy_trial()

  numeric_dir <- calculate_sri(dat, "trt", "ctrl",
    traits = c("up", "down"),
    direction = c(up = 1, down = -1), verbose = FALSE
  )
  char_dir <- calculate_sri(dat, "trt", "ctrl",
    traits = c("up", "down"),
    direction = c(up = "lower_is_better", down = "higher_is_better"),
    verbose = FALSE
  )

  expect_equal(numeric_dir$sri, char_dir$sri)
})
