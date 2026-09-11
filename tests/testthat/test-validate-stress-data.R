test_that("a clean design reports no issue", {
  chk <- validate_stress_data(toy_trial(), "trt", "ctrl",
    traits = c("up", "down"), by = "genotype"
  )

  expect_s3_class(chk, "plantstress_validation")
  expect_true(chk$ok)
  expect_equal(nrow(chk$issues), 0L)
  expect_equal(chk$stress_levels, "stress")
  expect_output(print(chk), "plantstress_validation")
})

test_that("design columns are never treated as traits", {
  dat <- toy_trial()
  dat$block <- rep(1:4, length.out = nrow(dat))

  chk <- validate_stress_data(dat, "trt", "ctrl", by = "genotype", verbose = FALSE)

  expect_false("genotype" %in% chk$traits)
  expect_true("block" %in% chk$traits)
})

test_that("hard design problems raise errors", {
  dat <- toy_trial()

  expect_error(validate_stress_data(dat, "trt", "nope"), class = "plantstressR_error")
  expect_error(
    validate_stress_data(dat[dat$trt == "ctrl", ], "trt", "ctrl"),
    class = "plantstressR_error"
  )
  expect_error(validate_stress_data(dat[0, ], "trt", "ctrl"), class = "plantstressR_error")

  no_control <- dat
  no_control$trt[no_control$genotype == "B"] <- "stress"
  expect_error(
    validate_stress_data(no_control, "trt", "ctrl", by = "genotype"),
    class = "plantstressR_error"
  )

  all_na <- dat
  all_na$up <- NA_real_
  expect_error(
    validate_stress_data(all_na, "trt", "ctrl", traits = "up"),
    class = "plantstressR_error"
  )
})

test_that("soft problems are collected and warned about", {
  dat <- toy_trial(n = 2)

  # Every issue is warned about, not only the first one.
  warned <- capture_warnings(
    chk <- validate_stress_data(dat, "trt", "ctrl",
      traits = c("up", "flat"), by = "genotype"
    )
  )

  expect_length(warned, 2L)
  expect_false(chk$ok)
  expect_true("replication" %in% chk$issues$check)
  expect_true("constant_trait" %in% chk$issues$check)
})

test_that("missingness above the threshold is flagged", {
  dat <- toy_trial()
  dat$up[1:12] <- NA_real_

  expect_warning(
    chk <- validate_stress_data(dat, "trt", "ctrl", traits = c("up", "down"), max_missing = 0.1),
    class = "plantstressR_warning"
  )

  expect_true("missing_values" %in% chk$issues$check)
})

test_that("verbose = FALSE keeps the issues without warning", {
  dat <- toy_trial(n = 2)

  expect_silent(
    chk <- validate_stress_data(dat, "trt", "ctrl",
      traits = "up", by = "genotype", verbose = FALSE
    )
  )
  expect_gt(nrow(chk$issues), 0L)
})

test_that("the design table counts every cell", {
  chk <- validate_stress_data(toy_trial(), "trt", "ctrl",
    traits = "up", by = "genotype"
  )

  expect_equal(nrow(chk$design), 4L)
  expect_true(all(chk$design$n == 8L))
})

test_that("the shipped data set is a valid trial", {
  chk <- validate_stress_data(
    brachiaria_stress,
    treatment = "drought_level",
    control = "control",
    traits = c("Fv_Fm", "A", "RWC", "shoot_biomass"),
    by = "genotype",
    verbose = FALSE
  )

  expect_setequal(chk$stress_levels, c("moderate", "severe"))
  expect_equal(chk$n_rows, 120L)
})

test_that("the simulator is reproducible and restores the RNG", {
  a <- simulate_brachiaria_stress(n = 40, seed = 7)
  b <- simulate_brachiaria_stress(n = 40, seed = 7)
  expect_equal(a, b)

  set.seed(42)
  before <- .Random.seed
  invisible(simulate_brachiaria_stress(n = 40, seed = 7))
  expect_identical(.Random.seed, before)

  expect_error(simulate_brachiaria_stress(n = 5), class = "rlang_error")
})
