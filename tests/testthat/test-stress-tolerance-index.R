# Yp = 10, 20, 30 and Ys = 8, 10, 27, so the trial means are Yp = 20 and
# Ys = 15 and every index can be checked by hand.
sti_trial <- function() {
  data.frame(
    genotype = rep(c("A", "B", "C"), each = 4),
    trt = rep(c("ctrl", "ctrl", "stress", "stress"), 3),
    yield = c(10, 10, 8, 8, 20, 20, 10, 10, 30, 30, 27, 27),
    stringsAsFactors = FALSE
  )
}

sti_value <- function(x, unit, index) x[[index]][x$unit == unit]

test_that("every index matches its published definition", {
  sti <- stress_tolerance_index(sti_trial(), "yield", "trt", "ctrl", "genotype")

  expect_equal(sti_value(sti, "A", "Yp"), 10)
  expect_equal(sti_value(sti, "A", "Ys"), 8)
  expect_equal(unique(sti$stress_intensity), 0.25)

  expect_equal(sti_value(sti, "B", "TOL"), 10)
  expect_equal(sti_value(sti, "C", "MP"), 28.5)
  expect_equal(sti_value(sti, "A", "GMP"), sqrt(80))
  expect_equal(sti_value(sti, "B", "HM"), 2 * 20 * 10 / 30)
  expect_equal(sti_value(sti, "A", "SSI"), 0.8)
  expect_equal(sti_value(sti, "C", "STI"), 30 * 27 / 400)
  expect_equal(sti_value(sti, "A", "YI"), 8 / 15)
  expect_equal(sti_value(sti, "B", "YSI"), 0.5)
  expect_equal(sti_value(sti, "C", "RDI"), 0.9 / 0.75)
  expect_equal(sti_value(sti, "B", "SSPI"), 25)
})

test_that("ranks orient every index towards tolerance", {
  sti <- stress_tolerance_index(sti_trial(), "yield", "trt", "ctrl", "genotype")

  # C keeps 90% of its potential from the highest potential of the trial, so it
  # must win; B loses half of the largest absolute amount and must lose.
  expect_equal(sti$unit[sti$rank_overall == 1], "C")
  expect_equal(sti$unit[sti$rank_overall == max(sti$rank_overall)], "B")
  expect_true(all(sti$rank_mean >= 1 & sti$rank_mean <= 3))
})

test_that("a subset of indices is honoured and changes the ranking basis", {
  sti <- stress_tolerance_index(sti_trial(), "yield", "trt", "ctrl", "genotype",
    indices = c("STI", "GMP")
  )

  expect_true(all(c("STI", "GMP") %in% names(sti)))
  expect_false("TOL" %in% names(sti))
  # Both retained indices reward potential, so the ranking follows them exactly.
  expect_equal(sti$unit[sti$rank_overall == 1], "C")
})

test_that("results are ordered and returned for every stress level", {
  dat <- rbind(
    sti_trial(),
    transform(sti_trial(), trt = ifelse(trt == "stress", "severe", "ctrl"))
  )

  sti <- stress_tolerance_index(dat, "yield", "trt", "ctrl", "genotype")

  expect_setequal(unique(sti$group), c("stress", "severe"))
  expect_equal(nrow(sti), 6L)
  for (g in unique(sti$group)) {
    expect_false(is.unsorted(sti$rank_overall[sti$group == g]))
  }
})

test_that("fun controls how the trait is summarised", {
  dat <- sti_trial()
  dat$yield[dat$genotype == "A" & dat$trt == "ctrl"] <- c(2, 18)

  by_mean <- stress_tolerance_index(dat, "yield", "trt", "ctrl", "genotype")
  by_median <- stress_tolerance_index(dat, "yield", "trt", "ctrl", "genotype",
    fun = stats::median
  )

  expect_equal(sti_value(by_mean, "A", "Yp"), 10)
  expect_equal(sti_value(by_median, "A", "Yp"), 10)

  dat$yield[dat$genotype == "A" & dat$trt == "ctrl"] <- c(2, 22)
  shifted <- stress_tolerance_index(dat, "yield", "trt", "ctrl", "genotype")
  expect_equal(sti_value(shifted, "A", "Yp"), 12)
})

test_that("missing values in the trait are ignored", {
  dat <- sti_trial()
  dat$yield[1] <- NA_real_

  sti <- stress_tolerance_index(dat, "yield", "trt", "ctrl", "genotype")

  expect_equal(sti_value(sti, "A", "Yp"), 10)
})

test_that("input is validated", {
  dat <- sti_trial()

  expect_error(
    stress_tolerance_index(dat, "yield", "trt", "nope", "genotype"),
    class = "plantstressR_error"
  )
  expect_error(
    stress_tolerance_index(dat, "genotype", "trt", "ctrl", "genotype"),
    class = "plantstressR_error"
  )
  expect_error(
    stress_tolerance_index(dat, "yield", "trt", "ctrl", "genotype", indices = "NOPE"),
    class = "plantstressR_error"
  )
  expect_error(
    stress_tolerance_index(dat[dat$genotype == "A", ], "yield", "trt", "ctrl", "genotype"),
    class = "plantstressR_error"
  )

  no_control <- dat
  no_control$trt[no_control$genotype == "B"] <- "stress"
  expect_error(
    stress_tolerance_index(no_control, "yield", "trt", "ctrl", "genotype"),
    class = "plantstressR_error"
  )
})

test_that("non-positive control performance is flagged", {
  dat <- sti_trial()
  dat$yield[dat$genotype == "A" & dat$trt == "ctrl"] <- 0

  expect_warning(
    stress_tolerance_index(dat, "yield", "trt", "ctrl", "genotype"),
    class = "plantstressR_warning"
  )
})

test_that("the shipped trial gives a usable selection table", {
  sti <- stress_tolerance_index(
    brachiaria_stress,
    trait = "shoot_biomass", treatment = "drought_level",
    control = "control", by = "genotype"
  )

  expect_s3_class(sti, "plantstress_sti")
  expect_equal(nrow(sti), 8L)
  expect_true(all(sti$stress_intensity > 0))
  expect_true(all(sti$Yp > sti$Ys))
  expect_output(print(sti), "plantstress_sti")
  expect_output(print(sti), "Stress intensity")
})
