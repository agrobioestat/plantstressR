test_that("the fast path reproduces the differences of the full computation", {
  # stress_index_ci(scale = "fixed") rebuilds each index from a difference of
  # means instead of calling calculate_sri(). Applied to the observed trial it
  # must land on exactly the numbers calculate_sri() reported.
  dat <- toy_block_trial()
  meta_of <- function(x) attr(x, "plantstress", exact = TRUE)

  for (blk in list(NULL, "block")) {
    sri <- calculate_sri(dat,
      treatment = "trt", control = "ctrl", traits = c("up", "down"),
      by = "genotype", block = blk, verbose = FALSE
    )
    diffs <- boot_differences(dat, meta_of(sri))

    key_obs <- paste(sri$unit, sri$group, sri$trait, sep = "\r")
    key_new <- paste(diffs$unit, diffs$group, diffs$trait, sep = "\r")
    expect_equal(
      diffs$difference[match(key_obs, key_new)],
      sri$difference,
      tolerance = 1e-8
    )
  }
})

test_that("stress_index_ci returns a coherent interval and rank range", {
  data(brachiaria_stress, envir = environment())
  sri <- calculate_sri(brachiaria_stress,
    treatment = "drought_level", control = "control",
    traits = c("Fv_Fm", "A", "RWC", "shoot_biomass"),
    by = "genotype", block = "block", verbose = FALSE
  )

  ci <- stress_index_ci(sri, n_boot = 30, seed = 42)
  expect_s3_class(ci, "plantstress_isi_ci")
  # Note what is *not* asserted: a percentile interval is not required to
  # contain the point estimate, and a skewed bootstrap distribution can legally
  # put `isi` outside its own limits.
  expect_true(all(ci$conf_low <= ci$conf_high))
  expect_true(all(is.finite(ci$conf_low) & is.finite(ci$conf_high)))
  expect_true(all(ci$rank_low <= ci$rank_high))
  expect_true(all(ci$p_best >= 0 & ci$p_best <= 1))
  expect_true(all(ci$n_ok <= 30))

  # Exactly one unit can be ranked first in each resample.
  totals <- tapply(ci$p_best, ci$group, sum)
  expect_true(all(abs(totals - 1) < 1e-8))
})

test_that("the observed index matches integrated_stress_index", {
  data(brachiaria_stress, envir = environment())
  sri <- calculate_sri(brachiaria_stress,
    treatment = "drought_level", control = "control",
    traits = c("Fv_Fm", "A", "RWC"), by = "genotype", verbose = FALSE
  )
  ci <- stress_index_ci(sri, n_boot = 10, seed = 7)
  isi <- integrated_stress_index(sri, rescale = FALSE)

  key <- paste(isi$unit, isi$group)
  expect_equal(ci$isi[match(key, paste(ci$unit, ci$group))], isi$isi,
    tolerance = 1e-10
  )
})

test_that("holding the scale fixed tames the upper tail", {
  # The control standard deviation is itself an estimate; resampling it lets it
  # collapse and sends the index up. The default must not inherit that.
  data(brachiaria_stress, envir = environment())
  sri <- calculate_sri(brachiaria_stress,
    treatment = "drought_level", control = "control",
    traits = c("Fv_Fm", "A", "RWC", "shoot_biomass"),
    by = "genotype", block = "block", verbose = FALSE
  )

  fixed <- stress_index_ci(sri, n_boot = 30, seed = 3, scale = "fixed")
  loose <- stress_index_ci(sri, n_boot = 30, seed = 3, scale = "resampled")

  width <- function(x) mean(x$conf_high - x$conf_low)
  expect_lt(width(fixed), width(loose))
  expect_equal(fixed$isi, loose$isi, tolerance = 1e-10)
})

test_that("the seed makes the interval reproducible and restores the stream", {
  data(brachiaria_stress, envir = environment())
  sri <- calculate_sri(brachiaria_stress,
    treatment = "drought_level", control = "control",
    traits = c("Fv_Fm", "A", "RWC"), by = "genotype", verbose = FALSE
  )

  set.seed(99)
  before <- stats::runif(1)
  set.seed(99)
  a <- stress_index_ci(sri, n_boot = 15, seed = 5)
  after <- stats::runif(1)
  b <- stress_index_ci(sri, n_boot = 15, seed = 5)

  expect_equal(a$conf_low, b$conf_low)
  expect_equal(before, after)
})

test_that("stress_index_ci rejects what it cannot resample", {
  data(brachiaria_stress, envir = environment())
  sri <- calculate_sri(brachiaria_stress,
    treatment = "drought_level", control = "control",
    traits = c("Fv_Fm", "A", "RWC"), by = "genotype", verbose = FALSE
  )

  expect_error(stress_index_ci(sri, n_boot = 1), class = "plantstressR_error")
  expect_error(stress_index_ci(as.data.frame(sri)), class = "plantstressR_error")

  stripped <- sri
  meta <- attr(stripped, "plantstress", exact = TRUE)
  meta$data <- NULL
  attr(stripped, "plantstress") <- meta
  expect_error(stress_index_ci(stripped, n_boot = 5), class = "plantstressR_error")
})

test_that("resampling preserves the layout of the trial", {
  dat <- toy_block_trial()
  set.seed(11)

  # Blocked: whole blocks are drawn, and duplicates are relabelled so that a
  # block taken twice counts as two strata.
  boot <- resample_trial(dat, "trt", "genotype", "block")
  expect_equal(nrow(boot), nrow(dat))
  expect_lte(length(unique(boot$block)), length(unique(dat$block)))

  # Unblocked: replication of every unit-by-treatment cell is preserved.
  boot2 <- resample_trial(dat, "trt", "genotype", NULL)
  expect_equal(
    as.vector(table(boot2$genotype, boot2$trt)),
    as.vector(table(dat$genotype, dat$trt))
  )
})

test_that("the interval has a plot method", {
  skip_if_not_installed("ggplot2")
  data(brachiaria_stress, envir = environment())
  sri <- calculate_sri(brachiaria_stress,
    treatment = "drought_level", control = "control",
    traits = c("Fv_Fm", "A", "RWC"), by = "genotype", verbose = FALSE
  )
  ci <- stress_index_ci(sri, n_boot = 10, seed = 2)

  expect_s3_class(plot(ci), "ggplot")
  expect_s3_class(plot(ci, top_n = 2), "ggplot")
  expect_error(plot(ci, top_n = 0), class = "plantstressR_error")
})
