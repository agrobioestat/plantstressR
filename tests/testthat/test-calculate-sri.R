test_that("the sign convention makes positive SRI mean damage", {
  dat <- toy_trial()

  sri <- calculate_sri(dat, "trt", "ctrl",
    traits = c("up", "down"),
    direction = c(up = "lower_is_better", down = "higher_is_better"),
    verbose = FALSE
  )

  # `up` rises under stress and is bad when high; `down` falls and is good when
  # high. Both must therefore come out positive.
  expect_true(all(sri$sri > 0))
  expect_true(all(sri$difference[sri$trait == "up"] > 0))
  expect_true(all(sri$difference[sri$trait == "down"] < 0))
})

test_that("flipping the declared direction flips the index", {
  dat <- toy_trial()

  a <- calculate_sri(dat, "trt", "ctrl",
    traits = "up",
    direction = c(up = "higher_is_better"), verbose = FALSE
  )
  b <- calculate_sri(dat, "trt", "ctrl",
    traits = "up",
    direction = c(up = "lower_is_better"), verbose = FALSE
  )

  expect_equal(a$sri, -b$sri)
  expect_equal(a$conf_low, -b$conf_high)
  expect_equal(a$conf_high, -b$conf_low)
})

test_that("Glass's delta matches its definition", {
  dat <- toy_trial()
  ctrl <- dat$up[dat$genotype == "A" & dat$trt == "ctrl"]
  stress <- dat$up[dat$genotype == "A" & dat$trt == "stress"]
  expected <- (mean(stress) - mean(ctrl)) / stats::sd(ctrl)

  sri <- calculate_sri(dat, "trt", "ctrl",
    traits = "up", by = "genotype",
    direction = c(up = "lower_is_better"), method = "glass", verbose = FALSE
  )

  expect_equal(sri$sri[sri$unit == "A"], expected)
})

test_that("the four methods agree on sign and order of magnitude", {
  dat <- toy_trial()
  methods <- c("glass", "cohen", "hedges", "relative")

  values <- vapply(methods, function(m) {
    calculate_sri(dat, "trt", "ctrl",
      traits = "up", direction = c(up = "lower_is_better"),
      method = m, verbose = FALSE
    )$sri
  }, numeric(1))

  expect_true(all(values > 0))
  expect_lt(values[["relative"]], values[["glass"]])
  expect_lt(abs(values[["hedges"]]), abs(values[["cohen"]]))
})

test_that("`by` computes each unit against its own control", {
  dat <- toy_trial()

  sri <- calculate_sri(dat, "trt", "ctrl",
    traits = c("up", "down"), by = "genotype", verbose = FALSE
  )

  expect_setequal(unique(sri$unit), c("A", "B"))
  expect_equal(nrow(sri), 4L)

  # Genotype A was pushed twice as hard as B on both traits.
  a <- sri[sri$unit == "A", ]
  b <- sri[sri$unit == "B", ]
  expect_true(all(abs(a$sri[order(a$trait)]) > abs(b$sri[order(b$trait)])))
})

test_that("a constant trait yields NA rather than an error", {
  dat <- toy_trial()

  expect_warning(
    sri <- calculate_sri(dat, "trt", "ctrl", traits = c("up", "flat")),
    class = "plantstressR_warning"
  )

  expect_true(is.na(sri$sri[sri$trait == "flat"]))
  expect_false(is.na(sri$sri[sri$trait == "up"]))
})

test_that("multiple stress levels each get their own comparison", {
  dat <- rbind(
    toy_trial(),
    transform(toy_trial(), trt = ifelse(trt == "stress", "severe", "ctrl"))
  )

  sri <- calculate_sri(dat, "trt", "ctrl", traits = "up", verbose = FALSE)

  expect_setequal(unique(sri$group), c("stress", "severe"))
})

test_that("design problems are rejected with classed errors", {
  dat <- toy_trial()

  expect_error(calculate_sri(dat, "trt", "absent"), class = "plantstressR_error")
  expect_error(calculate_sri(dat, "missing", "ctrl"), class = "plantstressR_error")
  expect_error(
    calculate_sri(dat, "trt", "ctrl", traits = "not_a_column"),
    class = "plantstressR_error"
  )
  expect_error(
    calculate_sri(dat, "trt", "ctrl", traits = "genotype"),
    class = "plantstressR_error"
  )
  expect_error(
    calculate_sri(dat, "trt", "ctrl", conf_level = 2),
    class = "plantstressR_error"
  )
})

test_that("print and summary methods work", {
  sri <- calculate_sri(toy_trial(), "trt", "ctrl",
    traits = c("up", "down"), by = "genotype", verbose = FALSE
  )

  expect_output(print(sri), "plantstress_sri")
  s <- summary(sri)
  expect_s3_class(s, "tbl_df")
  expect_true(all(c("mean_sri", "most_impaired") %in% names(s)))
})

test_that("se_sampling reduces to the textbook form under equal variances", {
  # Two groups built to have exactly the same spread, so the heteroscedastic
  # expression must collapse onto sqrt((n_c + n_s) / (n_c * n_s)).
  n <- 15
  spread <- as.numeric(scale(seq_len(n)))
  dat <- data.frame(
    trt = rep(c("ctrl", "stress"), each = n),
    y = c(10 + spread, 14 + spread)
  )

  sri <- calculate_sri(dat,
    treatment = "trt", control = "ctrl", traits = "y",
    method = "glass", verbose = FALSE
  )
  expect_equal(sri$se_sampling, sqrt((n + n) / (n * n)), tolerance = 1e-8)
})

test_that("se_sampling grows when the stress treatment inflates the variance", {
  # Same means, same replication; only the dispersion under stress differs. The
  # index is measured much less precisely in the erratic trait, and the
  # Welch-style numerator variance is what notices -- the homoscedastic form
  # would report both as equally precise while `p_value` used a Welch test.
  set.seed(11)
  n <- 20
  dat <- data.frame(
    trt = rep(c("ctrl", "stress"), each = n),
    steady = c(stats::rnorm(n, 10, 1), stats::rnorm(n, 14, 1)),
    erratic = c(stats::rnorm(n, 10, 1), stats::rnorm(n, 14, 6))
  )

  sri <- calculate_sri(dat,
    treatment = "trt", control = "ctrl",
    traits = c("steady", "erratic"), method = "glass", verbose = FALSE
  )
  steady <- sri$se_sampling[sri$trait == "steady"]
  erratic <- sri$se_sampling[sri$trait == "erratic"]
  expect_gt(erratic, 2 * steady)

  # And precision weighting must act on it.
  w <- stress_weights(integrated_stress_index(sri,
    weights = "precision", rescale = FALSE
  ))
  expect_gt(w$weight[w$trait == "steady"], w$weight[w$trait == "erratic"])
})
