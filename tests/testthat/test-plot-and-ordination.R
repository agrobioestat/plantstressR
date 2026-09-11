signature_fixture <- function() {
  calculate_sri(
    brachiaria_stress,
    treatment = "drought_level",
    control = "control",
    traits = c("Fv_Fm", "PIabs", "DI0_RC", "A", "gs", "RWC", "SPAD"),
    verbose = FALSE
  )
}

test_that("both signature views build a ggplot", {
  sri <- signature_fixture()

  expect_s3_class(plot_stress_signature(sri), "ggplot")
  expect_s3_class(plot_stress_signature(sri, type = "radar"), "ggplot")
})

test_that("trait filtering options narrow the panel", {
  sri <- signature_fixture()

  p <- plot_stress_signature(sri, top_n = 3)
  expect_equal(length(unique(as.character(p$data$trait))), 3L)

  p2 <- plot_stress_signature(sri, traits = c("Fv_Fm", "A", "RWC"))
  expect_setequal(as.character(unique(p2$data$trait)), c("Fv_Fm", "A", "RWC"))

  p3 <- plot_stress_signature(sri, groups = "severe")
  expect_equal(as.character(unique(p3$data$group)), "severe")
})

test_that("clustering only reorders the traits", {
  sri <- signature_fixture()

  clustered <- levels(plot_stress_signature(sri, cluster = TRUE)$data$trait)
  ranked <- levels(plot_stress_signature(sri, cluster = FALSE)$data$trait)

  expect_setequal(clustered, ranked)
})

test_that("the radar view needs at least three traits", {
  sri <- signature_fixture()

  expect_error(
    plot_stress_signature(sri, type = "radar", top_n = 2),
    class = "plantstressR_error"
  )
})

test_that("plot input is validated", {
  sri <- signature_fixture()

  expect_error(plot_stress_signature(data.frame(a = 1)), class = "plantstressR_error")
  expect_error(plot_stress_signature(sri, traits = "absent"), class = "plantstressR_error")
  expect_error(plot_stress_signature(sri, palette = "red"), class = "plantstressR_error")
  expect_error(plot_stress_signature(sri, top_n = 0), class = "plantstressR_error")
})

test_that("ordination decomposes the trait space", {
  ord <- stress_ordination(
    brachiaria_stress,
    treatment = "drought_level",
    traits = c("Fv_Fm", "PIabs", "A", "gs", "RWC", "SPAD", "shoot_biomass")
  )

  expect_s3_class(ord, "plantstress_ordination")
  expect_equal(nrow(ord$scores), nrow(brachiaria_stress))
  expect_equal(nrow(ord$loadings), 7L)
  expect_true(all(diff(ord$eigenvalues$eigenvalue) <= 0))
  expect_lte(max(ord$eigenvalues$cumulative_pct), 100)
  expect_s3_class(plot(ord), "ggplot")
  expect_output(print(ord), "plantstress_ordination")
})

test_that("the FactoMineR engine reproduces the stats decomposition", {
  skip_if_not_installed("FactoMineR")
  traits <- c("Fv_Fm", "PIabs", "A", "gs", "RWC", "SPAD")

  base <- stress_ordination(brachiaria_stress, traits = traits, engine = "stats")
  fm <- stress_ordination(brachiaria_stress, traits = traits, engine = "FactoMineR")

  expect_equal(fm$engine, "FactoMineR")
  expect_equal(base$eigenvalues$variance_pct, fm$eigenvalues$variance_pct,
    tolerance = 1e-6
  )
  expect_s3_class(plot(fm), "ggplot")
})

test_that("the FactoMineR engine falls back cleanly when absent", {
  skip_if(requireNamespace("FactoMineR", quietly = TRUE), "FactoMineR is installed")

  expect_warning(
    ord <- stress_ordination(brachiaria_stress,
      traits = c("Fv_Fm", "A", "gs", "RWC"), engine = "FactoMineR"
    ),
    class = "plantstressR_warning"
  )

  expect_equal(ord$engine, "stats")
})

test_that("ordination input is validated", {
  expect_error(
    stress_ordination(brachiaria_stress, traits = "Fv_Fm"),
    class = "plantstressR_error"
  )
  expect_error(
    stress_ordination(brachiaria_stress, traits = c("Fv_Fm", "A"), ncomp = 1),
    class = "plantstressR_error"
  )

  ord <- stress_ordination(brachiaria_stress, traits = c("Fv_Fm", "A", "gs"))
  expect_error(plot(ord, components = c(1, 5)), class = "plantstressR_error")
})

test_that("plot() on an SRI object draws the signature", {
  sri <- signature_fixture()

  expect_s3_class(plot(sri), "ggplot")
  expect_s3_class(plot(sri, type = "radar", top_n = 5), "ggplot")
})

test_that("both index views build a ggplot", {
  sri <- calculate_sri(brachiaria_stress,
    treatment = "drought_level", control = "control",
    traits = c("Fv_Fm", "PIabs", "A", "gs", "RWC", "shoot_biomass"),
    by = "genotype", verbose = FALSE
  )
  isi <- integrated_stress_index(sri, weights = "precision")

  expect_s3_class(plot(isi), "ggplot")
  expect_s3_class(plot(isi, type = "contribution"), "ggplot")

  pooled <- plot(isi, type = "contribution", top_n = 3)
  expect_true("other" %in% pooled$data$trait)
  expect_equal(length(unique(pooled$data$trait)), 4L)
})

test_that("index plot input is validated", {
  sri <- calculate_sri(brachiaria_stress,
    treatment = "drought_level", control = "control",
    traits = c("Fv_Fm", "A", "RWC"), by = "genotype", verbose = FALSE
  )
  isi <- integrated_stress_index(sri)

  expect_error(plot(isi, palette = "red"), class = "plantstressR_error")
  expect_error(plot(isi, type = "contribution", top_n = 0), class = "plantstressR_error")
})
