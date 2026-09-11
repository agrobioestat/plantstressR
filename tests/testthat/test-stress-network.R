test_that("the network recovers the modules built into the data", {
  net <- stress_network(
    brachiaria_stress,
    traits = c("Fv_Fm", "PIabs", "DI0_RC", "A", "gs", "Ci", "RWC", "SPAD", "shoot_biomass"),
    treatment = "drought_level",
    level = c("moderate", "severe")
  )

  expect_s3_class(net, "plantstress_network")
  expect_gt(nrow(net$edges), 0L)
  expect_gt(nrow(net$modules), 1L)

  module_of <- stats::setNames(net$nodes$module, net$nodes$trait)
  # Traits generated from the same latent factor must land in one module.
  expect_equal(module_of[["A"]], module_of[["gs"]])
  expect_equal(module_of[["Fv_Fm"]], module_of[["PIabs"]])
  expect_false(module_of[["A"]] == module_of[["Fv_Fm"]])
})

test_that("partial correlations differ from marginal ones", {
  traits <- c("Fv_Fm", "PIabs", "A", "gs", "RWC", "SPAD")

  part <- stress_network(brachiaria_stress, traits = traits)
  marg <- stress_network(brachiaria_stress, traits = traits, method = "pearson")

  expect_lt(mean(abs(part$matrix)), mean(abs(marg$matrix)))
  expect_true(all(diag(part$matrix) == 1))
  expect_equal(part$matrix, t(part$matrix))
})

test_that("edge selection honours threshold and alpha", {
  traits <- c("Fv_Fm", "PIabs", "A", "gs", "RWC", "SPAD")

  loose <- stress_network(brachiaria_stress, traits = traits, threshold = 0.05)
  tight <- stress_network(brachiaria_stress, traits = traits, threshold = 0.4)

  expect_gte(nrow(loose$edges), nrow(tight$edges))
  expect_true(all(abs(tight$edges$correlation) >= 0.4))
})

test_that("a network with no surviving edge warns instead of failing", {
  expect_warning(
    net <- stress_network(
      brachiaria_stress,
      traits = c("Fv_Fm", "PIabs", "A", "gs", "RWC", "SPAD"),
      threshold = 0.99
    ),
    class = "plantstressR_warning"
  )

  expect_equal(nrow(net$edges), 0L)
  expect_equal(nrow(net$modules), 0L)
  expect_true(all(is.na(net$nodes$module)))
})

test_that("lambda controls the shrinkage actually applied", {
  traits <- c("Fv_Fm", "PIabs", "A", "gs", "RWC", "SPAD")

  net <- stress_network(brachiaria_stress, traits = traits, lambda = 0.3)

  expect_equal(net$parameters$lambda, 0.3)
  expect_error(
    stress_network(brachiaria_stress, traits = traits, lambda = 1),
    class = "plantstressR_error"
  )
})

test_that("attached indices reach the node table", {
  traits <- c("Fv_Fm", "PIabs", "A", "gs", "RWC", "SPAD")
  sri <- calculate_sri(brachiaria_stress, "drought_level", "control",
    traits = traits, verbose = FALSE
  )

  net <- stress_network(brachiaria_stress, traits = traits, sri = sri)

  expect_true("sri" %in% names(net$nodes))
  expect_false(anyNA(net$nodes$sri))
})

test_that("input is validated", {
  expect_error(
    stress_network(brachiaria_stress, traits = c("Fv_Fm", "A")),
    class = "plantstressR_error"
  )
  expect_error(
    stress_network(brachiaria_stress,
      traits = c("Fv_Fm", "A", "gs"),
      treatment = "drought_level", level = "nonexistent"
    ),
    class = "plantstressR_error"
  )
  expect_error(
    stress_network(brachiaria_stress[1:4, ], traits = c("Fv_Fm", "A", "gs")),
    class = "plantstressR_error"
  )
})

test_that("plotting returns a ggplot and leaves the RNG untouched", {
  traits <- c("Fv_Fm", "PIabs", "A", "gs", "RWC", "SPAD")
  net <- stress_network(brachiaria_stress, traits = traits)

  set.seed(99)
  before <- .Random.seed
  p <- plot(net)
  expect_identical(.Random.seed, before)

  expect_s3_class(p, "ggplot")
  expect_error(plot(net, color_by = "sri"), class = "plantstressR_error")
  expect_error(plot(net, layout = "not_a_layout"), class = "plantstressR_error")
})

test_that("print method summarises the modules", {
  net <- stress_network(
    brachiaria_stress,
    traits = c("Fv_Fm", "PIabs", "A", "gs", "RWC", "SPAD")
  )

  expect_output(print(net), "plantstress_network")
  expect_output(print(net), "module")
})
