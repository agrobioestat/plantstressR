test_that("validate_stress_data checks the block layout", {
  dat <- toy_block_trial()

  chk <- validate_stress_data(dat,
    treatment = "trt", control = "ctrl",
    by = "genotype", block = "block", verbose = FALSE
  )
  expect_equal(chk$block, "block")
  expect_false("block" %in% chk$traits)

  one_block <- dat[dat$block == "B1", , drop = FALSE]
  expect_error(
    validate_stress_data(one_block,
      treatment = "trt", control = "ctrl",
      block = "block", verbose = FALSE
    ),
    class = "plantstressR_error"
  )

  # A treatment confined to a single block cannot be told apart from it.
  confounded <- dat
  confounded$block[confounded$trt == "stress"] <- "B1"
  expect_error(
    validate_stress_data(confounded,
      treatment = "trt", control = "ctrl",
      block = "block", verbose = FALSE
    ),
    class = "plantstressR_error"
  )
})

test_that("block adjustment removes the block spread but keeps the difference", {
  dat <- toy_block_trial(block_shift = 5)

  blocked <- calculate_sri(dat,
    treatment = "trt", control = "ctrl", traits = c("up", "down"),
    by = "genotype", block = "block", verbose = FALSE
  )
  unblocked <- calculate_sri(dat,
    treatment = "trt", control = "ctrl", traits = c("up", "down"),
    by = "genotype", verbose = FALSE
  )

  # The raw treatment difference is a property of the data, not of the model.
  expect_equal(
    blocked$difference[blocked$unit == "A" & blocked$trait == "up"],
    unblocked$difference[unblocked$unit == "A" & unblocked$trait == "up"],
    tolerance = 1e-8
  )

  # The scale, however, must lose the between-block variance.
  expect_lt(
    blocked$sd_control[blocked$unit == "A" & blocked$trait == "up"],
    unblocked$sd_control[unblocked$unit == "A" & unblocked$trait == "up"]
  )
  expect_gt(
    abs(blocked$sri[blocked$unit == "A" & blocked$trait == "up"]),
    abs(unblocked$sri[unblocked$unit == "A" & unblocked$trait == "up"])
  )
})

test_that("the sign convention survives the blocked path", {
  dat <- toy_block_trial()
  sri <- calculate_sri(dat,
    treatment = "trt", control = "ctrl", traits = c("up", "down"),
    by = "genotype", block = "block", verbose = FALSE
  )

  # `up` rises under stress and defaults to higher_is_better, so its index is
  # negative (an improvement); `down` falls, so its index is positive (damage).
  expect_lt(sri$sri[sri$unit == "A" & sri$trait == "up"], 0)
  expect_gt(sri$sri[sri$unit == "A" & sri$trait == "down"], 0)

  # Genotype A was pushed twice as hard as B.
  expect_gt(
    abs(sri$sri[sri$unit == "A" & sri$trait == "down"]),
    abs(sri$sri[sri$unit == "B" & sri$trait == "down"])
  )
})

test_that("a blocked fit reports the model's degrees of freedom", {
  dat <- toy_block_trial()
  sri <- calculate_sri(dat,
    treatment = "trt", control = "ctrl", traits = c("up", "down"),
    block = "block", verbose = FALSE
  )
  expect_true(all(is.finite(sri$p_value)))
  expect_true(all(sri$se_sampling > 0, na.rm = TRUE))
  expect_true(all(sri$se >= sri$se_sampling, na.rm = TRUE))
})

test_that("stress_network and stress_ordination accept a block", {
  data(brachiaria_stress, envir = environment())
  traits <- c("Fv_Fm", "PIabs", "A", "gs", "RWC", "SPAD", "shoot_biomass")

  net <- stress_network(brachiaria_stress,
    traits = traits, treatment = "drought_level", block = "block"
  )
  expect_s3_class(net, "plantstress_network")
  expect_equal(net$parameters$block, "block")
  expect_false("block" %in% net$nodes$trait)

  ord <- stress_ordination(brachiaria_stress,
    traits = traits, treatment = "drought_level", block = "block"
  )
  expect_s3_class(ord, "plantstress_ordination")
  expect_false("block" %in% ord$traits)
})

test_that("block_adjust leaves treatment means untouched", {
  dat <- toy_block_trial(block_shift = 7)
  adjusted <- block_adjust(dat$up, dat$trt, dat$block)

  expect_equal(
    tapply(adjusted, dat$trt, mean),
    tapply(dat$up, dat$trt, mean),
    tolerance = 1e-8
  )
  expect_lt(stats::sd(adjusted), stats::sd(dat$up))

  # Nothing to remove when a single block is present.
  single <- dat[dat$block == "B1", , drop = FALSE]
  expect_equal(block_adjust(single$up, single$trt, single$block), single$up)
})
