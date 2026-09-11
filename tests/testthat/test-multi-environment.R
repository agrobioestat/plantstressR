met_trial <- function() {
  dat <- toy_block_trial()
  dat$site <- ifelse(dat$block %in% c("B1", "B2"), "north", "south")
  # The south site pushes genotype B harder than the north one, so the two
  # sites genuinely disagree about B and agree about A.
  hit <- dat$site == "south" & dat$genotype == "B" & dat$trt == "stress"
  dat$down[hit] <- dat$down[hit] - 4
  dat
}

test_that("by accepts several columns and builds one label per cell", {
  dat <- met_trial()
  sri <- calculate_sri(dat,
    treatment = "trt", control = "ctrl", traits = c("up", "down"),
    by = c("genotype", "site"), verbose = FALSE
  )

  expect_setequal(
    unique(sri$unit),
    c("A | north", "A | south", "B | north", "B | south")
  )
  expect_equal(attr(sri, "plantstress")$by, c("genotype", "site"))

  # One column must behave exactly as before.
  single <- calculate_sri(dat,
    treatment = "trt", control = "ctrl", traits = c("up", "down"),
    by = "genotype", verbose = FALSE
  )
  expect_setequal(unique(single$unit), c("A", "B"))
})

test_that("validate_stress_data checks each grouping column", {
  dat <- met_trial()
  chk <- validate_stress_data(dat,
    treatment = "trt", control = "ctrl",
    by = c("genotype", "site"), verbose = FALSE
  )
  expect_equal(chk$by, c("genotype", "site"))
  expect_false(any(c("genotype", "site") %in% chk$traits))

  expect_error(
    validate_stress_data(dat,
      treatment = "trt", control = "ctrl",
      by = c("genotype", "nowhere"), verbose = FALSE
    ),
    class = "plantstressR_error"
  )
  expect_error(
    validate_stress_data(dat,
      treatment = "trt", control = "ctrl",
      by = c("genotype", "genotype"), verbose = FALSE
    ),
    class = "plantstressR_error"
  )

  # A cell without control plants cannot be indexed.
  broken <- dat[!(dat$site == "south" & dat$genotype == "B" & dat$trt == "ctrl"), ]
  expect_error(
    validate_stress_data(broken,
      treatment = "trt", control = "ctrl",
      by = c("genotype", "site"), verbose = FALSE
    ),
    class = "plantstressR_error"
  )
})

test_that("stress_stability splits the label back and summarises environments", {
  dat <- met_trial()
  isi <- integrated_stress_index(calculate_sri(dat,
    treatment = "trt", control = "ctrl", traits = c("up", "down"),
    by = c("genotype", "site"), verbose = FALSE
  ))

  stab <- stress_stability(isi)
  expect_s3_class(stab, "plantstress_stability")
  expect_setequal(stab$unit, c("A", "B"))
  expect_true(all(stab$n_env == 2L))
  expect_true(all(stab$rank_min <= stab$mean_rank))
  expect_true(all(stab$rank_max >= stab$mean_rank))

  # B is the genotype the two sites disagree about, so it must be the less
  # consistent of the two.
  expect_gt(
    stab$sd_isi[stab$unit == "B"],
    stab$sd_isi[stab$unit == "A"]
  )
})

test_that("ecovalence partitions the interaction", {
  dat <- met_trial()
  stab <- stress_stability(integrated_stress_index(calculate_sri(dat,
    treatment = "trt", control = "ctrl", traits = c("up", "down"),
    by = c("genotype", "site"), verbose = FALSE
  )))

  expect_true(all(stab$ecovalence >= 0))
  totals <- tapply(stab$ecovalence_pct, stab$group, sum)
  expect_true(all(abs(totals - 100) < 1e-8))
})

test_that("a trial with no genotype x environment interaction has zero ecovalence", {
  # Both genotypes shifted by the same amount at the second site.
  dat <- toy_block_trial()
  dat$site <- ifelse(dat$block %in% c("B1", "B2"), "north", "south")
  shift <- dat$site == "south"
  dat$up[shift] <- dat$up[shift] + 2
  dat$down[shift] <- dat$down[shift] + 2

  stab <- stress_stability(integrated_stress_index(
    calculate_sri(dat,
      treatment = "trt", control = "ctrl", traits = c("up", "down"),
      by = c("genotype", "site"), verbose = FALSE
    ),
    rescale = FALSE
  ))
  expect_true(all(stab$ecovalence < 1e-8))
})

test_that("stress_stability refuses what it cannot summarise", {
  dat <- met_trial()

  single <- integrated_stress_index(calculate_sri(dat,
    treatment = "trt", control = "ctrl", traits = c("up", "down"),
    by = "genotype", verbose = FALSE
  ))
  expect_error(stress_stability(single), class = "plantstressR_error")

  multi <- integrated_stress_index(calculate_sri(dat,
    treatment = "trt", control = "ctrl", traits = c("up", "down"),
    by = c("genotype", "site"), verbose = FALSE
  ))
  expect_error(stress_stability(multi, environment = "nowhere"),
    class = "plantstressR_error"
  )
  expect_error(stress_stability(as.data.frame(multi)),
    class = "plantstressR_error"
  )

  # Naming the genotype as the environment simply swaps the two roles.
  swapped <- stress_stability(multi, environment = "genotype")
  expect_setequal(swapped$unit, c("north", "south"))
})

test_that("the bootstrap follows a multi-environment design", {
  dat <- met_trial()
  sri <- calculate_sri(dat,
    treatment = "trt", control = "ctrl", traits = c("up", "down"),
    by = c("genotype", "site"), block = "block", verbose = FALSE
  )
  ci <- stress_index_ci(sri, n_boot = 15, seed = 4)

  expect_equal(nrow(ci), length(unique(sri$unit)) * length(unique(sri$group)))
  expect_true(all(ci$conf_low <= ci$isi))
})

test_that("stability has a plot method", {
  skip_if_not_installed("ggplot2")
  dat <- met_trial()
  stab <- stress_stability(integrated_stress_index(calculate_sri(dat,
    treatment = "trt", control = "ctrl", traits = c("up", "down"),
    by = c("genotype", "site"), verbose = FALSE
  )))
  expect_s3_class(plot(stab), "ggplot")
})
