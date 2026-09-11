test_that("the bundled trials have the shape their documentation claims", {
  data(brachiaria_stress, envir = environment())
  data(wheat_salinity, envir = environment())
  data(maize_heat, envir = environment())

  expect_equal(dim(brachiaria_stress), c(120L, 34L))
  expect_equal(dim(wheat_salinity), c(162L, 20L))
  expect_equal(dim(maize_heat), c(38L, 12L))

  # A stray `names` attribute on a trait column is easy to introduce when the
  # simulation indexes a named effect vector, and it survives into every result.
  for (dat in list(brachiaria_stress, wheat_salinity, maize_heat)) {
    expect_equal(sum(vapply(dat, function(x) !is.null(names(x)), logical(1))), 0L)
  }
})

test_that("wheat_salinity is a complete multi-environment layout", {
  data(wheat_salinity, envir = environment())

  cells <- table(wheat_salinity$site, wheat_salinity$salinity)
  expect_true(all(cells == 18L))
  expect_setequal(levels(wheat_salinity$salinity), c("control", "moderate", "severe"))
  expect_equal(length(unique(wheat_salinity$cultivar)), 6L)

  # Every cultivar is present at every site, which is what ecovalence needs.
  expect_true(all(table(wheat_salinity$cultivar, wheat_salinity$site) == 9L))

  # The nominal dose tracks the treatment class.
  expect_equal(
    as.vector(tapply(wheat_salinity$nacl_mM, wheat_salinity$salinity, unique)),
    c(0, 75, 150)
  )
})

test_that("the salinity dictionary resolves the ion traits unaided", {
  dirs <- trait_directions(c(
    "Fv_Fm", "A", "RWC", "Na", "Cl", "K", "MDA",
    "electrolyte_leakage", "grain_yield"
  ))
  expect_equal(unname(dirs[c("Na", "Cl", "MDA", "electrolyte_leakage")]),
    rep("lower_is_better", 4)
  )
  # Potassium is excluded from the shoot under salinity: falling K is damage,
  # so it keeps the default orientation.
  expect_equal(unname(dirs[c("K", "Fv_Fm", "grain_yield")]),
    rep("higher_is_better", 3)
  )
})

test_that("wheat_salinity carries a real genotype x environment interaction", {
  data(wheat_salinity, envir = environment())

  sri <- calculate_sri(wheat_salinity,
    treatment = "salinity", control = "control",
    traits = c("Fv_Fm", "A", "RWC", "Na", "K", "MDA", "grain_yield"),
    by = c("cultivar", "site"), block = "block", verbose = FALSE
  )
  stab <- stress_stability(integrated_stress_index(sri))

  expect_true(all(stab$n_env == 3L))
  expect_true(all(is.finite(stab$ecovalence)))
  # W4 is built to rank near the top at one site and near the bottom at
  # another; the spread of its ranks has to show that.
  w4 <- stab[stab$unit == "W4" & stab$group == "severe", ]
  expect_gt(w4$rank_max - w4$rank_min, 2)
})

test_that("maize_heat triggers every design check it was built to trigger", {
  data(maize_heat, envir = environment())

  chk <- validate_stress_data(maize_heat,
    treatment = "treatment", control = "control",
    by = "hybrid", verbose = FALSE
  )
  expect_setequal(
    chk$issues$check,
    c("replication", "missing_values", "constant_trait")
  )

  # The thin cell, the abandoned trait and the constant column, respectively.
  expect_equal(sum(maize_heat$hybrid == "H4" & maize_heat$treatment == "control"), 2L)
  expect_gt(mean(is.na(maize_heat$gs)), 0.3)
  expect_equal(stats::sd(maize_heat$plot_area), 0)
})

test_that("the messy trial still analyses once the constant column is dropped", {
  data(maize_heat, envir = environment())

  sri <- calculate_sri(maize_heat,
    treatment = "treatment", control = "control",
    traits = c("Fv_Fm", "F0", "A", "RWC", "kernel_number"),
    by = "hybrid", verbose = FALSE
  )
  expect_true(all(is.finite(sri$sri)))

  # F0 rises under heat and that is damage, so its index must be positive.
  expect_true(all(sri$sri[sri$trait == "F0"] > 0))
  expect_true(all(sri$sri[sri$trait == "Fv_Fm"] > 0))

  # Unequal replication is exactly when precision weighting earns its keep.
  w <- stress_weights(integrated_stress_index(sri, weights = "precision"))
  expect_equal(sum(w$weight), 1, tolerance = 1e-8)
})
