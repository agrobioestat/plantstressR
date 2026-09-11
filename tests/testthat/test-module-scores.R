module_fixture <- function() {
  traits <- c("Fv_Fm", "PIabs", "DI0_RC", "A", "gs", "Ci", "RWC", "SPAD")
  list(
    traits = traits,
    sri = calculate_sri(brachiaria_stress,
      treatment = "drought_level", control = "control",
      traits = traits, verbose = FALSE
    ),
    net = stress_network(brachiaria_stress,
      traits = traits, treatment = "drought_level",
      level = c("moderate", "severe")
    )
  )
}

test_that("a module score is the mean index of its traits", {
  f <- module_fixture()
  scores <- stress_module_scores(f$net, f$sri)

  module_of <- stats::setNames(f$net$nodes$module, f$net$nodes$trait)
  row <- scores[scores$group == "severe", ][1, ]
  members <- names(module_of)[!is.na(module_of) & module_of == row$module]
  expected <- mean(f$sri$sri[f$sri$group == "severe" & f$sri$trait %in% members])

  expect_equal(row$score, expected)
  expect_equal(row$n_traits, length(members))
  expect_equal(sort(strsplit(row$traits, ", ")[[1]]), sort(members))
})

test_that("the rms aggregation is a magnitude", {
  f <- module_fixture()
  signed <- stress_module_scores(f$net, f$sri, aggregate = "mean")
  magnitude <- stress_module_scores(f$net, f$sri, aggregate = "rms")

  expect_true(all(magnitude$score >= 0))
  expect_true(all(magnitude$score >= abs(signed$score) - 1e-8))
})

test_that("scores are ordered by severity within each unit and level", {
  f <- module_fixture()
  scores <- stress_module_scores(f$net, f$sri)

  for (g in unique(scores$group)) {
    expect_false(is.unsorted(rev(scores$score[scores$group == g])))
  }
  expect_true(all(scores$hub %in% f$traits))
})

test_that("traits outside the network are dropped with a warning", {
  f <- module_fixture()
  sri_extra <- calculate_sri(brachiaria_stress,
    treatment = "drought_level", control = "control",
    traits = c(f$traits, "root_biomass"), verbose = FALSE
  )

  expect_warning(
    scores <- stress_module_scores(f$net, sri_extra),
    class = "plantstressR_warning"
  )
  expect_false(any(grepl("root_biomass", scores$traits)))
})

test_that("unassigned traits are excluded unless requested", {
  f <- module_fixture()
  # A threshold this tight leaves every trait isolated and unassigned.
  sparse <- suppressWarnings(
    stress_network(brachiaria_stress, traits = f$traits, threshold = 0.99)
  )

  expect_error(stress_module_scores(sparse, f$sri), class = "plantstressR_error")

  pooled <- stress_module_scores(sparse, f$sri, include_unassigned = TRUE)
  expect_true(all(pooled$label == "unassigned"))
  expect_equal(nrow(pooled), length(unique(f$sri$group)))
})

test_that("input is validated", {
  f <- module_fixture()

  expect_error(stress_module_scores(f$sri, f$sri), class = "plantstressR_error")
  expect_error(stress_module_scores(f$net, data.frame(a = 1)), class = "plantstressR_error")

  other_panel <- calculate_sri(brachiaria_stress,
    treatment = "drought_level", control = "control",
    traits = c("N", "P", "K"), verbose = FALSE
  )
  expect_error(
    suppressWarnings(stress_module_scores(f$net, other_panel)),
    class = "plantstressR_error"
  )
})
