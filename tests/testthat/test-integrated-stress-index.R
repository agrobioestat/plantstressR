sri_fixture <- function() {
  calculate_sri(toy_trial(), "trt", "ctrl",
    traits = c("up", "down"), by = "genotype", verbose = FALSE
  )
}

test_that("equal weights reproduce the arithmetic mean of the indices", {
  sri <- sri_fixture()
  isi <- integrated_stress_index(sri, weights = "equal")

  expected <- vapply(split(sri$sri, sri$unit), mean, numeric(1))
  got <- stats::setNames(isi$isi, isi$unit)

  expect_equal(got[names(expected)], expected)
})

test_that("weights always sum to one", {
  sri <- sri_fixture()

  for (w in c("equal", "precision", "pca")) {
    isi <- suppressWarnings(integrated_stress_index(sri, weights = w))
    expect_equal(sum(stress_weights(isi)$weight), 1)
  }
})

test_that("manual weights are honoured and rescaled", {
  sri <- sri_fixture()
  isi <- integrated_stress_index(sri, weights = c(up = 3, down = 1))
  w <- stress_weights(isi)

  expect_equal(w$weight[w$trait == "up"], 0.75)
  expect_equal(w$weight[w$trait == "down"], 0.25)
})

test_that("ranking happens within a stress level", {
  moderate <- as.data.frame(sri_fixture())
  severe <- moderate
  severe$group <- "severe"
  severe$sri <- severe$sri * 3
  sri <- rbind(moderate, severe)

  isi <- integrated_stress_index(sri)

  ranks <- split(isi$rank, isi$group)
  expect_true(all(vapply(ranks, function(r) setequal(r, 1:2), logical(1))))
  expect_true(all(range(isi$isi_scaled, na.rm = TRUE) == c(0, 100)))
})

test_that("rank_by chooses the direction of the ranking", {
  sri <- sri_fixture()

  tol <- integrated_stress_index(sri, rank_by = "tolerance")
  sev <- integrated_stress_index(sri, rank_by = "severity")

  expect_equal(tol$unit[tol$rank == 1], sev$unit[sev$rank == max(sev$rank)])
})

test_that("rms aggregation never returns a negative index", {
  sri <- sri_fixture()
  sri$sri <- sri$sri * c(1, -1)

  isi <- integrated_stress_index(sri, aggregate = "rms")

  expect_true(all(isi$isi >= 0))
})

test_that("contributions decompose the index", {
  sri <- sri_fixture()
  isi <- integrated_stress_index(sri, weights = "equal")
  contrib <- stress_contributions(isi)

  totals <- vapply(split(contrib$contribution, contrib$unit), sum, numeric(1))
  expect_equal(unname(totals[isi$unit]), isi$isi)
})

test_that("invalid input is rejected", {
  expect_error(integrated_stress_index(data.frame(a = 1)), class = "plantstressR_error")
  expect_error(
    integrated_stress_index(sri_fixture(), weights = c(up = 1)),
    class = "plantstressR_error"
  )
  expect_error(stress_weights(data.frame()), class = "plantstressR_error")
})

test_that("print method works", {
  expect_output(print(integrated_stress_index(sri_fixture())), "plantstress_isi")
})
