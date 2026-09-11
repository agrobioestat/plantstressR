# Regenerates data/maize_heat.rda.
# Run from the package root with: source("data-raw/create_maize_heat.R")
#
# A deliberately awkward trial. Every other data set in the package is clean,
# which leaves validate_stress_data() with nothing to report and users with no
# idea what its warnings look like. This one carries, on purpose: a hybrid with
# almost no control plants, replication that varies from cell to cell, a trait
# that was recorded but never varied, and one trait with heavy missingness.

set.seed(404)

hybrids <- paste0("H", 1:4)
treatments <- c("control", "heat")

# Unbalanced on purpose: H4 lost most of its control plot.
replication <- list(
  H1 = c(control = 5, heat = 5),
  H2 = c(control = 4, heat = 6),
  H3 = c(control = 5, heat = 5),
  H4 = c(control = 2, heat = 6)
)

rows <- do.call(rbind, lapply(hybrids, function(h) {
  do.call(rbind, lapply(treatments, function(t) {
    data.frame(
      hybrid = h, treatment = t,
      rep = seq_len(replication[[h]][[t]]),
      stringsAsFactors = FALSE
    )
  }))
}))
n <- nrow(rows)

stressed <- rows$treatment == "heat"
tolerance <- c(H1 = 0.80, H2 = 0.50, H3 = 0.35, H4 = 0.60)
hit <- stressed * (1 - 0.7 * unname(tolerance[rows$hybrid]))

eps <- function(sd) stats::rnorm(n, 0, sd)

dat <- data.frame(
  plant_id = sprintf("M%02d", seq_len(n)),
  hybrid = rows$hybrid,
  treatment = rows$treatment,

  Fv_Fm = 0.80 - 0.16 * hit + eps(0.010),
  # F0 rises under heat: the dictionary knows this one is damage when it grows.
  F0 = 318 + 74 * hit + eps(9),
  A = 27.0 - 15.0 * hit + eps(1.4),
  gs = 0.330 - 0.165 * hit + eps(0.024),
  RWC = 90.0 - 19.0 * hit + eps(2.1),
  SPAD = 45.0 - 9.5 * hit + eps(1.6),
  kernel_number = 480 - 210 * hit + eps(26),
  shoot_biomass = 42.0 - 17.0 * hit + eps(2.4),

  # Recorded on every plant, identical on every plant. Real data sets are full
  # of columns like this, and they must not silently become traits.
  plot_area = 1.5,
  stringsAsFactors = FALSE
)

dat$Fv_Fm <- pmin(pmax(dat$Fv_Fm, 0.05), 0.85)
dat$RWC <- pmin(pmax(dat$RWC, 1), 100)
for (nm in c("F0", "A", "gs", "SPAD", "kernel_number", "shoot_biomass")) {
  dat[[nm]] <- pmax(dat[[nm]], 0.001)
}

# One trait the operator gave up on halfway through the campaign.
set.seed(405)
dat$gs[sample(seq_len(n), ceiling(n * 0.35))] <- NA_real_
# And the ordinary scatter of lost readings elsewhere.
for (nm in c("SPAD", "kernel_number")) {
  dat[[nm]][sample(seq_len(n), 2)] <- NA_real_
}

maize_heat <- tibble::as_tibble(dat)

save(
  maize_heat,
  file = file.path("data", "maize_heat.rda"),
  compress = "bzip2",
  version = 3
)
