# Regenerates data/wheat_salinity.rda.
# Run from the package root with: source("data-raw/create_wheat_salinity.R")
#
# A synthetic multi-environment salinity trial. Its job in the package is to be
# the data set the multi-environment functions were written for: six cultivars
# grown at three sites, in complete blocks within each site, so that
# stress_stability() has a real genotype-by-environment table instead of one
# faked by splitting blocks. The ion-relation traits (Na, Cl, K) also exercise
# the direction dictionary, which knows that rising sodium is damage.

set.seed(2024)

cultivars <- paste0("W", 1:6)
sites <- c("coastal", "inland", "upland")
blocks <- paste0("B", 1:3)
levels_salt <- c("control", "moderate", "severe")
nacl <- c(control = 0, moderate = 75, severe = 150)

design <- expand.grid(
  block = blocks,
  cultivar = cultivars,
  salinity = levels_salt,
  site = sites,
  stringsAsFactors = FALSE
)
n <- nrow(design)

# Dose on a 0-1 scale drives every response.
dose <- unname(nacl[design$salinity]) / 150

# Cultivar tolerance: 0 = susceptible, 1 = tolerant. W1 and W2 are the good
# material, W5 and W6 the poor material.
tol_base <- c(W1 = 0.85, W2 = 0.75, W3 = 0.55, W4 = 0.45, W5 = 0.25, W6 = 0.15)

# Genotype by environment: W2 holds up inland but collapses on the coast, and
# W4 is the opposite. Without this the stability table would be empty of news.
gxe <- matrix(0, nrow = length(cultivars), ncol = length(sites),
              dimnames = list(cultivars, sites))
gxe["W2", "coastal"] <- -0.30
gxe["W2", "upland"] <- 0.08
gxe["W4", "coastal"] <- 0.28
gxe["W4", "inland"] <- -0.12

tol <- unname(tol_base[design$cultivar]) +
  gxe[cbind(design$cultivar, design$site)]

# Sites differ in how hard the same nominal dose bites.
site_severity <- c(coastal = 1.25, inland = 1.00, upland = 0.80)
severity <- dose * unname(site_severity[design$site]) * (1 - 0.6 * tol)

# Blocks are nuisance strata nested in the site.
block_key <- paste(design$site, design$block, sep = "_")
block_eff <- stats::setNames(
  stats::rnorm(length(unique(block_key)), 0, 0.5),
  unique(block_key)
)
be <- unname(block_eff[block_key])

eps <- function(sd) stats::rnorm(n, 0, sd)

# Latent factors shared by functionally related traits, so that the trait
# network of this data set has genuine modules.
f_photo <- eps(1)
f_stomata <- eps(1)
f_ion <- eps(1)
f_damage <- eps(1)
f_growth <- eps(1)

dat <- data.frame(
  plot_id = sprintf("P%03d", seq_len(n)),
  cultivar = design$cultivar,
  site = design$site,
  block = design$block,
  salinity = factor(design$salinity, levels = levels_salt),
  nacl_mM = unname(nacl[design$salinity]),

  Fv_Fm = 0.81 - 0.13 * severity + 0.011 * f_photo + eps(0.007),
  PIabs = 4.10 - 2.60 * severity + 0.300 * f_photo + eps(0.170),

  A = 24.0 - 12.5 * severity + 1.60 * f_stomata + 0.75 * f_photo + eps(1.05),
  gs = 0.340 - 0.190 * severity + 0.031 * f_stomata + eps(0.019),
  E = 6.40 - 3.30 * severity + 0.320 * f_stomata + eps(0.210),

  RWC = 89.0 - 21.0 * severity + 2.10 * f_growth + eps(1.60),
  SPAD = 44.0 - 11.0 * severity + 1.55 * f_photo + eps(1.25),

  # Ion relations: sodium and chloride accumulate, potassium is excluded.
  Na = 0.85 + 6.30 * severity + 0.210 * f_ion + eps(0.140),
  Cl = 1.10 + 7.10 * severity + 0.260 * f_ion + eps(0.180),
  K = 3.90 - 1.55 * severity - 0.115 * f_ion + eps(0.095),

  # Oxidative damage.
  MDA = 4.20 + 9.40 * severity + 0.330 * f_damage + eps(0.240),
  electrolyte_leakage = 12.0 + 38.0 * severity + 1.45 * f_damage + eps(1.10),

  shoot_biomass = 38.0 - 19.5 * severity + 2.90 * f_growth +
    0.85 * f_stomata + eps(1.90),
  grain_yield = 6.30 - 3.60 * severity + 0.500 * f_growth + eps(0.360),
  stringsAsFactors = FALSE
)

trait_cols <- setdiff(
  names(dat)[vapply(dat, is.numeric, logical(1))],
  "nacl_mM"
)

# The block shift acts on each trait's own scale.
for (nm in trait_cols) {
  dat[[nm]] <- dat[[nm]] + be * stats::sd(dat[[nm]])
}

# Keep every trait inside its physically possible range.
for (nm in trait_cols) {
  dat[[nm]] <- if (nm == "Fv_Fm") {
    pmin(pmax(dat[[nm]], 0.05), 0.85)
  } else if (nm %in% c("RWC", "electrolyte_leakage")) {
    pmin(pmax(dat[[nm]], 1), 100)
  } else {
    pmax(dat[[nm]], 0.001)
  }
}

# A little missingness, as any real trial has.
set.seed(2025)
n_missing <- floor(n * length(trait_cols) * 0.02)
miss_rows <- sample(seq_len(n), n_missing, replace = TRUE)
miss_cols <- sample(trait_cols, n_missing, replace = TRUE)
for (i in seq_len(n_missing)) {
  dat[[miss_cols[i]]][miss_rows[i]] <- NA_real_
}

wheat_salinity <- tibble::as_tibble(dat)

save(
  wheat_salinity,
  file = file.path("data", "wheat_salinity.rda"),
  compress = "bzip2",
  version = 3
)
