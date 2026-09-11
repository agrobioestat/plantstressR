#' Simulate a Control Versus Drought Phenotyping Trial
#'
#' @description
#' Generates a synthetic control-versus-stress trait table with biologically
#' plausible drought responses across chlorophyll a fluorescence (OJIP), gas
#' exchange, leaf chlorophyll index, water status, mineral nutrition, growth and
#' biomass traits. It exists so that the examples, tests and vignette of the
#' package can run on a reproducible data set with a known stress structure.
#'
#' Trait means shift linearly with the drought level, genotype effects are drawn
#' once per genotype, and about three per cent of the numeric cells are set to
#' `NA` so that the missing-value handling of the package is exercised.
#'
#' @param n Number of plants to simulate (at least 30).
#' @param seed Random seed passed to [set.seed()].
#'
#' @return A tibble with `n` rows and 33 columns; see [brachiaria_stress] for the
#'   column dictionary.
#'
#' @seealso [brachiaria_stress], [calculate_sri()]
#'
#' @examples
#' dat <- simulate_brachiaria_stress(n = 60, seed = 123)
#' table(dat$drought_level)
#' @export
simulate_brachiaria_stress <- function(n = 120, seed = 123) {
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 30) {
    rlang::abort("`n` must be a numeric scalar >= 30.")
  }
  n <- as.integer(n)

  # The caller's random stream is restored when the function returns.
  ps_local_seed(seed)

  genotype_levels <- paste0("G", 1:4)
  drought_levels <- c("control", "moderate", "severe")
  treatment_levels <- c("well_watered", "moderate_deficit", "severe_deficit")

  drought_level <- sample(drought_levels, size = n, replace = TRUE, prob = c(0.35, 0.35, 0.30))
  drought_num <- as.integer(factor(drought_level, levels = drought_levels)) - 1L
  genotype <- sample(genotype_levels, size = n, replace = TRUE)
  treatment <- treatment_levels[match(drought_level, drought_levels)]

  g_eff <- stats::rnorm(length(genotype_levels), mean = 0, sd = 0.25)
  names(g_eff) <- genotype_levels
  ge <- g_eff[genotype]
  eps <- function(sd) stats::rnorm(n, mean = 0, sd = sd)

  # Latent physiological factors shared by functionally related traits. They are
  # what makes the trait network of this data set modular: without them, traits
  # would be conditionally independent given the drought level and every partial
  # correlation would be zero by construction.
  f_photo <- eps(1)
  f_stomata <- eps(1)
  f_water <- eps(1)
  f_nutrition <- eps(1)
  f_growth <- eps(1)

  Fo <- 320 + 8 * drought_num + 4 * ge - 8.0 * f_photo + eps(7)
  Fm <- 1250 - 90 * drought_num + 8 * ge + 26.0 * f_photo + eps(20)
  Fv_Fm <- 0.82 - 0.035 * drought_num + 0.003 * ge + 0.010 * f_photo + eps(0.006)
  ABS_RC <- 1.75 + 0.08 * drought_num + 0.02 * ge - 0.060 * f_photo + eps(0.045)
  TR0_RC <- 1.35 - 0.11 * drought_num + 0.02 * ge + 0.048 * f_photo + eps(0.032)
  ET0_RC <- 1.05 - 0.20 * drought_num + 0.02 * ge + 0.062 * f_photo + eps(0.045)
  DI0_RC <- 0.25 + 0.12 * drought_num + 0.01 * ge - 0.030 * f_photo + eps(0.022)
  PIabs <- 3.8 - 1.0 * drought_num + 0.06 * ge + 0.290 * f_photo + eps(0.180)

  A <- 26 - 6.5 * drought_num + 0.6 * ge + 1.55 * f_stomata + 0.80 * f_photo + eps(1.1)
  gs <- 0.36 - 0.12 * drought_num + 0.01 * ge + 0.030 * f_stomata + eps(0.020)
  E <- 6.8 - 1.9 * drought_num + 0.07 * ge + 0.330 * f_stomata + eps(0.230)
  Ci <- 255 + 12 * drought_num + 1.2 * ge - 13.0 * f_stomata + eps(9)
  WUE <- 3.5 + 0.30 * drought_num + 0.05 * ge - 0.230 * f_stomata + eps(0.160)

  SPAD <- 42 - 4.0 * drought_num + 0.5 * ge + 1.60 * f_water + 0.85 * f_photo + eps(1.3)
  RWC <- 88 - 12.5 * drought_num + 0.7 * ge + 2.30 * f_water + eps(1.7)
  leaf_water_potential <- -0.55 - 0.42 * drought_num + 0.01 * ge +
    0.085 * f_water + eps(0.055)

  N <- 25.5 - 2.4 * drought_num + 0.2 * ge + 0.75 * f_nutrition + eps(0.55)
  P <- 3.3 - 0.45 * drought_num + 0.03 * ge + 0.135 * f_nutrition + eps(0.095)
  K <- 18.5 - 2.0 * drought_num + 0.15 * ge + 0.68 * f_nutrition + eps(0.50)
  Ca <- 7.8 - 0.6 * drought_num + 0.05 * ge + 0.330 * f_nutrition + eps(0.240)
  Mg <- 3.9 - 0.3 * drought_num + 0.04 * ge + 0.165 * f_nutrition + eps(0.120)

  plant_height <- 72 - 13 * drought_num + 1.8 * ge + 3.7 * f_growth + eps(2.6)
  tiller_number <- 9.2 - 1.8 * drought_num + 0.2 * ge + 0.74 * f_growth + eps(0.52)
  leaf_area <- 320 - 85 * drought_num + 6 * ge + 21.0 * f_growth + eps(14)
  root_length <- 27 - 3.8 * drought_num + 0.4 * ge + 1.60 * f_growth + eps(1.15)

  shoot_biomass <- 46 - 13 * drought_num + 1.2 * ge + 3.1 * f_growth +
    0.9 * f_stomata + eps(2.0)
  root_biomass <- 21 - 5.2 * drought_num + 0.8 * ge + 1.55 * f_growth + eps(1.05)
  total_biomass <- shoot_biomass + root_biomass + eps(1.2)
  yield_proxy <- 9.8 - 3.3 * drought_num + 0.2 * ge + 0.75 * f_growth + eps(0.50)

  dat <- tibble::tibble(
    sample_id = sprintf("S%03d", seq_len(n)),
    genotype = genotype,
    treatment = treatment,
    drought_level = factor(drought_level, levels = drought_levels),
    Fo = Fo,
    Fm = Fm,
    Fv_Fm = Fv_Fm,
    ABS_RC = ABS_RC,
    TR0_RC = TR0_RC,
    ET0_RC = ET0_RC,
    DI0_RC = DI0_RC,
    PIabs = PIabs,
    A = A,
    gs = gs,
    E = E,
    Ci = Ci,
    WUE = WUE,
    SPAD = SPAD,
    RWC = RWC,
    leaf_water_potential = leaf_water_potential,
    N = N,
    P = P,
    K = K,
    Ca = Ca,
    Mg = Mg,
    plant_height = plant_height,
    tiller_number = tiller_number,
    leaf_area = leaf_area,
    root_length = root_length,
    shoot_biomass = shoot_biomass,
    root_biomass = root_biomass,
    total_biomass = total_biomass,
    yield_proxy = yield_proxy
  )

  numeric_cols <- names(dat)[vapply(dat, is.numeric, logical(1))]
  numeric_cols <- setdiff(numeric_cols, "sample_id")
  for (nm in numeric_cols) {
    vals <- dat[[nm]]
    if (nm == "leaf_water_potential") {
      dat[[nm]] <- pmin(vals, -0.05)
    } else if (nm %in% c("Fv_Fm", "gs", "P", "Mg", "Ca")) {
      dat[[nm]] <- pmax(vals, 0.001)
    } else {
      dat[[nm]] <- pmax(vals, 0)
    }
  }

  set.seed(seed + 1L)
  n_missing <- max(1L, floor(n * length(numeric_cols) * 0.03))
  miss_rows <- sample(seq_len(n), n_missing, replace = TRUE)
  miss_cols <- sample(numeric_cols, n_missing, replace = TRUE)
  for (i in seq_len(n_missing)) {
    dat[[miss_cols[i]]][miss_rows[i]] <- NA_real_
  }

  dat
}

#' Synthetic Brachiaria Drought Phenotyping Trial
#'
#' @description
#' A synthetic control-versus-drought phenotyping trial used throughout the
#' documentation of `plantstressR`. Four genotypes were measured under three
#' water regimes for 29 physiological, growth and productivity traits.
#'
#' @format A tibble with 120 rows and 33 columns:
#' \describe{
#'   \item{sample_id}{Sample identifier}
#'   \item{genotype}{Genotype label}
#'   \item{treatment}{Water treatment label}
#'   \item{drought_level}{Stress class (`control`, `moderate`, `severe`)}
#'   \item{Fo, Fm, Fv_Fm, ABS_RC, TR0_RC, ET0_RC, DI0_RC, PIabs}{OJIP traits}
#'   \item{A, gs, E, Ci, WUE}{Gas exchange traits}
#'   \item{SPAD}{Leaf chlorophyll index}
#'   \item{RWC, leaf_water_potential}{Water status traits}
#'   \item{N, P, K, Ca, Mg}{Mineral nutrition traits}
#'   \item{plant_height, tiller_number, leaf_area, root_length}{Growth traits}
#'   \item{shoot_biomass, root_biomass, total_biomass, yield_proxy}{Productivity traits}
#' }
#'
#' @source Simulated with `simulate_brachiaria_stress()`.
#' @keywords datasets
#' @name brachiaria_stress
#' @docType data
NULL
