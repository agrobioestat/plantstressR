#' Synthetic Wheat Salinity Trial Across Three Sites
#'
#' @description
#' A synthetic multi-environment salinity experiment: six cultivars grown at
#' three sites, under three levels of sodium chloride, in three complete blocks
#' within each site. One plot per combination, 162 plots in all.
#'
#' This is the data set the multi-environment functions were written for. It has
#' a genuine genotype-by-environment table, so [stress_stability()] has
#' something to report. Cultivar `W2` is reliably among the best wherever it is
#' grown. `W4` is the opposite: under severe salinity it ranks second on the
#' upland site and seventeenth on the coastal one, and carries about a fifth of
#' that level's interaction on its own. Which cultivar is the least dependable
#' changes with the dose, which is itself the lesson -- at the moderate level it
#' is `W1` that reshuffles hardest.
#'
#' The ion-relation traits also exercise [trait_directions()]: `Na` and `Cl`
#' accumulate under salinity and are damage when they rise, while `K` is
#' excluded from the shoot and is damage when it falls. The dictionary resolves
#' all three without being told.
#'
#' @format A tibble with 162 rows and 20 columns:
#' \describe{
#'   \item{plot_id}{Plot identifier}
#'   \item{cultivar}{Cultivar label (`W1`-`W6`)}
#'   \item{site}{Site (`coastal`, `inland`, `upland`), differing in how hard the
#'     same nominal dose bites}
#'   \item{block}{Complete block within the site (`B1`-`B3`)}
#'   \item{salinity}{Treatment class (`control`, `moderate`, `severe`)}
#'   \item{nacl_mM}{Nominal dose in mmol per liter (0, 75, 150). A numeric
#'     **design** column: name it in `traits` only if you mean to analyze the
#'     dose as a trait}
#'   \item{Fv_Fm, PIabs}{Chlorophyll fluorescence}
#'   \item{A, gs, E}{Gas exchange}
#'   \item{RWC}{Relative water content}
#'   \item{SPAD}{Leaf chlorophyll index}
#'   \item{Na, Cl, K}{Shoot ion concentrations}
#'   \item{MDA, electrolyte_leakage}{Membrane damage}
#'   \item{shoot_biomass, grain_yield}{Productivity}
#' }
#'
#' @source Simulated by `data-raw/create_wheat_salinity.R` in the package
#'   repository. No plants were harmed, and none were grown.
#' @keywords datasets
#' @name wheat_salinity
#' @docType data
#' @seealso [stress_stability()], [brachiaria_stress], [maize_heat]
#' @examples
#' data(wheat_salinity)
#'
#' traits <- c("Fv_Fm", "A", "RWC", "Na", "K", "MDA", "grain_yield")
#'
#' sri <- calculate_sri(wheat_salinity,
#'   treatment = "salinity", control = "control", traits = traits,
#'   by = c("cultivar", "site"), block = "block", verbose = FALSE
#' )
#'
#' stress_stability(integrated_stress_index(sri))
NULL

#' Synthetic Maize Heat Trial With the Usual Problems
#'
#' @description
#' A small heat-stress experiment that is deliberately awkward, so that
#' [validate_stress_data()] has something to say. Every other data set in the
#' package is clean, which leaves users with no idea what the warnings look like
#' or what to do about them.
#'
#' Four hybrids, one heat treatment against a control, 38 plants. Built into it:
#'
#' \itemize{
#'   \item `H4` kept only two control plants, below the replication every other
#'     cell has.
#'   \item `gs` was abandoned partway through the campaign and is 37 per cent
#'     missing.
#'   \item `plot_area` was recorded on every plant and never varied. It is
#'     numeric, so it will be picked up as a trait unless `traits` is given
#'     explicitly, and its index is `NA` by construction.
#'   \item Replication differs from cell to cell, which is what makes
#'     `weights = "precision"` do something other than `"equal"`.
#' }
#'
#' `F0` is included because it is one of the traits whose *increase* means
#' damage; the dictionary in [trait_directions()] knows it.
#'
#' @format A tibble with 38 rows and 12 columns:
#' \describe{
#'   \item{plant_id}{Plant identifier}
#'   \item{hybrid}{Hybrid label (`H1`-`H4`)}
#'   \item{treatment}{`control` or `heat`}
#'   \item{Fv_Fm, F0}{Chlorophyll fluorescence}
#'   \item{A, gs}{Gas exchange, `gs` heavily incomplete}
#'   \item{RWC, SPAD}{Water status and chlorophyll index}
#'   \item{kernel_number, shoot_biomass}{Productivity}
#'   \item{plot_area}{Constant, and left in on purpose}
#' }
#'
#' @source Simulated by `data-raw/create_maize_heat.R` in the package
#'   repository.
#' @keywords datasets
#' @name maize_heat
#' @docType data
#' @seealso [validate_stress_data()], [brachiaria_stress], [wheat_salinity]
#' @examples
#' data(maize_heat)
#'
#' # Three complaints, and each one is worth acting on.
#' validate_stress_data(maize_heat,
#'   treatment = "treatment", control = "control",
#'   by = "hybrid", verbose = FALSE
#' )
#'
#' # Leaving out the constant column is the fix for the third.
#' calculate_sri(maize_heat,
#'   treatment = "treatment", control = "control",
#'   traits = c("Fv_Fm", "F0", "A", "RWC", "kernel_number"),
#'   by = "hybrid", verbose = FALSE
#' )
NULL
