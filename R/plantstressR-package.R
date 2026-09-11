#' @keywords internal
#' @aliases plantstressR-package
#'
#' @description
#' `plantstressR` answers a single question: *what is the physiological stress
#' signature of a plant, and how do its traits interact?* The workflow is
#' deliberately short and always starts from one tidy table of physiological
#' traits measured in a control treatment and in one or more stress treatments.
#'
#' 1. [validate_stress_data()] checks the design (control level, replication,
#'    numeric traits, missing values) before anything is computed.
#' 2. [calculate_sri()] converts every trait into a signed, standardized
#'    Stress Response Index (SRI) relative to the control treatment.
#' 3. [integrated_stress_index()] collapses the per-trait indices into one
#'    weighted integrated stress index (ISI) used to rank genotypes or
#'    treatments.
#' 4. [plot_stress_signature()] draws the signature as a heat map or a radar
#'    plot.
#' 5. [stress_network()] estimates a regularized partial-correlation network
#'    among traits and detects the modules that respond to stress together, and
#'    [stress_module_scores()] scores those modules with the trait indices.
#' 6. [stress_ordination()] places samples in a multivariate trait space when a
#'    principal component view is required.
#'
#' [stress_tolerance_index()] sits beside this workflow: it computes the
#' classical productivity-based selection indices (`STI`, `SSI`, `GMP`, `TOL`
#' and others) from a single yield trait, so that a physiological ranking can be
#' checked against the agronomic one.
#'
#' @importFrom rlang .data
"_PACKAGE"
