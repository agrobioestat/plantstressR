#' Score the Network Modules of a Stress Signature
#'
#' @description
#' Joins the two halves of the analysis. [stress_network()] finds the modules of
#' traits that move together; [calculate_sri()] measures how far each trait moved.
#' `stress_module_scores()` averages the second within the first, turning a
#' trait-level signature into a handful of module-level statements such as
#' "the photochemical module was impaired by 2.4 standard deviations under
#' severe deficit, the stomatal module by 1.1".
#'
#' @details
#' Only traits present in both objects are used; traits that were dropped from
#' the network, or that carry no module because they were isolated, are reported
#' in the `module` column as `NA` and excluded from the scores unless
#' `include_unassigned = TRUE`, in which case they are pooled into a module
#' labeled `"unassigned"`.
#'
#' `aggregate = "mean"` keeps the sign, so a module can come out negative when
#' stress improved its traits. `aggregate = "rms"` returns the magnitude of the
#' disturbance and never cancels.
#'
#' @param network A `plantstress_network` object from [stress_network()].
#' @param sri A `plantstress_sri` object from [calculate_sri()], or a data frame
#'   with the columns `unit`, `group`, `trait` and `sri`.
#' @param aggregate `"mean"` (signed) or `"rms"` (magnitude).
#' @param include_unassigned Logical. Keep traits that were not placed in any
#'   module, pooled under the label `"unassigned"`.
#'
#' @return A tibble with columns `unit`, `group`, `module`, `label`, `n_traits`,
#'   `hub`, `score` and `traits`, sorted by decreasing score within each unit and
#'   stress level.
#'
#' @seealso [stress_network()], [calculate_sri()], [integrated_stress_index()]
#'
#' @examples
#' data(brachiaria_stress)
#'
#' traits <- c("Fv_Fm", "PIabs", "DI0_RC", "A", "gs", "Ci", "RWC", "SPAD")
#'
#' sri <- calculate_sri(brachiaria_stress,
#'   treatment = "drought_level", control = "control",
#'   traits = traits, verbose = FALSE
#' )
#'
#' net <- stress_network(brachiaria_stress,
#'   traits = traits, treatment = "drought_level",
#'   level = c("moderate", "severe")
#' )
#'
#' stress_module_scores(net, sri)
#' @export
stress_module_scores <- function(network,
                                 sri,
                                 aggregate = c("mean", "rms"),
                                 include_unassigned = FALSE) {
  aggregate <- match.arg(aggregate)

  if (!inherits(network, "plantstress_network")) {
    ps_abort("`network` must be a `plantstress_network` object from `stress_network()`.")
  }
  if (!is.data.frame(sri)) {
    ps_abort("`sri` must be a `plantstress_sri` object or a data frame.")
  }
  required <- c("unit", "group", "trait", "sri")
  missing_cols <- setdiff(required, names(sri))
  if (length(missing_cols) > 0L) {
    ps_abort(paste0(
      "`sri` is missing required column(s): ", paste(missing_cols, collapse = ", "),
      ". Build it with `calculate_sri()`."
    ))
  }
  if (!is.logical(include_unassigned) || length(include_unassigned) != 1L ||
    is.na(include_unassigned)) {
    ps_abort("`include_unassigned` must be `TRUE` or `FALSE`.")
  }

  nodes <- network$nodes[c("trait", "module", "strength")]
  x <- tibble::as_tibble(sri)[required]

  shared <- intersect(unique(x$trait), nodes$trait)
  if (length(shared) == 0L) {
    ps_abort("`network` and `sri` share no trait; they were built on different panels.")
  }
  dropped <- setdiff(unique(x$trait), nodes$trait)
  if (length(dropped) > 0L) {
    ps_warn(paste0(
      "Trait(s) absent from the network and ignored: ",
      paste(dropped, collapse = ", "), "."
    ))
  }

  x <- dplyr::inner_join(x, nodes, by = "trait")
  x <- x[!is.na(x$sri), , drop = FALSE]

  if (isTRUE(include_unassigned)) {
    x$label <- ifelse(is.na(x$module), "unassigned", paste0("M", x$module))
  } else {
    x <- x[!is.na(x$module), , drop = FALSE]
    x$label <- paste0("M", x$module)
  }
  if (nrow(x) == 0L) {
    ps_abort("No trait with both a module and a usable index; nothing to score.")
  }

  out <- dplyr::summarise(
    dplyr::group_by(x, .data$unit, .data$group, .data$label),
    module = .data$module[1],
    n_traits = dplyr::n(),
    hub = .data$trait[which.max(.data$strength)][1],
    score = if (aggregate == "mean") {
      mean(.data$sri)
    } else {
      sqrt(mean(.data$sri^2))
    },
    traits = paste(sort(.data$trait), collapse = ", "),
    .groups = "drop"
  )

  out <- out[c("unit", "group", "module", "label", "n_traits", "hub", "score", "traits")]
  dplyr::arrange(out, .data$unit, .data$group, dplyr::desc(.data$score))
}
