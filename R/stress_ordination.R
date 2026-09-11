#' Multivariate Ordination of the Trait Space
#'
#' @description
#' Places individual plants in the multivariate trait space with a principal
#' component analysis, so that the separation between control and stress
#' treatments can be inspected at the sample level. It complements
#' [calculate_sri()], which summarizes treatments, and [stress_network()], which
#' summarizes traits.
#'
#' @details
#' Traits are centered and, by default, scaled to unit variance, which is
#' mandatory whenever the panel mixes units. Missing values are imputed with the
#' trait mean before decomposition; traits with zero variance are dropped with a
#' warning because they carry no ordination information.
#'
#' The `"stats"` engine uses [stats::prcomp()] and has no dependency beyond base
#' R. The `"FactoMineR"` engine calls `FactoMineR::PCA()` when that package is
#' installed, which is convenient if the rest of the analysis is being done with
#' the FactoMineR and factoextra toolchain; results are returned in the same
#' shape either way, and the engine silently falls back to `"stats"` when
#' FactoMineR is not available.
#'
#' @param data A data frame with one row per experimental unit.
#' @param treatment Optional name of a column used to color and group samples.
#' @param traits Character vector of trait columns. Defaults to every numeric
#'   column that is not a design column.
#' @param block Optional name of a block (replicate) column whose additive
#'   effect is removed from every trait before the decomposition, so that the
#'   first components describe treatment and genotype rather than the layout of
#'   the trial.
#' @param ncomp Number of components to retain.
#' @param scale Logical. Scale traits to unit variance.
#' @param engine `"stats"` (default) or `"FactoMineR"`.
#'
#' @return An object of class `"plantstress_ordination"`: a list with `scores`,
#'   `loadings`, `eigenvalues`, `traits`, `treatment` and `engine`.
#'
#' @references
#' Le S., Josse J., Husson F. (2008). FactoMineR: an R package for multivariate
#' analysis. \doi{10.18637/jss.v025.i01}
#'
#' @seealso [calculate_sri()], [stress_network()]
#'
#' @examples
#' data(brachiaria_stress)
#' ord <- stress_ordination(
#'   brachiaria_stress,
#'   treatment = "drought_level",
#'   traits = c("Fv_Fm", "PIabs", "A", "gs", "RWC", "SPAD", "shoot_biomass")
#' )
#' ord
#' plot(ord)
#' @export
stress_ordination <- function(data,
                              treatment = NULL,
                              traits = NULL,
                              block = NULL,
                              ncomp = 2L,
                              scale = TRUE,
                              engine = c("stats", "FactoMineR")) {
  engine <- match.arg(engine)
  check_data(data)
  if (!is.null(treatment)) check_column(data, treatment, "treatment")
  if (!is.null(block)) check_column(data, block, "block")
  if (!is.numeric(ncomp) || length(ncomp) != 1L || is.na(ncomp) || ncomp < 2) {
    ps_abort("`ncomp` must be a numeric scalar >= 2.")
  }
  ncomp <- as.integer(ncomp)

  traits <- resolve_traits(data, traits, exclude = c(treatment, block))
  if (!is.null(block)) {
    data <- block_adjust_traits(data, traits, treatment, block)
  }
  mat <- impute_column_means(as.matrix(data[traits]))

  sds <- apply(mat, 2, stats::sd)
  constant <- !is.finite(sds) | sds == 0
  if (any(constant)) {
    ps_warn(paste0(
      "Dropping trait(s) with zero variance: ",
      paste(traits[constant], collapse = ", "), "."
    ))
    mat <- mat[, !constant, drop = FALSE]
    traits <- traits[!constant]
  }
  if (ncol(mat) < 2L) {
    ps_abort("At least two varying traits are required for an ordination.")
  }
  if (nrow(mat) < 3L) {
    ps_abort("At least three observations are required for an ordination.")
  }

  max_comp <- min(ncomp, ncol(mat), nrow(mat) - 1L)
  if (max_comp < ncomp) {
    ps_warn(paste0("`ncomp` reduced to ", max_comp, " given the data dimensions."))
  }

  if (engine == "FactoMineR" && !requireNamespace("FactoMineR", quietly = TRUE)) {
    ps_warn('Package "FactoMineR" is not installed; using the "stats" engine.')
    engine <- "stats"
  }

  if (engine == "FactoMineR") {
    fit <- FactoMineR::PCA(
      as.data.frame(mat),
      scale.unit = scale, ncp = max_comp, graph = FALSE
    )
    scores <- fit$ind$coord[, seq_len(max_comp), drop = FALSE]
    loadings <- fit$var$coord[, seq_len(max_comp), drop = FALSE]
    eig <- fit$eig[seq_len(max_comp), , drop = FALSE]
    eigenvalues <- tibble::tibble(
      component = paste0("PC", seq_len(max_comp)),
      eigenvalue = as.numeric(eig[, 1]),
      variance_pct = as.numeric(eig[, 2]),
      cumulative_pct = as.numeric(eig[, 3])
    )
  } else {
    fit <- stats::prcomp(mat, center = TRUE, scale. = scale)
    scores <- fit$x[, seq_len(max_comp), drop = FALSE]
    loadings <- fit$rotation[, seq_len(max_comp), drop = FALSE]
    variance <- fit$sdev^2
    eigenvalues <- tibble::tibble(
      component = paste0("PC", seq_len(max_comp)),
      eigenvalue = variance[seq_len(max_comp)],
      variance_pct = 100 * variance[seq_len(max_comp)] / sum(variance),
      cumulative_pct = cumsum(100 * variance / sum(variance))[seq_len(max_comp)]
    )
  }

  colnames(scores) <- paste0("PC", seq_len(max_comp))
  colnames(loadings) <- paste0("PC", seq_len(max_comp))

  scores_tbl <- tibble::as_tibble(scores)
  scores_tbl$group <- if (is.null(treatment)) "all" else as.character(data[[treatment]])

  loadings_tbl <- tibble::as_tibble(loadings)
  loadings_tbl$trait <- rownames(loadings)
  loadings_tbl <- loadings_tbl[c("trait", setdiff(names(loadings_tbl), "trait"))]

  out <- list(
    scores = scores_tbl,
    loadings = loadings_tbl,
    eigenvalues = eigenvalues,
    traits = traits,
    treatment = treatment,
    engine = engine,
    scaled = scale
  )
  class(out) <- "plantstress_ordination"
  out
}

#' @param x A `plantstress_ordination` object.
#' @param ... Ignored.
#' @rdname stress_ordination
#' @export
print.plantstress_ordination <- function(x, ...) {
  cat("<plantstress_ordination>\n")
  cat("  Engine:     ", x$engine, "\n", sep = "")
  cat("  Samples:    ", nrow(x$scores), "\n", sep = "")
  cat("  Traits:     ", length(x$traits), "\n", sep = "")
  cat("  Treatment:  ", x$treatment %||% "<none>", "\n", sep = "")
  cat("  Variance explained:\n")
  for (i in seq_len(nrow(x$eigenvalues))) {
    cat("   - ", x$eigenvalues$component[i], ": ",
      sprintf("%.1f%%", x$eigenvalues$variance_pct[i]),
      " (cumulative ", sprintf("%.1f%%", x$eigenvalues$cumulative_pct[i]), ")\n",
      sep = ""
    )
  }
  invisible(x)
}

#' @param components Integer vector of length two giving the components to draw.
#' @param loadings Logical. Draw trait loading arrows.
#' @param ellipse Logical. Draw a normal confidence ellipse per group.
#' @rdname stress_ordination
#' @export
plot.plantstress_ordination <- function(x,
                                        components = c(1L, 2L),
                                        loadings = TRUE,
                                        ellipse = TRUE,
                                        ...) {
  if (!is.numeric(components) || length(components) != 2L || anyNA(components)) {
    ps_abort("`components` must be a numeric vector of length two.")
  }
  comp <- paste0("PC", as.integer(components))
  missing_comp <- setdiff(comp, names(x$scores))
  if (length(missing_comp) > 0L) {
    ps_abort(paste0(
      "Component(s) not available: ", paste(missing_comp, collapse = ", "),
      ". Refit with a larger `ncomp`."
    ))
  }

  scores <- x$scores
  scores$.x <- scores[[comp[1]]]
  scores$.y <- scores[[comp[2]]]

  idx <- match(comp, x$eigenvalues$component)
  lab <- sprintf("%s (%.1f%%)", comp, x$eigenvalues$variance_pct[idx])

  p <- ggplot2::ggplot(
    scores,
    ggplot2::aes(x = .data$.x, y = .data$.y, colour = .data$group, fill = .data$group)
  ) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey85") +
    ggplot2::geom_hline(yintercept = 0, colour = "grey85") +
    ggplot2::geom_point(size = 2, alpha = 0.8, na.rm = TRUE)

  if (isTRUE(ellipse) && min(table(scores$group)) >= 4L) {
    p <- p + ggplot2::stat_ellipse(
      ggplot2::aes(group = .data$group), type = "norm", alpha = 0.25,
      geom = "polygon", show.legend = FALSE
    )
  }

  if (isTRUE(loadings)) {
    scale_factor <- 0.9 * max(abs(c(scores$.x, scores$.y)), na.rm = TRUE) /
      max(abs(c(x$loadings[[comp[1]]], x$loadings[[comp[2]]])), na.rm = TRUE)
    ld <- x$loadings
    ld$.x <- ld[[comp[1]]] * scale_factor
    ld$.y <- ld[[comp[2]]] * scale_factor
    p <- p +
      ggplot2::geom_segment(
        data = ld,
        ggplot2::aes(x = 0, y = 0, xend = .data$.x, yend = .data$.y),
        inherit.aes = FALSE, colour = "grey30",
        arrow = ggplot2::arrow(length = ggplot2::unit(0.15, "cm"))
      ) +
      ggplot2::geom_text(
        data = ld,
        ggplot2::aes(x = .data$.x, y = .data$.y, label = .data$trait),
        inherit.aes = FALSE, colour = "grey20", size = 3, vjust = -0.6,
        fontface = "italic"
      )
  }

  p +
    ggplot2::labs(
      x = lab[1], y = lab[2],
      colour = x$treatment %||% "Group", fill = x$treatment %||% "Group",
      title = "Ordination of the physiological trait space"
    ) +
    ggplot2::theme_minimal(base_size = 11)
}
