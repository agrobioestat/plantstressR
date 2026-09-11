# Build the weight vector used to collapse per-trait indices into one number.
sri_weights <- function(x, weights) {
  traits <- unique(x$trait)

  if (is.numeric(weights)) {
    if (is.null(names(weights))) {
      ps_abort("A numeric `weights` vector must be named after the traits.")
    }
    if (any(!is.finite(weights)) || any(weights < 0)) {
      ps_abort("`weights` must be finite and non-negative.")
    }
    missing_w <- setdiff(traits, names(weights))
    if (length(missing_w) > 0L) {
      ps_abort(paste0(
        "`weights` is missing entries for trait(s): ",
        paste(missing_w, collapse = ", "), "."
      ))
    }
    w <- weights[traits]
    scheme <- "manual"
  } else {
    scheme <- match.arg(weights, c("equal", "precision", "pca"))
    w <- switch(scheme,
      equal = stats::setNames(rep(1, length(traits)), traits),
      precision = {
        v <- vapply(
          traits,
          function(tr) {
            se <- x$se[x$trait == tr]
            se <- se[is.finite(se) & se > 0]
            if (length(se) == 0L) NA_real_ else 1 / mean(se^2)
          },
          numeric(1)
        )
        if (all(is.na(v))) {
          ps_warn('No usable standard errors; falling back to "equal" weights.')
          scheme <- "equal"
          v <- stats::setNames(rep(1, length(traits)), traits)
        }
        v[is.na(v)] <- 0
        v
      },
      pca = {
        mat <- impute_column_means(sri_matrix(x))
        keep <- apply(mat, 2, function(col) is.finite(stats::sd(col)) && stats::sd(col) > 0)
        if (nrow(mat) < 3L || sum(keep) < 2L) {
          ps_warn(paste0(
            'Not enough units x traits variation for PCA weights; ',
            'falling back to "equal" weights.'
          ))
          scheme <- "equal"
          stats::setNames(rep(1, length(traits)), traits)
        } else {
          pc <- stats::prcomp(mat[, keep, drop = FALSE], center = TRUE, scale. = TRUE)
          v <- stats::setNames(rep(0, length(traits)), traits)
          v[colnames(mat)[keep]] <- abs(pc$rotation[, 1])
          v
        }
      }
    )
  }

  total <- sum(w, na.rm = TRUE)
  if (!is.finite(total) || total <= 0) {
    ps_abort("Weights sum to zero; the integrated index cannot be computed.")
  }
  list(
    weights = tibble::tibble(trait = traits, weight = unname(w[traits] / total)),
    scheme = scheme
  )
}

#' Integrated Stress Index (ISI) and Ranking
#'
#' @description
#' Collapses the per-trait indices produced by [calculate_sri()] into a single
#' weighted index per experimental unit and stress level, and ranks the units.
#' Because every trait index is already signed and dimensionless, the aggregate
#' is a physiologically interpretable severity score: `0` means "indistinguishable
#' from the control", positive values mean impairment, negative values mean the
#' unit performed better under stress than under control conditions.
#'
#' @details
#' Four weighting schemes are available:
#'
#' \describe{
#'   \item{`"equal"`}{Every trait contributes the same amount. Use it when the
#'     trait panel was chosen a priori and no trait should dominate.}
#'   \item{`"precision"`}{Weight proportional to \eqn{1/\overline{SE^2}}, so that
#'     noisy traits with few replicates are down-weighted.}
#'   \item{`"pca"`}{Absolute loadings of the first principal component of the
#'     unit-by-trait index matrix; traits that carry the dominant axis of stress
#'     variation weigh more. Requires at least three units and two varying
#'     traits, and falls back to `"equal"` otherwise.}
#'   \item{a named numeric vector}{Expert weights, one per trait; internally
#'     rescaled to sum to one.}
#' }
#'
#' `aggregate = "mean"` returns the weighted mean of the indices and preserves
#' the sign. `aggregate = "rms"` returns the weighted quadratic mean, a pure
#' magnitude of physiological disturbance that does not let impaired and
#' improved traits cancel out.
#'
#' @param sri A `plantstress_sri` object from [calculate_sri()], or any data
#'   frame with the columns `unit`, `group`, `trait`, `sri` and `se`.
#' @param weights `"equal"`, `"precision"`, `"pca"`, or a named numeric vector of
#'   trait weights.
#' @param aggregate `"mean"` (signed weighted mean) or `"rms"` (weighted
#'   quadratic mean).
#' @param rank_by `"tolerance"` (rank 1 = lowest index = most tolerant, the
#'   default for selection work) or `"severity"` (rank 1 = highest index).
#'   Ranks are computed within each stress level.
#' @param rescale Logical. Add `isi_scaled`, the index linearly rescaled to
#'   `0-100` within each stress level, so that `0` marks the most tolerant and
#'   `100` the most affected unit of that level.
#' @param na_rm Logical. Drop missing trait indices instead of propagating them.
#'
#' @return An object of class `"plantstress_isi"`, a tibble with columns `unit`,
#'   `group`, `isi`, `isi_scaled`, `n_traits` and `rank`. The trait weights and
#'   the per-trait contributions are stored as attributes and returned by
#'   [stress_weights()] and [stress_contributions()].
#'
#' @references
#' Rosielle A.A., Hamblin J. (1981). Theoretical aspects of selection for yield
#' in stress and non-stress environments. \doi{10.2135/cropsci1981.0011183X002100060033x}
#'
#' @seealso [calculate_sri()], [plot_stress_signature()]
#'
#' @examples
#' data(brachiaria_stress)
#'
#' sri <- calculate_sri(
#'   brachiaria_stress,
#'   treatment = "drought_level",
#'   control = "control",
#'   traits = c("Fv_Fm", "PIabs", "A", "gs", "RWC", "shoot_biomass"),
#'   by = "genotype",
#'   verbose = FALSE
#' )
#'
#' isi <- integrated_stress_index(sri, weights = "precision")
#' isi
#'
#' stress_weights(isi)
#' head(stress_contributions(isi))
#' @export
integrated_stress_index <- function(sri,
                                    weights = c("equal", "precision", "pca"),
                                    aggregate = c("mean", "rms"),
                                    rank_by = c("tolerance", "severity"),
                                    rescale = TRUE,
                                    na_rm = TRUE) {
  aggregate <- match.arg(aggregate)
  rank_by <- match.arg(rank_by)

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
  if (!"se" %in% names(sri)) sri$se <- NA_real_

  x <- tibble::as_tibble(sri)
  if (isTRUE(na_rm)) {
    x <- x[!is.na(x$sri), , drop = FALSE]
  }
  if (nrow(x) == 0L) {
    ps_abort("No usable stress response indices are available.")
  }

  w <- sri_weights(x, weights)
  x <- dplyr::left_join(x, w$weights, by = "trait")

  x$contribution <- if (aggregate == "mean") {
    x$weight * x$sri
  } else {
    x$weight * x$sri^2
  }

  isi <- dplyr::summarise(
    dplyr::group_by(x, .data$unit, .data$group),
    n_traits = dplyr::n(),
    total_weight = sum(.data$weight, na.rm = TRUE),
    value = sum(.data$contribution, na.rm = TRUE),
    .groups = "drop"
  )
  isi$isi <- if (aggregate == "mean") {
    isi$value / isi$total_weight
  } else {
    sqrt(isi$value / isi$total_weight)
  }

  # Ranking and rescaling happen *within* a stress level: a moderate and a
  # severe treatment are different environments and are not comparable on one
  # scale.
  isi$isi_scaled <- NA_real_
  if (isTRUE(rescale)) {
    isi$isi_scaled <- stats::ave(isi$isi, isi$group, FUN = function(v) {
      rng <- range(v, na.rm = TRUE)
      if (length(v) < 2L || !is.finite(diff(rng)) || diff(rng) <= 0) {
        return(rep(NA_real_, length(v)))
      }
      100 * (v - rng[1]) / diff(rng)
    })
    if (all(is.na(isi$isi_scaled))) {
      ps_warn("`isi_scaled` needs at least two differing units per stress level; returning NA.")
    }
  }

  sign_rank <- if (rank_by == "tolerance") 1 else -1
  isi$rank <- stats::ave(isi$isi, isi$group, FUN = function(v) {
    rank(sign_rank * v, ties.method = "min", na.last = "keep")
  })

  contributions <- x[c("unit", "group", "trait", "sri", "weight", "contribution")]
  contributions <- dplyr::arrange(
    contributions, .data$unit, .data$group, dplyr::desc(abs(.data$contribution))
  )

  out <- dplyr::arrange(
    isi[c("unit", "group", "isi", "isi_scaled", "n_traits", "rank")],
    .data$group, .data$rank
  )

  structure(
    out,
    class = c("plantstress_isi", class(tibble::tibble())),
    plantstress = list(
      weights = w$weights,
      weighting = w$scheme,
      aggregate = aggregate,
      rank_by = rank_by,
      contributions = contributions
    )
  )
}

#' @param x A `plantstress_isi` object.
#' @param ... Ignored.
#' @rdname integrated_stress_index
#' @export
print.plantstress_isi <- function(x, ...) {
  meta <- sri_meta(x)
  cat("<plantstress_isi>\n")
  cat("  Weighting: ", meta$weighting %||% "?", "\n", sep = "")
  cat("  Aggregate: ", meta$aggregate %||% "?", "\n", sep = "")
  cat("  Ranking:   ", meta$rank_by %||% "?",
    " (rank 1 = ", if (identical(meta$rank_by, "severity")) "most stressed" else "most tolerant",
    ")\n",
    sep = ""
  )
  print(tibble::as_tibble(x))
  invisible(x)
}

#' Weights and Trait Contributions of an Integrated Stress Index
#'
#' @description
#' Accessors for the two by-products of [integrated_stress_index()]: the trait
#' weights actually used, and how much each trait contributed to each unit's
#' index. Both are the natural inputs for reporting why a genotype was ranked
#' where it was.
#'
#' @param x A `plantstress_isi` object.
#'
#' @return A tibble. `stress_weights()` has columns `trait` and `weight` (summing
#'   to one); `stress_contributions()` has columns `unit`, `group`, `trait`,
#'   `sri`, `weight` and `contribution`.
#'
#' @seealso [integrated_stress_index()]
#'
#' @examples
#' data(brachiaria_stress)
#' sri <- calculate_sri(brachiaria_stress,
#'   treatment = "drought_level", control = "control",
#'   traits = c("Fv_Fm", "A", "RWC"), by = "genotype", verbose = FALSE
#' )
#' isi <- integrated_stress_index(sri)
#' stress_weights(isi)
#' @export
stress_weights <- function(x) {
  if (!inherits(x, "plantstress_isi")) {
    ps_abort("`x` must be a `plantstress_isi` object from `integrated_stress_index()`.")
  }
  sri_meta(x)$weights
}

#' @rdname stress_weights
#' @export
stress_contributions <- function(x) {
  if (!inherits(x, "plantstress_isi")) {
    ps_abort("`x` must be a `plantstress_isi` object from `integrated_stress_index()`.")
  }
  sri_meta(x)$contributions
}
