#' Stability of the Stress Response Across Environments
#'
#' @description
#' Summarizes a multi-environment trial: how severely each genotype was
#' affected on average, and how consistently. A genotype that ranks first at one
#' site and last at the next is not the same proposition as one that ranks
#' second everywhere, and a single mean hides the difference.
#'
#' Run [calculate_sri()] with `by = c("genotype", "site")` -- or whatever names
#' the genotype and the environment -- then [integrated_stress_index()], then
#' this function.
#'
#' @details
#' The table is read per stress level, because a moderate and a severe
#' treatment are different environments in themselves and are never pooled.
#'
#' `sd_isi` and `cv_isi` describe the spread of a genotype's index across
#' environments. `rank_min`, `rank_max` and `mean_rank` describe the spread of
#' its position, which is what selection actually acts on.
#'
#' `ecovalence` is Wricke's contribution to the genotype-by-environment
#' interaction,
#' \deqn{W_i = \sum_j (X_{ij} - \bar{X}_{i.} - \bar{X}_{.j} + \bar{X}_{..})^2,}
#' computed on the integrated index rather than on yield. It is zero for a
#' genotype whose response tracks the environment mean exactly, and large for
#' one that reshuffles. `ecovalence_pct` expresses it as a share of the total
#' interaction, so the genotypes responsible for the instability of the trial
#' can be named. Both need at least two environments and two genotypes, and are
#' `NA` otherwise.
#'
#' Genotypes missing from some environments are kept, with `n_env` reporting how
#' many they appeared in, but their ecovalence is `NA`: an unbalanced
#' genotype-by-environment table has no meaningful interaction decomposition.
#'
#' @param isi A `plantstress_isi` object from [integrated_stress_index()],
#'   computed from an SRI whose `by` named at least two columns.
#' @param environment Which of the `by` columns is the environment. Defaults to
#'   the last one. The remaining columns identify the genotype.
#' @param rank_by `"tolerance"` (rank 1 = least affected, the default) or
#'   `"severity"`. Ranks are recomputed within each environment and stress
#'   level, so they describe the standing of a genotype among the material it
#'   was actually grown beside.
#'
#' @return An object of class `"plantstress_stability"`, a tibble with one row
#'   per genotype and stress level and columns `unit`, `group`, `n_env`,
#'   `mean_isi`, `sd_isi`, `cv_isi`, `mean_rank`, `rank_min`, `rank_max`,
#'   `ecovalence` and `ecovalence_pct`, ordered by `group` and `mean_isi`.
#'
#' @references
#' Wricke G. (1962). Uber eine Methode zur Erfassung der okologischen Streubreite
#' in Feldversuchen. Zeitschrift fur Pflanzenzuchtung 47, 92-96.
#'
#' Finlay K.W., Wilkinson G.N. (1963). The analysis of adaptation in a plant
#' breeding programme. \doi{10.1071/AR9630742}
#'
#' @seealso [integrated_stress_index()], [calculate_sri()]
#'
#' @examples
#' data(brachiaria_stress)
#'
#' # Two sites, by splitting the blocks of the bundled trial.
#' met <- brachiaria_stress
#' met$site <- ifelse(met$block %in% c("B1", "B2"), "north", "south")
#'
#' sri <- calculate_sri(met,
#'   treatment = "drought_level", control = "control",
#'   traits = c("Fv_Fm", "A", "RWC", "shoot_biomass"),
#'   by = c("genotype", "site"), verbose = FALSE
#' )
#'
#' stress_stability(integrated_stress_index(sri))
#' @export
stress_stability <- function(isi,
                             environment = NULL,
                             rank_by = c("tolerance", "severity")) {
  rank_by <- match.arg(rank_by)
  if (!inherits(isi, "plantstress_isi")) {
    ps_abort("`isi` must be a `plantstress_isi` object from `integrated_stress_index()`.")
  }

  by <- sri_meta(isi)$by
  if (is.null(by) || length(by) < 2L) {
    ps_abort(paste0(
      "Stability needs at least two grouping columns. Recompute the index with ",
      "`calculate_sri(by = c(\"genotype\", \"site\"))`, naming the genotype and ",
      "the environment."
    ))
  }

  if (is.null(environment)) {
    environment <- by[length(by)]
  }
  check_string(environment, "environment")
  if (!environment %in% by) {
    ps_abort(paste0(
      "`environment` must be one of the grouping columns (",
      paste(by, collapse = ", "), "); got `", environment, "`."
    ))
  }
  geno_cols <- setdiff(by, environment)

  df <- tibble::as_tibble(isi)
  parts <- split_unit_labels(df$unit, by)
  df$.env <- parts[[environment]]
  df$.geno <- if (length(geno_cols) == 1L) {
    parts[[geno_cols]]
  } else {
    do.call(paste, c(unname(as.list(parts[geno_cols])), list(sep = PS_UNIT_SEP)))
  }
  df <- df[!is.na(df$isi), , drop = FALSE]
  if (nrow(df) == 0L) {
    ps_abort("No usable integrated indices are available.")
  }
  if (length(unique(df$.env)) < 2L) {
    ps_warn(paste0(
      "Only one level of `", environment, "` is present; the spread across ",
      "environments cannot be estimated."
    ))
  }

  # Standing within the material each genotype was actually grown beside.
  sign_rank <- if (rank_by == "tolerance") 1 else -1
  df$.rank <- stats::ave(
    df$isi, paste(df$group, df$.env),
    FUN = function(v) rank(sign_rank * v, ties.method = "min", na.last = "keep")
  )

  eco <- ecovalence_by_group(df)

  out <- dplyr::summarise(
    dplyr::group_by(df, .data$group, unit = .data$.geno),
    n_env = dplyr::n_distinct(.data$.env),
    mean_isi = mean(.data$isi),
    sd_isi = if (dplyr::n() > 1L) stats::sd(.data$isi) else NA_real_,
    mean_rank = mean(.data$.rank),
    rank_min = min(.data$.rank),
    rank_max = max(.data$.rank),
    .groups = "drop"
  )
  out$cv_isi <- ifelse(
    is.finite(out$mean_isi) & abs(out$mean_isi) > .Machine$double.eps^0.5,
    100 * out$sd_isi / abs(out$mean_isi),
    NA_real_
  )
  out <- dplyr::left_join(out, eco, by = c("group", "unit"))

  out <- out[c(
    "unit", "group", "n_env", "mean_isi", "sd_isi", "cv_isi",
    "mean_rank", "rank_min", "rank_max", "ecovalence", "ecovalence_pct"
  )]
  out <- dplyr::arrange(out, .data$group, .data$mean_isi)

  structure(
    out,
    class = c("plantstress_stability", class(tibble::tibble())),
    plantstress = list(
      by = by,
      environment = environment,
      genotype = geno_cols,
      rank_by = rank_by,
      n_env = length(unique(df$.env))
    )
  )
}

# Wricke's ecovalence per stress level, on the genotype-by-environment table of
# integrated indices. Returns NA for every genotype of a level whose table is
# not complete, since the decomposition assumes it is.
ecovalence_by_group <- function(df) {
  out <- list()
  for (g in unique(df$group)) {
    sub <- df[df$group == g, , drop = FALSE]
    genos <- unique(sub$.geno)
    envs <- unique(sub$.env)
    empty <- tibble::tibble(
      group = g, unit = genos,
      ecovalence = NA_real_, ecovalence_pct = NA_real_
    )

    if (length(genos) < 2L || length(envs) < 2L ||
      nrow(sub) != length(genos) * length(envs) ||
      anyDuplicated(paste(sub$.geno, sub$.env)) > 0L) {
      out[[length(out) + 1L]] <- empty
      next
    }

    mat <- matrix(NA_real_,
      nrow = length(genos), ncol = length(envs),
      dimnames = list(genos, envs)
    )
    mat[cbind(match(sub$.geno, genos), match(sub$.env, envs))] <- sub$isi
    if (anyNA(mat)) {
      out[[length(out) + 1L]] <- empty
      next
    }

    resid <- mat - rowMeans(mat)[row(mat)] - colMeans(mat)[col(mat)] + mean(mat)
    w <- rowSums(resid^2)
    total <- sum(w)
    out[[length(out) + 1L]] <- tibble::tibble(
      group = g, unit = genos,
      ecovalence = unname(w),
      ecovalence_pct = if (total > 0) unname(100 * w / total) else NA_real_
    )
  }
  dplyr::bind_rows(out)
}

#' @param x A `plantstress_stability` object.
#' @param ... Ignored.
#' @rdname stress_stability
#' @export
print.plantstress_stability <- function(x, ...) {
  meta <- sri_meta(x)
  cat("<plantstress_stability>\n")
  cat("  Genotype:     ", paste(meta$genotype %||% "?", collapse = ", "), "\n", sep = "")
  cat("  Environment:  ", meta$environment %||% "?",
    " (", meta$n_env %||% "?", " levels)\n",
    sep = ""
  )
  cat("  Ranking:      ", meta$rank_by %||% "?",
    " (rank 1 = ", if (identical(meta$rank_by, "severity")) "most stressed" else "most tolerant",
    ", within each environment)\n",
    sep = ""
  )
  cat("  ecovalence = share of the genotype x environment interaction\n")
  print(tibble::as_tibble(x))
  invisible(x)
}

#' @rdname stress_stability
#' @export
plot.plantstress_stability <- function(x, ...) {
  df <- tibble::as_tibble(x)
  df <- df[!is.na(df$sd_isi), , drop = FALSE]
  if (nrow(df) == 0L) {
    ps_abort("Nothing to draw: the spread across environments is undefined.")
  }

  ggplot2::ggplot(df, ggplot2::aes(x = .data$mean_isi, y = .data$sd_isi)) +
    ggplot2::geom_point(ggplot2::aes(size = .data$n_env),
      shape = 21, fill = "#2C7BB6", colour = "grey25", alpha = 0.85
    ) +
    ggplot2::geom_text(ggplot2::aes(label = .data$unit),
      vjust = -1.1, size = 3.1, colour = "grey30"
    ) +
    ggplot2::facet_wrap(~ .data$group) +
    ggplot2::scale_size_continuous(range = c(2.5, 5), name = "environments") +
    ggplot2::labs(
      x = "Mean integrated stress index across environments",
      y = "Standard deviation across environments",
      title = "Severity against consistency",
      subtitle = "Bottom left = least affected and most predictable"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"))
}
