# Resample a trial while preserving its layout: whole blocks when the design is
# blocked, individual plants inside each design cell otherwise. Duplicated
# blocks are relabelled, because a block drawn twice is two replicates of the
# same stratum, not one stratum with twice the plants.
resample_trial <- function(data, treatment, by, block) {
  unit <- unit_labels(data, by)
  idx <- integer(0)
  new_block <- character(0)

  if (is.null(block)) {
    key <- paste(unit, as.character(data[[treatment]]), sep = "\r")
    for (k in unique(key)) {
      rows <- which(key == k)
      idx <- c(idx, sample(rows, length(rows), replace = TRUE))
    }
  } else {
    blk <- as.character(data[[block]])
    for (u in unique(unit)) {
      rows_u <- which(unit == u)
      blocks_u <- unique(stats::na.omit(blk[rows_u]))
      if (length(blocks_u) < 2L) {
        idx <- c(idx, rows_u)
        new_block <- c(new_block, blk[rows_u])
        next
      }
      picked <- sample(blocks_u, length(blocks_u), replace = TRUE)
      for (i in seq_along(picked)) {
        rows_i <- rows_u[!is.na(blk[rows_u]) & blk[rows_u] == picked[i]]
        idx <- c(idx, rows_i)
        new_block <- c(new_block, rep(paste0("rb", i), length(rows_i)))
      }
    }
  }

  out <- data[idx, , drop = FALSE]
  if (!is.null(block)) {
    out[[block]] <- new_block
  }
  out
}

# Difference between each stress level and the control for one trait of one
# unit: adjusted means from the additive model when the trial is blocked, plain
# group means otherwise.
cell_differences <- function(v, trt, blk, control, stress_levels) {
  out <- stats::setNames(rep(NA_real_, length(stress_levels)), stress_levels)
  ok <- !is.na(v) & !is.na(trt)
  if (!is.null(blk)) ok <- ok & !is.na(blk)
  if (sum(ok) < 3L) {
    return(out)
  }
  v <- v[ok]
  trt <- as.character(trt[ok])
  if (!control %in% trt) {
    return(out)
  }
  present <- intersect(stress_levels, unique(trt))
  if (length(present) == 0L) {
    return(out)
  }

  if (is.null(blk)) {
    m <- tapply(v, trt, mean)
    for (g in present) out[[g]] <- unname(m[[g]] - m[[control]])
    return(out)
  }

  blk <- as.character(blk[ok])
  keep <- trt %in% c(control, present)
  f_trt <- stats::relevel(factor(trt[keep]), ref = control)
  f_blk <- factor(blk[keep])
  fit <- if (nlevels(f_blk) < 2L) {
    try(stats::lm(v[keep] ~ f_trt), silent = TRUE)
  } else {
    try(stats::lm(v[keep] ~ f_trt + f_blk), silent = TRUE)
  }
  if (inherits(fit, "try-error")) {
    return(out)
  }
  cf <- stats::coef(fit)
  for (g in present) {
    slot <- paste0("f_trt", g)
    if (slot %in% names(cf) && is.finite(cf[[slot]])) out[[g]] <- unname(cf[[slot]])
  }
  out
}

# Every control-versus-stress difference of a resampled trial, in one table.
boot_differences <- function(data, meta) {
  trt <- as.character(data[[meta$treatment]])
  unit <- unit_labels(data, meta$by)
  blk <- if (is.null(meta$block)) NULL else as.character(data[[meta$block]])
  stress_levels <- setdiff(unique(stats::na.omit(trt)), meta$control)
  if (length(stress_levels) == 0L) {
    return(NULL)
  }

  units <- unique(stats::na.omit(unit))
  n_out <- length(units) * length(stress_levels) * length(meta$traits)
  out_unit <- character(n_out)
  out_group <- character(n_out)
  out_trait <- character(n_out)
  out_diff <- rep(NA_real_, n_out)
  k <- 0L

  for (u in units) {
    in_u <- !is.na(unit) & unit == u
    for (tr in meta$traits) {
      diffs <- cell_differences(
        data[[tr]][in_u], trt[in_u],
        if (is.null(blk)) NULL else blk[in_u],
        meta$control, stress_levels
      )
      for (g in stress_levels) {
        k <- k + 1L
        out_unit[k] <- u
        out_group[k] <- g
        out_trait[k] <- tr
        out_diff[k] <- diffs[[g]]
      }
    }
  }

  tibble::tibble(
    unit = out_unit, group = out_group, trait = out_trait, difference = out_diff
  )
}

# One bootstrap replicate. With `ruler` supplied the observed scaling is reused
# and only the response is resampled; otherwise the whole chain is recomputed.
boot_replicate <- function(meta, weights, aggregate, rank_by, ruler = NULL) {
  boot_data <- resample_trial(meta$data, meta$treatment, meta$by, meta$block)

  if (is.null(ruler)) {
    sri_b <- try(
      suppressWarnings(calculate_sri(
        boot_data,
        treatment = meta$treatment, control = meta$control, traits = meta$traits,
        by = meta$by, block = meta$block, direction = meta$direction,
        method = meta$method, conf_level = meta$conf_level,
        p_adjust = meta$p_adjust, verbose = FALSE
      )),
      silent = TRUE
    )
    if (inherits(sri_b, "try-error")) {
      return(NULL)
    }
  } else {
    diffs <- try(boot_differences(boot_data, meta), silent = TRUE)
    if (inherits(diffs, "try-error") || is.null(diffs)) {
      return(NULL)
    }
    pos <- match(
      paste(diffs$unit, diffs$group, diffs$trait, sep = "\r"), ruler$key
    )
    diffs$sri <- ruler$sign[pos] * diffs$difference / ruler$scale[pos]
    diffs$se_sampling <- ruler$se[pos]
    sri_b <- diffs[!is.na(diffs$sri), , drop = FALSE]
    if (nrow(sri_b) == 0L) {
      return(NULL)
    }
  }

  isi_b <- try(
    suppressWarnings(integrated_stress_index(
      sri_b,
      weights = weights, aggregate = aggregate,
      rank_by = rank_by, rescale = FALSE
    )),
    silent = TRUE
  )
  if (inherits(isi_b, "try-error")) {
    return(NULL)
  }
  tibble::as_tibble(isi_b)[c("unit", "group", "isi", "rank")]
}

#' Bootstrap Confidence Limits and Rank Stability for the Integrated Index
#'
#' @description
#' Puts an uncertainty statement on the ranking. A table of genotypes ordered by
#' their Integrated Stress Index looks decisive, but the order is an estimate
#' like any other: rerun the trial and it would come out somewhat differently.
#' `stress_index_ci()` resamples the trial, recomputes the whole chain --
#' indices, weights, aggregation and ranking -- on each resample, and reports
#' how much of the ordering survives.
#'
#' @details
#' The resampling respects the layout of the experiment. In a blocked trial
#' whole blocks are drawn with replacement, since the block is the unit that was
#' randomized; otherwise plants are drawn with replacement inside each
#' unit-by-treatment cell, which keeps the replication of every cell intact.
#'
#' `scale` decides whether the ruler is resampled along with what it measures.
#' A standardized index is a difference divided by a standard deviation, and
#' that denominator is itself an estimate: resample a handful of control plants
#' and it can collapse, which sends the index to implausible values and drags
#' the upper limit of the interval with it. On the bundled trial the smallest
#' control standard deviation is `0.0027` for `Fv/Fm` and falls to `0.0006` in
#' the worst resample, a fivefold shrinkage of the ruler that inflates every
#' index measured against it.
#'
#' With `scale = "fixed"`, the default, each trait keeps the scaling it had in
#' the observed trial, together with the replication-based precision weights, so
#' the interval reflects uncertainty about the plants' response rather than
#' about the yardstick. With `scale = "resampled"` the whole chain is recomputed
#' from scratch on every draw; that is the stricter reading of the bootstrap,
#' but expect long upper tails whenever a trait has a small control variance.
#' `"pca"` weights are recomputed from the resampled indices either way.
#'
#' `p_best` is the proportion of resamples in which a unit came out ranked
#' first. Read it as the evidence that this genotype, and not its neighbour in
#' the table, is the most tolerant one of that stress level: a `p_best` of 0.30
#' at the top of a ranking means the winner is barely distinguishable from the
#' rest.
#'
#' Replicates in which the resampled trial cannot be analysed at all -- a design
#' cell that lost its control plants, say -- are discarded, and `n_ok` reports
#' how many of the `n_boot` attempts contributed.
#'
#' @param sri A `plantstress_sri` object from [calculate_sri()]. It carries the
#'   trial and the settings it was computed with, which is what makes the
#'   resampling possible.
#' @param n_boot Number of bootstrap resamples. The default of 500 is enough for
#'   the interval; a few thousand steady the extreme quantiles. Cost grows
#'   linearly, and a blocked design fits one model per trait and unit in every
#'   replicate.
#' @param weights,aggregate,rank_by Passed to [integrated_stress_index()] and
#'   applied identically to the observed trial and to every resample.
#' @param scale `"fixed"` (default) keeps each trait's scaling standard
#'   deviation at its observed value, so that only the response is resampled;
#'   `"resampled"` re-estimates it in every draw. See Details.
#' @param conf_level Confidence level for the percentile interval.
#' @param seed Optional seed. The caller's random stream is restored afterwards.
#'
#' @return An object of class `"plantstress_isi_ci"`, a tibble with one row per
#'   unit and stress level and columns `unit`, `group`, `isi`, `conf_low`,
#'   `conf_high`, `rank`, `rank_low`, `rank_high`, `p_best` and `n_ok`. `rank`
#'   is the observed rank; `rank_low` and `rank_high` bracket the ranks the unit
#'   took across the resamples.
#'
#' @references
#' Efron B., Tibshirani R.J. (1993). An Introduction to the Bootstrap. Chapman
#' and Hall.
#'
#' @seealso [integrated_stress_index()], [calculate_sri()]
#'
#' @examples
#' data(brachiaria_stress)
#'
#' sri <- calculate_sri(
#'   brachiaria_stress,
#'   treatment = "drought_level",
#'   control = "control",
#'   traits = c("Fv_Fm", "A", "RWC", "shoot_biomass"),
#'   by = "genotype",
#'   block = "block",
#'   verbose = FALSE
#' )
#'
#' # n_boot is kept small here so the example runs quickly.
#' ci <- stress_index_ci(sri, n_boot = 25, seed = 1)
#' ci
#' @export
stress_index_ci <- function(sri,
                            n_boot = 500,
                            weights = c("equal", "precision", "pca"),
                            aggregate = c("mean", "rms"),
                            rank_by = c("tolerance", "severity"),
                            scale = c("fixed", "resampled"),
                            conf_level = 0.95,
                            seed = NULL) {
  if (!inherits(sri, "plantstress_sri")) {
    ps_abort("`sri` must be a `plantstress_sri` object from `calculate_sri()`.")
  }
  if (!is.numeric(n_boot) || length(n_boot) != 1L || is.na(n_boot) || n_boot < 2) {
    ps_abort("`n_boot` must be a numeric scalar >= 2.")
  }
  n_boot <- as.integer(n_boot)
  if (!is.numeric(weights)) weights <- match.arg(weights)
  aggregate <- match.arg(aggregate)
  rank_by <- match.arg(rank_by)
  scale <- match.arg(scale)
  check_prob(conf_level, "conf_level", lower = 0.5, upper = 0.9999)

  meta <- sri_meta(sri)
  if (is.null(meta$data)) {
    ps_abort(paste0(
      "This `plantstress_sri` object does not carry its trial, so it cannot be ",
      "resampled. Recompute it with `calculate_sri()` from this version of the ",
      "package."
    ))
  }

  observed <- integrated_stress_index(
    sri,
    weights = weights, aggregate = aggregate,
    rank_by = rank_by, rescale = FALSE
  )

  # The ruler: the scaling each trait was measured against in the observed
  # trial, recovered from `difference / effect` so that Hedges' correction and
  # the "relative" denominator come along with it.
  ruler <- NULL
  if (identical(scale, "fixed")) {
    obs <- tibble::as_tibble(sri)
    scale_obs <- obs$difference / obs$effect
    scale_obs[!is.finite(scale_obs) | scale_obs == 0] <- NA_real_
    ruler <- list(
      key = paste(obs$unit, obs$group, obs$trait, sep = "\r"),
      scale = scale_obs,
      sign = ifelse(obs$direction == "lower_is_better", 1, -1),
      se = obs$se_sampling
    )
  }

  ps_local_seed(seed)
  draws <- vector("list", n_boot)
  for (b in seq_len(n_boot)) {
    draws[[b]] <- boot_replicate(meta, weights, aggregate, rank_by, ruler)
  }
  draws <- draws[!vapply(draws, is.null, logical(1))]
  if (length(draws) == 0L) {
    ps_abort("No bootstrap replicate could be analysed; check the design of the trial.")
  }
  if (length(draws) < n_boot / 2) {
    ps_warn(paste0(
      "Only ", length(draws), " of ", n_boot, " bootstrap replicates were usable; ",
      "the interval is based on few resamples."
    ))
  }

  boots <- dplyr::bind_rows(draws)
  lo <- (1 - conf_level) / 2
  hi <- 1 - lo

  summary_boot <- dplyr::summarise(
    dplyr::group_by(boots, .data$unit, .data$group),
    conf_low = unname(stats::quantile(.data$isi, lo, na.rm = TRUE)),
    conf_high = unname(stats::quantile(.data$isi, hi, na.rm = TRUE)),
    rank_low = unname(stats::quantile(.data$rank, lo, na.rm = TRUE, type = 1)),
    rank_high = unname(stats::quantile(.data$rank, hi, na.rm = TRUE, type = 1)),
    p_best = mean(.data$rank == 1, na.rm = TRUE),
    n_ok = dplyr::n(),
    .groups = "drop"
  )

  out <- dplyr::left_join(
    tibble::as_tibble(observed)[c("unit", "group", "isi", "rank")],
    summary_boot,
    by = c("unit", "group")
  )
  out <- out[c(
    "unit", "group", "isi", "conf_low", "conf_high",
    "rank", "rank_low", "rank_high", "p_best", "n_ok"
  )]
  out <- dplyr::arrange(out, .data$group, .data$rank)

  structure(
    out,
    class = c("plantstress_isi_ci", class(tibble::tibble())),
    plantstress = list(
      n_boot = n_boot,
      n_ok = length(draws),
      conf_level = conf_level,
      weighting = if (is.numeric(weights)) "manual" else weights,
      aggregate = aggregate,
      rank_by = rank_by,
      scale = scale
    )
  )
}

#' @param x A `plantstress_isi_ci` object.
#' @param ... Ignored.
#' @rdname stress_index_ci
#' @export
print.plantstress_isi_ci <- function(x, ...) {
  meta <- sri_meta(x)
  cat("<plantstress_isi_ci>\n")
  cat("  Resamples: ", meta$n_ok %||% "?", " of ", meta$n_boot %||% "?",
    " usable\n",
    sep = ""
  )
  cat("  Interval:  ", format(100 * (meta$conf_level %||% NA)), "% percentile\n", sep = "")
  cat("  Ranking:   ", meta$rank_by %||% "?",
    " (rank 1 = ", if (identical(meta$rank_by, "severity")) "most stressed" else "most tolerant",
    ")\n",
    sep = ""
  )
  cat("  p_best = share of resamples in which the unit ranked first\n")
  print(tibble::as_tibble(x))
  invisible(x)
}

#' @param top_n Draw only the `top_n` best-ranked units of each stress level.
#' @rdname stress_index_ci
#' @export
plot.plantstress_isi_ci <- function(x, top_n = NULL, ...) {
  df <- tibble::as_tibble(x)
  if (!is.null(top_n)) {
    if (!is.numeric(top_n) || length(top_n) != 1L || is.na(top_n) || top_n < 1) {
      ps_abort("`top_n` must be a numeric scalar >= 1.")
    }
    df <- dplyr::filter(
      dplyr::group_by(df, .data$group),
      .data$rank <= as.integer(top_n)
    )
    df <- dplyr::ungroup(df)
  }

  df$label <- factor(df$unit, levels = unique(df$unit[order(-df$isi)]))
  meta <- sri_meta(x)

  ggplot2::ggplot(df, ggplot2::aes(x = .data$isi, y = .data$label)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = .data$conf_low, xmax = .data$conf_high),
      height = 0.22, colour = "grey35"
    ) +
    ggplot2::geom_point(ggplot2::aes(fill = .data$p_best),
      shape = 21, size = 3.4, colour = "grey25"
    ) +
    ggplot2::scale_fill_gradient(
      low = "#FFFFFF", high = "#2C7BB6", limits = c(0, 1),
      name = "P(rank 1)"
    ) +
    ggplot2::facet_wrap(~ .data$group, scales = "free_y") +
    ggplot2::labs(
      x = "Integrated Stress Index",
      y = NULL,
      title = "Integrated stress index with bootstrap limits",
      subtitle = paste0(
        format(100 * (meta$conf_level %||% NA)), "% percentile interval from ",
        meta$n_ok %||% "?", " resamples; higher index = more affected"
      )
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold")
    )
}
