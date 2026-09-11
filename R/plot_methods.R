#' Plot a Stress Signature
#'
#' @description
#' Shortcut for [plot_stress_signature()], so that `plot()` does the expected
#' thing on the object returned by [calculate_sri()].
#'
#' @param x A `plantstress_sri` object.
#' @param ... Passed to [plot_stress_signature()].
#'
#' @return A [ggplot2::ggplot()] object.
#'
#' @seealso [plot_stress_signature()]
#'
#' @examples
#' data(brachiaria_stress)
#' sri <- calculate_sri(brachiaria_stress,
#'   treatment = "drought_level", control = "control",
#'   traits = c("Fv_Fm", "PIabs", "A", "gs", "RWC"), verbose = FALSE
#' )
#' plot(sri)
#' @export
plot.plantstress_sri <- function(x, ...) {
  plot_stress_signature(x, ...)
}

#' Plot an Integrated Stress Index
#'
#' @description
#' Two views of the ranking produced by [integrated_stress_index()].
#'
#' `type = "ranking"` draws a lollipop chart of the index per unit, with the
#' rank printed inside each point and units ordered by their average index, so
#' that a faceted plot keeps one genotype on one row. It is the figure to show
#' when the question is *which genotype should be selected*.
#'
#' `type = "contribution"` decomposes each bar into the weighted contribution of
#' every trait. It is the figure to show when the question is *why* a genotype
#' ended up where it did, and it is the honest companion to any ranking built
#' from weighted traits.
#'
#' @param x A `plantstress_isi` object.
#' @param type `"ranking"` (default) or `"contribution"`.
#' @param top_n Optional. For `"contribution"`, keep only the `top_n` traits with
#'   the largest mean absolute contribution and pool the rest into `"other"`.
#' @param palette Character vector of length two with the low and high colors of
#'   the ranking gradient.
#' @param title,subtitle Optional plot annotations.
#' @param ... Ignored.
#'
#' @return A [ggplot2::ggplot()] object.
#'
#' @seealso [integrated_stress_index()], [stress_contributions()]
#'
#' @examples
#' data(brachiaria_stress)
#' sri <- calculate_sri(brachiaria_stress,
#'   treatment = "drought_level", control = "control",
#'   traits = c("Fv_Fm", "PIabs", "A", "gs", "RWC", "shoot_biomass"),
#'   by = "genotype", verbose = FALSE
#' )
#' isi <- integrated_stress_index(sri, weights = "precision")
#'
#' plot(isi)
#' plot(isi, type = "contribution")
#' @export
plot.plantstress_isi <- function(x,
                                 type = c("ranking", "contribution"),
                                 top_n = NULL,
                                 palette = c("#2C7BB6", "#D7191C"),
                                 title = NULL,
                                 subtitle = NULL,
                                 ...) {
  type <- match.arg(type)
  if (!inherits(x, "plantstress_isi")) {
    ps_abort("`x` must be a `plantstress_isi` object from `integrated_stress_index()`.")
  }
  if (!is.character(palette) || length(palette) != 2L) {
    ps_abort("`palette` must be a character vector of two colours (low, high).")
  }

  meta <- sri_meta(x)
  multi_group <- length(unique(x$group)) > 1L

  if (type == "ranking") {
    df <- tibble::as_tibble(x)
    # Units are ordered once, by their average index across stress levels, so
    # that a faceted plot keeps the same row for the same genotype.
    df$.unit <- stats::reorder(factor(df$unit), -df$isi, FUN = mean)

    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(x = .data$isi, y = .data$.unit, colour = .data$isi)
    ) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.4) +
      ggplot2::geom_segment(
        ggplot2::aes(x = 0, xend = .data$isi, y = .data$.unit, yend = .data$.unit),
        linewidth = 0.7
      ) +
      ggplot2::geom_point(size = 3.4) +
      ggplot2::geom_text(
        ggplot2::aes(label = .data$rank),
        colour = "white", size = 2.3, fontface = "bold"
      ) +
      ggplot2::scale_colour_gradient(low = palette[1], high = palette[2], name = "ISI") +
      ggplot2::labs(
        x = "Integrated stress index", y = NULL,
        title = title %||% "Integrated stress index",
        subtitle = subtitle %||% paste0(
          meta$weighting %||% "weighted", " weights; rank 1 = ",
          if (identical(meta$rank_by, "severity")) "most stressed" else "most tolerant"
        )
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
  } else {
    contrib <- stress_contributions(x)
    if (is.null(contrib) || nrow(contrib) == 0L) {
      ps_abort("No trait contributions are stored in `x`.")
    }
    df <- contrib

    if (!is.null(top_n)) {
      if (!is.numeric(top_n) || length(top_n) != 1L || is.na(top_n) || top_n < 1) {
        ps_abort("`top_n` must be a positive numeric scalar or `NULL`.")
      }
      strength <- vapply(
        split(abs(df$contribution), df$trait),
        function(v) mean(v, na.rm = TRUE), numeric(1)
      )
      keep <- names(sort(strength, decreasing = TRUE))[seq_len(min(top_n, length(strength)))]
      df$trait <- ifelse(df$trait %in% keep, df$trait, "other")
      df <- dplyr::summarise(
        dplyr::group_by(df, .data$unit, .data$group, .data$trait),
        contribution = sum(.data$contribution, na.rm = TRUE), .groups = "drop"
      )
    }

    order_by <- dplyr::summarise(
      dplyr::group_by(df, .data$unit),
      total = sum(.data$contribution, na.rm = TRUE), .groups = "drop"
    )
    df$.unit <- factor(df$unit, levels = order_by$unit[order(-order_by$total)])

    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(x = .data$contribution, y = .data$.unit, fill = .data$trait)
    ) +
      ggplot2::geom_col(width = 0.7, colour = "white", linewidth = 0.2) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.4) +
      ggplot2::scale_fill_viridis_d(name = "Trait", option = "D") +
      ggplot2::labs(
        x = "Weighted contribution to the index", y = NULL,
        title = title %||% "What drives the integrated stress index",
        subtitle = subtitle %||% "Bar length is the index; segments are the weighted traits"
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
  }

  if (multi_group) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$group))
  }
  p
}
