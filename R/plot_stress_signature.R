# Shared preparation for both signature displays: subset, rank and order the
# traits so that the two plot types always show the same panel of traits.
prepare_signature <- function(x, traits, groups, units, top_n, significant_only,
                              alpha, cluster) {
  if (!is.data.frame(x)) {
    ps_abort("`x` must be a `plantstress_sri` object or a data frame.")
  }
  required <- c("unit", "group", "trait", "sri")
  missing_cols <- setdiff(required, names(x))
  if (length(missing_cols) > 0L) {
    ps_abort(paste0(
      "`x` is missing required column(s): ", paste(missing_cols, collapse = ", "),
      ". Build it with `calculate_sri()`."
    ))
  }

  df <- tibble::as_tibble(x)
  if (!"p_adj" %in% names(df)) df$p_adj <- NA_real_

  if (!is.null(traits)) {
    unknown <- setdiff(traits, df$trait)
    if (length(unknown) > 0L) {
      ps_abort(paste0("Trait(s) absent from `x`: ", paste(unknown, collapse = ", "), "."))
    }
    df <- df[df$trait %in% traits, , drop = FALSE]
  }
  if (!is.null(groups)) {
    unknown <- setdiff(groups, df$group)
    if (length(unknown) > 0L) {
      ps_abort(paste0("Stress level(s) absent from `x`: ", paste(unknown, collapse = ", "), "."))
    }
    df <- df[df$group %in% groups, , drop = FALSE]
  }
  if (!is.null(units)) {
    unknown <- setdiff(units, df$unit)
    if (length(unknown) > 0L) {
      ps_abort(paste0("Unit(s) absent from `x`: ", paste(unknown, collapse = ", "), "."))
    }
    df <- df[df$unit %in% units, , drop = FALSE]
  }

  if (isTRUE(significant_only)) {
    check_prob(alpha, "alpha")
    keep_traits <- unique(df$trait[!is.na(df$p_adj) & df$p_adj <= alpha])
    if (length(keep_traits) == 0L) {
      ps_abort("No trait reached the significance threshold; relax `alpha` or set `significant_only = FALSE`.")
    }
    df <- df[df$trait %in% keep_traits, , drop = FALSE]
  }

  if (!is.null(top_n)) {
    if (!is.numeric(top_n) || length(top_n) != 1L || is.na(top_n) || top_n < 1) {
      ps_abort("`top_n` must be a positive numeric scalar or `NULL`.")
    }
    strength <- vapply(
      split(abs(df$sri), df$trait),
      function(v) max(v, na.rm = TRUE),
      numeric(1)
    )
    strength <- strength[is.finite(strength)]
    keep_traits <- names(sort(strength, decreasing = TRUE))[seq_len(min(top_n, length(strength)))]
    df <- df[df$trait %in% keep_traits, , drop = FALSE]
  }

  if (nrow(df) == 0L) {
    ps_abort("No rows left after filtering.")
  }

  trait_order <- if (isTRUE(cluster)) {
    cluster_order(sri_matrix(df))
  } else {
    mean_sri <- vapply(split(df$sri, df$trait), mean, numeric(1), na.rm = TRUE)
    names(sort(mean_sri, decreasing = TRUE))
  }
  trait_order <- trait_order[trait_order %in% df$trait]

  df$trait <- factor(df$trait, levels = trait_order)
  df$group <- factor(df$group, levels = unique(df$group))
  df$unit <- factor(df$unit, levels = unique(df$unit))
  df$label <- ifelse(
    is.na(df$p_adj), "",
    ifelse(df$p_adj <= 0.001, "***",
      ifelse(df$p_adj <= 0.01, "**",
        ifelse(df$p_adj <= 0.05, "*", "")
      )
    )
  )
  df
}

#' Plot the Physiological Stress Signature
#'
#' @description
#' Draws the multivariate stress signature held in a `plantstress_sri` object as
#' either a diverging heat map (traits by stress level) or a radar plot (one
#' closed polygon per stress level). Both displays share the same color
#' convention: positive values (warm) mean the trait was impaired by stress,
#' negative values (cool) mean it was enhanced, and zero means the stress
#' treatment was indistinguishable from the control.
#'
#' @details
#' The heat map is the default because it scales to dozens of traits and several
#' stress levels, and because trait rows can be ordered by hierarchical
#' clustering, which makes co-responding traits visually adjacent. The radar
#' plot is intended for a small, curated panel (typically 5 to 12 traits) in
#' presentations and manuscript figures; it requires at least three traits.
#'
#' When `x` contains several units (for example genotypes), the plot is
#' faceted by unit.
#'
#' @param x A `plantstress_sri` object from [calculate_sri()], or a data frame
#'   with the columns `unit`, `group`, `trait`, `sri` and, optionally, `p_adj`.
#' @param type `"heatmap"` (default) or `"radar"`.
#' @param traits Optional character vector restricting and not reordering the
#'   traits shown.
#' @param groups Optional character vector restricting the stress levels shown.
#' @param units Optional character vector restricting the units shown.
#' @param top_n Optional. Keep only the `top_n` traits with the largest absolute
#'   index.
#' @param significant_only Logical. Keep only traits significant at `alpha` in at
#'   least one cell.
#' @param alpha Adjusted p-value threshold used by `significant_only` and by the
#'   significance stars.
#' @param cluster Logical. Order traits by hierarchical clustering of their index
#'   profiles instead of by mean index.
#' @param show_significance Logical. Print significance stars on the heat map.
#' @param limits Optional numeric vector of length two giving the fill limits.
#'   By default the scale spans the observed range with the neutral color
#'   anchored at zero. Set it explicitly, for example `c(-8, 8)`, to make several
#'   figures directly comparable.
#' @param palette Character vector of length three with the low, mid and high
#'   colors.
#' @param title,subtitle Optional plot annotations.
#'
#' @return A [ggplot2::ggplot()] object, which can be further modified with the
#'   usual `+` syntax.
#'
#' @seealso [calculate_sri()], [integrated_stress_index()], [stress_network()]
#'
#' @examples
#' data(brachiaria_stress)
#'
#' sri <- calculate_sri(
#'   brachiaria_stress,
#'   treatment = "drought_level",
#'   control = "control",
#'   traits = c("Fv_Fm", "PIabs", "DI0_RC", "A", "gs", "RWC", "SPAD", "shoot_biomass"),
#'   verbose = FALSE
#' )
#'
#' plot_stress_signature(sri)
#' plot_stress_signature(sri, type = "radar", top_n = 6)
#' @export
plot_stress_signature <- function(x,
                                  type = c("heatmap", "radar"),
                                  traits = NULL,
                                  groups = NULL,
                                  units = NULL,
                                  top_n = NULL,
                                  significant_only = FALSE,
                                  alpha = 0.05,
                                  cluster = TRUE,
                                  show_significance = TRUE,
                                  limits = NULL,
                                  palette = c("#2C7BB6", "grey95", "#D7191C"),
                                  title = NULL,
                                  subtitle = NULL) {
  type <- match.arg(type)

  if (!is.character(palette) || length(palette) != 3L) {
    ps_abort("`palette` must be a character vector of three colours (low, mid, high).")
  }

  df <- prepare_signature(
    x = x, traits = traits, groups = groups, units = units, top_n = top_n,
    significant_only = significant_only, alpha = alpha, cluster = cluster
  )

  multi_unit <- nlevels(df$unit) > 1L
  if (!is.null(limits) && (!is.numeric(limits) || length(limits) != 2L || anyNA(limits))) {
    ps_abort("`limits` must be a numeric vector of length two or `NULL`.")
  }
  # With `limits = NULL`, ggplot2 anchors the neutral color at zero and
  # stretches each side of the palette over the observed range. Forcing
  # symmetry instead would spend half the palette on a sign the data may never
  # take, flattening the contrast between a mild and a severe response.

  if (type == "heatmap") {
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(x = .data$group, y = .data$trait, fill = .data$sri)
    ) +
      ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
      ggplot2::scale_fill_gradient2(
        low = palette[1], mid = palette[2], high = palette[3],
        midpoint = 0, limits = limits, name = "SRI"
      ) +
      ggplot2::labs(
        x = "Stress level", y = NULL,
        title = title %||% "Physiological stress signature",
        subtitle = subtitle %||% "Positive index = trait impaired relative to control"
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_text(face = "italic")
      )

    if (isTRUE(show_significance)) {
      p <- p + ggplot2::geom_text(
        ggplot2::aes(label = .data$label),
        size = 3, colour = "grey15", na.rm = TRUE
      )
    }
  } else {
    if (nlevels(df$trait) < 3L) {
      ps_abort("A radar plot needs at least three traits.")
    }
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(
        x = .data$trait, y = .data$sri,
        group = .data$group, colour = .data$group, fill = .data$group
      )
    ) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.4) +
      ggplot2::geom_polygon(alpha = 0.12, linewidth = 0.7, na.rm = TRUE) +
      ggplot2::geom_point(size = 1.8, na.rm = TRUE) +
      ggplot2::coord_polar() +
      ggplot2::labs(
        x = NULL, y = "Stress response index",
        colour = "Stress level", fill = "Stress level",
        title = title %||% "Physiological stress signature",
        subtitle = subtitle %||% "Distance from the zero circle = deviation from the control"
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(face = "italic"),
        panel.grid.minor = ggplot2::element_blank()
      )
  }

  if (multi_unit) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$unit))
  }
  p
}
