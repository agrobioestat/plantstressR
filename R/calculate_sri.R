# Standardized effect of a stress group relative to its own control group.
# Returns the *raw* (unoriented) effect; orientation by trait direction is
# applied by calculate_sri() so that the sign convention lives in one place.
sri_effect <- function(x_control, x_stress, method, conf_level) {
  x_control <- x_control[!is.na(x_control)]
  x_stress <- x_stress[!is.na(x_stress)]

  n_c <- length(x_control)
  n_s <- length(x_stress)
  empty <- list(
    n_control = n_c, n_stress = n_s,
    mean_control = NA_real_, sd_control = NA_real_,
    mean_stress = NA_real_, sd_stress = NA_real_,
    difference = NA_real_, effect = NA_real_, se = NA_real_,
    conf_low = NA_real_, conf_high = NA_real_,
    statistic = NA_real_, p_value = NA_real_
  )
  if (n_c < 2L || n_s < 2L) {
    return(empty)
  }

  m_c <- mean(x_control)
  m_s <- mean(x_stress)
  sd_c <- stats::sd(x_control)
  sd_s <- stats::sd(x_stress)
  diff <- m_s - m_c

  eff <- NA_real_
  se <- NA_real_

  if (method == "relative") {
    if (is.finite(m_c) && abs(m_c) > .Machine$double.eps^0.5) {
      eff <- diff / abs(m_c)
      se_c <- sd_c / sqrt(n_c)
      se_s <- sd_s / sqrt(n_s)
      # Delta method for the ratio of two independent means.
      se <- sqrt(se_s^2 / m_c^2 + (m_s^2 * se_c^2) / m_c^4)
    }
  } else {
    denom <- switch(method,
      glass = sd_c,
      cohen = sqrt(((n_c - 1) * sd_c^2 + (n_s - 1) * sd_s^2) / (n_c + n_s - 2)),
      hedges = sqrt(((n_c - 1) * sd_c^2 + (n_s - 1) * sd_s^2) / (n_c + n_s - 2))
    )
    if (is.finite(denom) && denom > 0) {
      d <- diff / denom
      df_var <- if (method == "glass") n_c - 1 else n_c + n_s - 2
      se_d <- sqrt((n_c + n_s) / (n_c * n_s) + d^2 / (2 * df_var))
      if (method == "hedges") {
        # Small-sample bias correction (Hedges' g).
        j <- 1 - 3 / (4 * (n_c + n_s) - 9)
        d <- j * d
        se_d <- j * se_d
      }
      eff <- d
      se <- se_d
    }
  }

  z <- z_multiplier(conf_level)
  tt <- try(stats::t.test(x_stress, x_control, var.equal = FALSE), silent = TRUE)
  statistic <- if (inherits(tt, "try-error")) NA_real_ else unname(tt$statistic)
  p_value <- if (inherits(tt, "try-error")) NA_real_ else tt$p.value

  list(
    n_control = n_c, n_stress = n_s,
    mean_control = m_c, sd_control = sd_c,
    mean_stress = m_s, sd_stress = sd_s,
    difference = diff, effect = eff, se = se,
    conf_low = eff - z * se, conf_high = eff + z * se,
    statistic = statistic, p_value = p_value
  )
}

#' Stress Response Index (SRI) per Trait
#'
#' @description
#' Converts each physiological trait into a Stress Response Index: the
#' standardized difference between a stress treatment and its own control,
#' re-signed so that **positive values always mean stress-induced impairment**,
#' whatever the direction of the trait.
#'
#' The magnitude is a standardized effect size, which puts traits measured in
#' incomparable units (`Fv/Fm`, micromoles of CO2 per square meter per second,
#' grams of biomass) on a common, dimensionless scale and makes them additive in
#' [integrated_stress_index()].
#'
#' @details
#' With \eqn{m_s} and \eqn{m_c} the stress and control means and \eqn{s} the
#' scaling standard deviation, the available `method` values are:
#'
#' \describe{
#'   \item{`"glass"`}{\eqn{(m_s - m_c) / s_c}, Glass's delta. The control group
#'     alone sets the scale, so a treatment that inflates variance does not
#'     shrink the index. This is the default.}
#'   \item{`"cohen"`}{\eqn{(m_s - m_c) / s_p}, with the pooled standard
#'     deviation.}
#'   \item{`"hedges"`}{Cohen's d with the small-sample bias correction, useful
#'     with few replicates per plot.}
#'   \item{`"relative"`}{\eqn{(m_s - m_c) / |m_c|}, the classical relative
#'     change used in the stress-tolerance selection literature; interpretable
#'     as a proportion but only defined for ratio-scale traits with a non-zero
#'     control mean.}
#' }
#'
#' The raw effect is then multiplied by \eqn{-1} for traits declared
#' `"higher_is_better"` and left unchanged for traits declared
#' `"lower_is_better"` (see [trait_directions()]). Confidence limits are normal
#' approximations built from the effect-size standard error; the reported
#' `p_value` comes from a Welch two-sample t-test on the raw trait and is
#' adjusted across all tests with `p_adjust`.
#'
#' @param data A data frame with one row per experimental unit, containing the
#'   treatment column, the trait columns and, optionally, a grouping column.
#' @param treatment Name of the treatment column.
#' @param control Level of `treatment` treated as the control. Every other level
#'   is treated as a stress level and compared against it.
#' @param traits Character vector of trait columns. Defaults to every numeric
#'   column that is not a design column.
#' @param by Optional grouping column (typically genotype). Indices are computed
#'   within each level, against that level's own control.
#' @param direction How each trait responds to stress. Either `"auto"` (the
#'   default, resolved by [trait_directions()]), a single string applied to all
#'   traits, or a named vector of `"higher_is_better"` / `"lower_is_better"`
#'   (equivalently `-1` / `1`).
#' @param method Effect-size definition: `"glass"`, `"cohen"`, `"hedges"` or
#'   `"relative"`.
#' @param conf_level Confidence level for the interval around the index.
#' @param p_adjust Multiple-testing adjustment passed to [stats::p.adjust()].
#' @param verbose Logical. Emit design warnings from [validate_stress_data()].
#'
#' @return An object of class `"plantstress_sri"`, a tibble with one row per
#'   unit, stress level and trait, and columns `unit`, `group`, `trait`,
#'   `direction`, `n_control`, `n_stress`, `mean_control`, `sd_control`,
#'   `mean_stress`, `sd_stress`, `difference`, `effect`, `sri`, `se`,
#'   `conf_low`, `conf_high`, `statistic`, `p_value` and `p_adj`. `unit` is
#'   `"overall"` when `by` is `NULL`.
#'
#' @references
#' Fischer R.A., Maurer R. (1978). Drought resistance in spring wheat cultivars.
#' I. Grain yield responses. \doi{10.1071/AR9780897}
#'
#' Hedges L.V. (1981). Distribution theory for Glass's estimator of effect size
#' and related estimators. \doi{10.3102/10769986006002107}
#'
#' @seealso [integrated_stress_index()], [plot_stress_signature()],
#'   [trait_directions()], [validate_stress_data()]
#'
#' @examples
#' data(brachiaria_stress)
#'
#' sri <- calculate_sri(
#'   brachiaria_stress,
#'   treatment = "drought_level",
#'   control = "control",
#'   traits = c("Fv_Fm", "PIabs", "A", "gs", "RWC", "SPAD", "shoot_biomass")
#' )
#' sri
#'
#' # Per genotype, with an explicit direction for a trait the dictionary
#' # cannot know about
#' sri_g <- calculate_sri(
#'   brachiaria_stress,
#'   treatment = "drought_level",
#'   control = "control",
#'   traits = c("Fv_Fm", "A", "RWC", "DI0_RC"),
#'   by = "genotype",
#'   direction = c(DI0_RC = "lower_is_better")
#' )
#' summary(sri_g)
#' @export
calculate_sri <- function(data,
                          treatment,
                          control,
                          traits = NULL,
                          by = NULL,
                          direction = "auto",
                          method = c("glass", "cohen", "hedges", "relative"),
                          conf_level = 0.95,
                          p_adjust = "BH",
                          verbose = TRUE) {
  method <- match.arg(method)
  check_prob(conf_level, "conf_level", lower = 0.5, upper = 0.9999)
  check_p_adjust(p_adjust)

  chk <- validate_stress_data(
    data = data,
    treatment = treatment,
    control = control,
    traits = traits,
    by = by,
    verbose = verbose
  )

  traits <- chk$traits
  control <- chk$control
  stress_levels <- chk$stress_levels
  direction <- resolve_direction(direction, traits)
  ori <- direction_sign(direction)

  trt <- as.character(data[[treatment]])
  unit <- if (is.null(by)) rep("overall", nrow(data)) else as.character(data[[by]])
  units <- unique(stats::na.omit(unit))

  rows <- list()
  for (u in units) {
    in_unit <- !is.na(unit) & unit == u
    ctrl_rows <- in_unit & !is.na(trt) & trt == control
    for (g in stress_levels) {
      stress_rows <- in_unit & !is.na(trt) & trt == g
      if (!any(stress_rows)) next
      for (tr in traits) {
        st <- sri_effect(
          x_control = data[[tr]][ctrl_rows],
          x_stress = data[[tr]][stress_rows],
          method = method,
          conf_level = conf_level
        )
        s <- ori[[tr]]
        lo <- st$conf_low
        hi <- st$conf_high
        rows[[length(rows) + 1L]] <- tibble::tibble(
          unit = u,
          group = g,
          trait = tr,
          direction = unname(direction[[tr]]),
          n_control = st$n_control,
          n_stress = st$n_stress,
          mean_control = st$mean_control,
          sd_control = st$sd_control,
          mean_stress = st$mean_stress,
          sd_stress = st$sd_stress,
          difference = st$difference,
          effect = st$effect,
          sri = s * st$effect,
          se = st$se,
          conf_low = if (s > 0) lo else -hi,
          conf_high = if (s > 0) hi else -lo,
          statistic = st$statistic,
          p_value = st$p_value
        )
      }
    }
  }

  if (length(rows) == 0L) {
    ps_abort("No control/stress pair could be compared. Check `treatment`, `control` and `by`.")
  }

  out <- dplyr::bind_rows(rows)
  out$p_adj <- stats::p.adjust(out$p_value, method = p_adjust)
  out <- dplyr::arrange(out, .data$unit, .data$group, dplyr::desc(abs(.data$sri)))

  if (all(is.na(out$sri))) {
    ps_warn("All stress response indices are NA; check replication and trait variance.")
  }

  structure(
    out,
    class = c("plantstress_sri", class(tibble::tibble())),
    plantstress = list(
      treatment = treatment,
      control = control,
      by = by,
      traits = traits,
      direction = direction,
      method = method,
      conf_level = conf_level,
      p_adjust = p_adjust
    )
  )
}

sri_meta <- function(x) attr(x, "plantstress", exact = TRUE) %||% list()

#' @param x A `plantstress_sri` object.
#' @param n Number of rows to display.
#' @param ... Ignored.
#' @rdname calculate_sri
#' @export
print.plantstress_sri <- function(x, n = 10, ...) {
  meta <- sri_meta(x)
  cat("<plantstress_sri>\n")
  cat("  Method:    ", meta$method %||% "?", " (positive SRI = stress damage)\n", sep = "")
  cat("  Control:   ", meta$control %||% "?", "\n", sep = "")
  cat("  Groups:    ", paste(unique(x$group), collapse = ", "), "\n", sep = "")
  cat("  Units:     ", if (is.null(meta$by)) "overall" else meta$by, "\n", sep = "")
  cat("  Traits:    ", length(unique(x$trait)), "\n", sep = "")
  print(tibble::as_tibble(x)[c("unit", "group", "trait", "sri", "conf_low", "conf_high", "p_adj")],
    n = n
  )
  invisible(x)
}

#' @param object A `plantstress_sri` object.
#' @param alpha Significance threshold applied to the adjusted p-values.
#' @rdname calculate_sri
#' @export
summary.plantstress_sri <- function(object, alpha = 0.05, ...) {
  check_prob(alpha, "alpha")
  df <- tibble::as_tibble(object)
  df$significant <- !is.na(df$p_adj) & df$p_adj <= alpha
  out <- dplyr::summarise(
    dplyr::group_by(df, .data$unit, .data$group),
    n_traits = dplyr::n(),
    n_significant = sum(.data$significant),
    mean_sri = mean(.data$sri, na.rm = TRUE),
    mean_abs_sri = mean(abs(.data$sri), na.rm = TRUE),
    most_impaired = .data$trait[which.max(.data$sri)][1],
    most_resilient = .data$trait[which.min(.data$sri)][1],
    .groups = "drop"
  )
  out
}
