# Definition and selection direction of every classical tolerance index.
# `higher` = TRUE means a larger value identifies a more tolerant unit.
.ps_sti_registry <- list(
  TOL = list(higher = FALSE, label = "Tolerance"),
  MP = list(higher = TRUE, label = "Mean productivity"),
  GMP = list(higher = TRUE, label = "Geometric mean productivity"),
  HM = list(higher = TRUE, label = "Harmonic mean"),
  SSI = list(higher = FALSE, label = "Stress susceptibility index"),
  STI = list(higher = TRUE, label = "Stress tolerance index"),
  YI = list(higher = TRUE, label = "Yield index"),
  YSI = list(higher = TRUE, label = "Yield stability index"),
  RDI = list(higher = TRUE, label = "Relative drought index"),
  SSPI = list(higher = FALSE, label = "Stress susceptibility percentage index")
)

compute_tolerance_index <- function(name, yp, ys, mean_yp, mean_ys) {
  switch(name,
    TOL = yp - ys,
    MP = (yp + ys) / 2,
    GMP = sqrt(yp * ys),
    HM = 2 * yp * ys / (yp + ys),
    SSI = (1 - ys / yp) / (1 - mean_ys / mean_yp),
    STI = (yp * ys) / mean_yp^2,
    YI = ys / mean_ys,
    YSI = ys / yp,
    RDI = (ys / yp) / (mean_ys / mean_yp),
    SSPI = 100 * (yp - ys) / (2 * mean_yp)
  )
}

#' Classical Stress Tolerance and Susceptibility Indices
#'
#' @description
#' Computes the productivity-based selection indices that plant breeders and
#' stress physiologists have used for decades, from the performance of each
#' genotype under control (\eqn{Y_p}) and under stress (\eqn{Y_s}) conditions.
#'
#' These indices answer a different question from [calculate_sri()]. The Stress
#' Response Index describes the *physiological signature* across a whole trait
#' panel; the tolerance indices summarize a *single productivity trait* and are
#' the currency in which selection decisions are usually written. Both views are
#' offered so that a physiological ranking can be checked against the agronomic
#' one.
#'
#' @details
#' Let \eqn{\bar{Y}_p} and \eqn{\bar{Y}_s} be the trial means across units.
#' The available indices are:
#'
#' \describe{
#'   \item{`TOL`}{\eqn{Y_p - Y_s}. Yield loss. Lower is better.}
#'   \item{`MP`}{\eqn{(Y_p + Y_s)/2}. Mean productivity. Higher is better.}
#'   \item{`GMP`}{\eqn{\sqrt{Y_p Y_s}}. Geometric mean productivity, less biased
#'     towards the high-potential end than `MP`. Higher is better.}
#'   \item{`HM`}{\eqn{2 Y_p Y_s / (Y_p + Y_s)}. Harmonic mean. Higher is better.}
#'   \item{`SSI`}{\eqn{(1 - Y_s/Y_p) / (1 - \bar{Y}_s/\bar{Y}_p)}. Stress
#'     susceptibility, scaled by the stress intensity of the trial. Lower is
#'     better.}
#'   \item{`STI`}{\eqn{Y_p Y_s / \bar{Y}_p^2}. Identifies units that perform
#'     well in both environments. Higher is better.}
#'   \item{`YI`}{\eqn{Y_s / \bar{Y}_s}. Yield index. Higher is better.}
#'   \item{`YSI`}{\eqn{Y_s / Y_p}. Yield stability. Higher is better.}
#'   \item{`RDI`}{\eqn{(Y_s/Y_p) / (\bar{Y}_s/\bar{Y}_p)}. Relative drought
#'     index. Higher is better.}
#'   \item{`SSPI`}{\eqn{100 (Y_p - Y_s) / (2 \bar{Y}_p)}. Yield loss as a
#'     percentage of the trial potential. Lower is better.}
#' }
#'
#' No single index is correct on its own: `TOL` and `SSI` reward stability even
#' when it comes with low potential, while `MP` rewards potential even when the
#' loss is large. `rank_mean`, the average rank across the requested indices
#' after orienting each one towards tolerance, is therefore reported as well, and
#' `rank_overall` ranks the units on it.
#'
#' The stress intensity of each stress level, \eqn{1 - \bar{Y}_s/\bar{Y}_p}, is
#' returned in the `stress_intensity` column; indices that divide by it are `NA`
#' when a stress level did not reduce the trial mean at all.
#'
#' @param data A data frame with one row per experimental unit.
#' @param trait Name of the single productivity trait (yield, biomass) to
#'   summarize.
#' @param treatment Name of the treatment column.
#' @param control Level of `treatment` used as the non-stress environment.
#' @param by Name of the column identifying the units being compared, typically
#'   genotype. Required: every index is relative to the trial means.
#' @param indices Character vector of indices to compute. Defaults to all ten.
#' @param fun Function used to summarize `trait` within each unit and treatment.
#'
#' @return An object of class `"plantstress_sti"`, a tibble with columns `unit`,
#'   `group`, `Yp`, `Ys`, `stress_intensity`, one column per requested index,
#'   `rank_mean` and `rank_overall`.
#'
#' @references
#' Fischer R.A., Maurer R. (1978). Drought resistance in spring wheat cultivars.
#' I. Grain yield responses. \doi{10.1071/AR9780897}
#'
#' Rosielle A.A., Hamblin J. (1981). Theoretical aspects of selection for yield
#' in stress and non-stress environments.
#' \doi{10.2135/cropsci1981.0011183X002100060033x}
#'
#' Fernandez G.C.J. (1992). Effective selection criteria for assessing plant
#' stress tolerance. In: Adaptation of Food Crops to Temperature and Water Stress.
#' AVRDC, Taiwan, 257-270.
#'
#' @seealso [calculate_sri()], [integrated_stress_index()]
#'
#' @examples
#' data(brachiaria_stress)
#'
#' sti <- stress_tolerance_index(
#'   brachiaria_stress,
#'   trait = "shoot_biomass",
#'   treatment = "drought_level",
#'   control = "control",
#'   by = "genotype"
#' )
#' sti
#'
#' # A compact, defensible selection set
#' stress_tolerance_index(
#'   brachiaria_stress,
#'   trait = "shoot_biomass",
#'   treatment = "drought_level",
#'   control = "control",
#'   by = "genotype",
#'   indices = c("STI", "GMP", "SSI")
#' )
#' @export
stress_tolerance_index <- function(data,
                                   trait,
                                   treatment,
                                   control,
                                   by,
                                   indices = names(.ps_sti_registry),
                                   fun = mean) {
  check_data(data)
  check_column(data, treatment, "treatment")
  check_column(data, by, "by")
  check_column(data, trait, "trait")
  if (!is.numeric(data[[trait]])) {
    ps_abort(paste0("Trait column `", trait, "` must be numeric."))
  }
  if (!is.function(fun)) {
    ps_abort("`fun` must be a function, for example `mean` or `median`.")
  }
  if (!is.character(indices) || length(indices) == 0L) {
    ps_abort("`indices` must be a non-empty character vector.")
  }
  unknown <- setdiff(indices, names(.ps_sti_registry))
  if (length(unknown) > 0L) {
    ps_abort(paste0(
      "Unknown index/indices: ", paste(unknown, collapse = ", "),
      ". Available: ", paste(names(.ps_sti_registry), collapse = ", "), "."
    ))
  }
  indices <- unique(indices)

  trt <- as.character(data[[treatment]])
  control <- as.character(control)
  levels_present <- sort(unique(stats::na.omit(trt)))
  if (!control %in% levels_present) {
    ps_abort(paste0(
      "Control level `", control, "` was not found in column `", treatment, "`."
    ))
  }
  stress_levels <- setdiff(levels_present, control)
  if (length(stress_levels) == 0L) {
    ps_abort("At least one stress level is required.")
  }

  unit <- as.character(data[[by]])
  units <- sort(unique(stats::na.omit(unit)))
  if (length(units) < 2L) {
    ps_abort(paste0(
      "`by` must identify at least two units; the indices are relative to the ",
      "trial means and are undefined for a single unit."
    ))
  }

  summarise_cell <- function(u, g) {
    v <- data[[trait]][!is.na(unit) & unit == u & !is.na(trt) & trt == g]
    v <- v[!is.na(v)]
    if (length(v) == 0L) {
      return(NA_real_)
    }
    out <- fun(v)
    if (!is.numeric(out) || length(out) != 1L) {
      ps_abort("`fun` must return a single numeric value.")
    }
    as.numeric(out)
  }

  yp <- vapply(units, summarise_cell, numeric(1), g = control)
  missing_control <- units[is.na(yp)]
  if (length(missing_control) > 0L) {
    ps_abort(paste0(
      "Unit(s) without usable control observations: ",
      paste(missing_control, collapse = ", "), "."
    ))
  }
  if (any(yp <= 0)) {
    ps_warn(paste0(
      "Non-positive control performance for unit(s): ",
      paste(units[yp <= 0], collapse = ", "),
      ". Ratio-based indices will be NA or unstable."
    ))
  }

  rows <- list()
  for (g in stress_levels) {
    ys <- vapply(units, summarise_cell, numeric(1), g = g)
    keep <- !is.na(ys)
    if (!any(keep)) next

    mean_yp <- mean(yp[keep])
    mean_ys <- mean(ys[keep])
    si <- 1 - mean_ys / mean_yp

    # `vapply()` names its result after `units`; strip the names so that the
    # returned columns are plain numeric vectors.
    tbl <- tibble::tibble(
      unit = unname(units[keep]),
      group = g,
      Yp = unname(yp[keep]),
      Ys = unname(ys[keep]),
      stress_intensity = si
    )

    for (nm in indices) {
      value <- unname(compute_tolerance_index(nm, tbl$Yp, tbl$Ys, mean_yp, mean_ys))
      value[!is.finite(value)] <- NA_real_
      tbl[[nm]] <- value
    }

    # Orient every index towards tolerance before averaging the ranks, so that
    # a low TOL and a high STI both push a genotype to the top of the list.
    rank_mat <- vapply(
      indices,
      function(nm) {
        v <- tbl[[nm]]
        if (all(is.na(v))) {
          return(rep(NA_real_, nrow(tbl)))
        }
        sign_nm <- if (isTRUE(.ps_sti_registry[[nm]]$higher)) -1 else 1
        rank(sign_nm * v, ties.method = "min", na.last = "keep")
      },
      numeric(nrow(tbl))
    )
    rank_mat <- matrix(rank_mat, nrow = nrow(tbl))

    tbl$rank_mean <- rowMeans(rank_mat, na.rm = TRUE)
    tbl$rank_mean[is.nan(tbl$rank_mean)] <- NA_real_
    tbl$rank_overall <- rank(tbl$rank_mean, ties.method = "min", na.last = "keep")

    rows[[length(rows) + 1L]] <- dplyr::arrange(tbl, .data$rank_overall)
  }

  if (length(rows) == 0L) {
    ps_abort("No stress level had usable observations.")
  }

  out <- dplyr::bind_rows(rows)

  structure(
    out,
    class = c("plantstress_sti", class(tibble::tibble())),
    plantstress = list(
      trait = trait,
      treatment = treatment,
      control = control,
      by = by,
      indices = indices
    )
  )
}

#' @param x A `plantstress_sti` object.
#' @param ... Ignored.
#' @rdname stress_tolerance_index
#' @export
print.plantstress_sti <- function(x, ...) {
  meta <- sri_meta(x)
  cat("<plantstress_sti>\n")
  cat("  Trait:   ", meta$trait %||% "?", "\n", sep = "")
  cat("  Units:   ", length(unique(x$unit)), " (", meta$by %||% "?", ")\n", sep = "")
  cat("  Indices: ", paste(meta$indices, collapse = ", "), "\n", sep = "")
  si <- unique(x[c("group", "stress_intensity")])
  for (i in seq_len(nrow(si))) {
    cat("  Stress intensity [", si$group[i], "]: ",
      sprintf("%.3f", si$stress_intensity[i]), "\n",
      sep = ""
    )
  }
  cat("  rank_overall 1 = most tolerant on the average of the indices\n")
  print(tibble::as_tibble(x))
  invisible(x)
}
