#' Validate a Control Versus Stress Trait Table
#'
#' @description
#' Runs the design checks that every other `plantstressR` function relies on,
#' before any index is computed: presence and type of the design columns,
#' existence of the control level, replication of each design cell, constant or
#' fully missing traits, and the amount of missing data per trait.
#'
#' Hard problems (missing columns, an unknown control level, no numeric traits)
#' raise an error. Soft problems (low replication, constant traits, high
#' missingness) are collected and returned so that the caller can decide what to
#' do, and are emitted as warnings when `verbose = TRUE`.
#'
#' @param data A data frame with one row per experimental unit.
#' @param treatment Name of the column holding the treatment factor.
#' @param control Level of `treatment` used as the control (reference) group.
#' @param traits Character vector of trait columns. Defaults to every numeric
#'   column that is not a design column.
#' @param by Optional name of a grouping column (typically genotype, cultivar or
#'   site). When supplied, every level must contain both control and stress
#'   observations.
#' @param block Optional name of a block (replicate) column. When supplied, the
#'   design is additionally checked for enough block levels, for treatments
#'   confined to a single block (which cannot be separated from it) and for
#'   empty treatment-by-block cells.
#' @param min_replicates Minimum number of observations expected per design
#'   cell. Cells below this threshold are flagged.
#' @param max_missing Maximum acceptable proportion of missing values per trait.
#' @param verbose Logical. Emit the collected issues as warnings.
#'
#' @return Invisibly, an object of class `"plantstress_validation"`: a list with
#'   `ok`, `issues` (a tibble with columns `severity`, `check` and `message`),
#'   `design` (replication per design cell), `traits`, `treatment`, `control`,
#'   `stress_levels`, `by` and `block`.
#'
#' @seealso [calculate_sri()]
#'
#' @examples
#' data(brachiaria_stress)
#' chk <- validate_stress_data(
#'   brachiaria_stress,
#'   treatment = "drought_level",
#'   control = "control",
#'   traits = c("Fv_Fm", "A", "RWC", "shoot_biomass"),
#'   by = "genotype"
#' )
#' chk
#' @export
validate_stress_data <- function(data,
                                 treatment,
                                 control,
                                 traits = NULL,
                                 by = NULL,
                                 block = NULL,
                                 min_replicates = 3L,
                                 max_missing = 0.2,
                                 verbose = TRUE) {
  check_data(data)
  check_column(data, treatment, "treatment")
  if (!is.null(by)) check_column(data, by, "by")
  if (!is.null(block)) check_column(data, block, "block")
  check_prob(max_missing, "max_missing")
  if (!is.numeric(min_replicates) || length(min_replicates) != 1L ||
    is.na(min_replicates) || min_replicates < 2) {
    ps_abort("`min_replicates` must be a numeric scalar >= 2.")
  }

  trt <- as.character(data[[treatment]])
  if (length(control) != 1L || is.na(control)) {
    ps_abort("`control` must be a single non-missing treatment level.")
  }
  control <- as.character(control)
  levels_present <- sort(unique(stats::na.omit(trt)))
  if (!control %in% levels_present) {
    ps_abort(paste0(
      "Control level `", control, "` was not found in column `", treatment,
      "`. Levels present: ", paste(levels_present, collapse = ", "), "."
    ))
  }
  stress_levels <- setdiff(levels_present, control)
  if (length(stress_levels) == 0L) {
    ps_abort(paste0(
      "Column `", treatment, "` contains only the control level; ",
      "at least one stress level is required."
    ))
  }

  traits <- resolve_traits(data, traits, exclude = c(treatment, by, block))

  issues <- list()
  add_issue <- function(severity, check, message) {
    issues[[length(issues) + 1L]] <<- tibble::tibble(
      severity = severity, check = check, message = message
    )
  }

  if (anyNA(trt)) {
    add_issue(
      "warning", "treatment_na",
      paste0(
        sum(is.na(trt)), " row(s) have a missing `", treatment,
        "` value and are ignored."
      )
    )
  }

  unit <- if (is.null(by)) rep("overall", nrow(data)) else as.character(data[[by]])
  design <- tibble::as_tibble(
    as.data.frame(table(unit = unit, group = trt), stringsAsFactors = FALSE)
  )
  names(design)[names(design) == "Freq"] <- "n"
  design <- dplyr::arrange(design, .data$unit, .data$group)

  low <- design[design$n < min_replicates, , drop = FALSE]
  if (nrow(low) > 0L) {
    add_issue(
      "warning", "replication",
      paste0(
        nrow(low), " design cell(s) have fewer than ", min_replicates,
        " observations: ",
        paste(paste0(low$unit, "/", low$group, " (n=", low$n, ")"), collapse = "; "),
        "."
      )
    )
  }

  if (!is.null(by)) {
    without_control <- setdiff(unique(unit), unique(unit[trt %in% control]))
    if (length(without_control) > 0L) {
      add_issue(
        "error", "missing_control",
        paste0(
          "Level(s) of `", by, "` without control observations: ",
          paste(without_control, collapse = ", "), "."
        )
      )
    }
  }

  if (!is.null(block)) {
    blk <- as.character(data[[block]])
    if (anyNA(blk)) {
      add_issue(
        "warning", "block_na",
        paste0(
          sum(is.na(blk)), " row(s) have a missing `", block,
          "` value and are ignored by the block adjustment."
        )
      )
    }
    n_blocks <- length(unique(stats::na.omit(blk)))
    if (n_blocks < 2L) {
      add_issue(
        "error", "block_levels",
        paste0(
          "Column `", block, "` has ", n_blocks,
          " level(s); at least two are required to estimate block effects."
        )
      )
    } else {
      # A treatment confined to a single block is not separable from it.
      per_trt <- tapply(blk, trt, function(b) length(unique(stats::na.omit(b))))
      confined <- names(per_trt)[!is.na(per_trt) & per_trt < 2L]
      if (length(confined) > 0L) {
        add_issue(
          "error", "block_confounded",
          paste0(
            "Treatment level(s) present in a single block, so treatment and ",
            "block cannot be separated: ", paste(confined, collapse = ", "), "."
          )
        )
      }
      cells <- table(trt, blk)
      if (any(cells == 0L)) {
        add_issue(
          "warning", "block_incomplete",
          paste0(
            sum(cells == 0L), " treatment x block cell(s) are empty: the design ",
            "is not a complete block layout and the adjusted means rely on the ",
            "additive model to fill the gaps."
          )
        )
      }
    }
  }

  miss_rate <- vapply(data[traits], function(x) mean(is.na(x)), numeric(1))
  all_na <- names(miss_rate)[miss_rate == 1]
  if (length(all_na) > 0L) {
    ps_abort(paste0(
      "Trait column(s) are entirely missing: ", paste(all_na, collapse = ", "), "."
    ))
  }
  high_miss <- names(miss_rate)[miss_rate > max_missing]
  if (length(high_miss) > 0L) {
    add_issue(
      "warning", "missing_values",
      paste0(
        "Trait(s) above ", round(100 * max_missing), "% missing values: ",
        paste(sprintf("%s (%.0f%%)", high_miss, 100 * miss_rate[high_miss]),
          collapse = ", "
        ), "."
      )
    )
  }

  control_rows <- which(trt %in% control)
  constant <- vapply(
    data[traits],
    function(x) {
      v <- stats::sd(x[control_rows], na.rm = TRUE)
      !is.finite(v) || v == 0
    },
    logical(1)
  )
  if (any(constant)) {
    add_issue(
      "warning", "constant_trait",
      paste0(
        "Trait(s) with zero or undefined variance in the control group ",
        "(standardized indices will be NA): ",
        paste(traits[constant], collapse = ", "), "."
      )
    )
  }

  issues_tbl <- if (length(issues) == 0L) {
    tibble::tibble(severity = character(), check = character(), message = character())
  } else {
    dplyr::bind_rows(issues)
  }

  hard <- issues_tbl[issues_tbl$severity == "error", , drop = FALSE]
  if (nrow(hard) > 0L) {
    ps_abort(paste(hard$message, collapse = " "))
  }

  if (isTRUE(verbose) && nrow(issues_tbl) > 0L) {
    for (msg in issues_tbl$message) ps_warn(msg)
  }

  out <- list(
    ok = nrow(issues_tbl) == 0L,
    issues = issues_tbl,
    design = design,
    traits = traits,
    treatment = treatment,
    control = control,
    stress_levels = stress_levels,
    by = by,
    block = block,
    n_rows = nrow(data)
  )
  class(out) <- "plantstress_validation"
  invisible(out)
}

#' @param x A `plantstress_validation` object.
#' @param ... Ignored.
#' @rdname validate_stress_data
#' @export
print.plantstress_validation <- function(x, ...) {
  cat("<plantstress_validation>\n")
  cat("  Rows:          ", x$n_rows, "\n", sep = "")
  cat("  Treatment:     ", x$treatment, " (control = ", x$control, ")\n", sep = "")
  cat("  Stress levels: ", paste(x$stress_levels, collapse = ", "), "\n", sep = "")
  cat("  Grouping:      ", x$by %||% "<none>", "\n", sep = "")
  cat("  Block:         ", x$block %||% "<none>", "\n", sep = "")
  cat("  Traits:        ", length(x$traits), "\n", sep = "")
  if (nrow(x$issues) == 0L) {
    cat("  Status:        no issues detected\n")
  } else {
    cat("  Status:        ", nrow(x$issues), " issue(s)\n", sep = "")
    for (i in seq_len(nrow(x$issues))) {
      cat("   - [", x$issues$check[i], "] ", x$issues$message[i], "\n", sep = "")
    }
  }
  invisible(x)
}
