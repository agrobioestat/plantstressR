# Internal helpers -------------------------------------------------------
# Kept in a single file so that no exported function depends on package
# internals that are hard to trace during review.

`%||%` <- function(x, y) if (is.null(x)) y else x

# Consistent, classed errors make the package easy to test with
# testthat::expect_error(class = "plantstressR_error").
ps_abort <- function(message, class = NULL, call = rlang::caller_env()) {
  rlang::abort(message, class = c(class, "plantstressR_error"), call = call)
}

ps_warn <- function(message, class = NULL) {
  rlang::warn(message, class = c(class, "plantstressR_warning"))
}

check_string <- function(x, arg) {
  if (!is.character(x) || length(x) != 1L || is.na(x)) {
    ps_abort(paste0("`", arg, "` must be a single non-missing character string."))
  }
  invisible(x)
}

check_data <- function(data) {
  if (!is.data.frame(data)) {
    ps_abort("`data` must be a data frame or tibble.")
  }
  if (nrow(data) == 0L) {
    ps_abort("`data` has zero rows.")
  }
  invisible(data)
}

check_column <- function(data, column, arg) {
  check_string(column, arg)
  if (!column %in% names(data)) {
    ps_abort(paste0("Column `", column, "` (argument `", arg, "`) was not found in `data`."))
  }
  invisible(column)
}

check_prob <- function(x, arg, lower = 0, upper = 1) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < lower || x > upper) {
    ps_abort(paste0("`", arg, "` must be a numeric scalar between ", lower, " and ", upper, "."))
  }
  invisible(x)
}

check_p_adjust <- function(method) {
  check_string(method, "p_adjust")
  if (!method %in% stats::p.adjust.methods) {
    ps_abort(paste0(
      "`p_adjust` must be one of: ",
      paste(stats::p.adjust.methods, collapse = ", "), "."
    ))
  }
  invisible(method)
}

# Resolve which columns are physiological traits. Columns used as design
# factors are always excluded, even when they happen to be numeric.
resolve_traits <- function(data, traits = NULL, exclude = character()) {
  exclude <- unique(stats::na.omit(exclude))
  if (is.null(traits)) {
    numeric_cols <- names(data)[vapply(data, is.numeric, logical(1))]
    traits <- setdiff(numeric_cols, exclude)
  } else {
    if (!is.character(traits)) {
      ps_abort("`traits` must be a character vector of column names or `NULL`.")
    }
    missing_cols <- setdiff(traits, names(data))
    if (length(missing_cols) > 0L) {
      ps_abort(paste0(
        "Trait column(s) not found in `data`: ",
        paste(missing_cols, collapse = ", "), "."
      ))
    }
    not_numeric <- traits[!vapply(data[traits], is.numeric, logical(1))]
    if (length(not_numeric) > 0L) {
      ps_abort(paste0(
        "Trait column(s) are not numeric: ",
        paste(not_numeric, collapse = ", "), "."
      ))
    }
    traits <- setdiff(traits, exclude)
  }
  if (length(traits) == 0L) {
    ps_abort("No numeric trait columns are available after excluding design columns.")
  }
  traits
}

# Long SRI table -> wide numeric matrix (rows = unit x group, cols = trait).
sri_matrix <- function(x, value = "sri") {
  row_key <- paste(x$unit, x$group, sep = " | ")
  rows <- unique(row_key)
  cols <- unique(x$trait)
  mat <- matrix(
    NA_real_,
    nrow = length(rows), ncol = length(cols),
    dimnames = list(rows, cols)
  )
  idx <- cbind(match(row_key, rows), match(x$trait, cols))
  mat[idx] <- x[[value]]
  mat
}

# Mean imputation keeps prcomp() usable when a genotype x trait cell is
# missing; it is never used for inference, only for weighting and ordering.
impute_column_means <- function(mat) {
  for (j in seq_len(ncol(mat))) {
    na_j <- is.na(mat[, j])
    if (any(na_j)) {
      mu <- mean(mat[, j], na.rm = TRUE)
      mat[na_j, j] <- if (is.finite(mu)) mu else 0
    }
  }
  mat
}

# Hierarchical ordering of matrix columns; falls back to the input order when
# clustering is impossible (fewer than three columns or degenerate distances).
cluster_order <- function(mat) {
  if (is.null(mat) || ncol(mat) < 3L || nrow(mat) < 2L) {
    return(colnames(mat))
  }
  d <- try(stats::dist(t(impute_column_means(mat))), silent = TRUE)
  if (inherits(d, "try-error") || anyNA(d) || all(d == 0)) {
    return(colnames(mat))
  }
  hc <- try(stats::hclust(d, method = "average"), silent = TRUE)
  if (inherits(hc, "try-error")) {
    return(colnames(mat))
  }
  colnames(mat)[hc$order]
}

# Format a two-sided normal confidence multiplier.
z_multiplier <- function(conf_level) {
  stats::qnorm(1 - (1 - conf_level) / 2)
}

# Block (randomized complete block) support ------------------------------
# Blocks are nuisance strata: they shift a whole replicate up or down without
# changing what the treatment did. Every block-aware routine in the package
# removes the same additive, mean-centred block effect estimated from
# `value ~ treatment + block`, so that treatment means are preserved exactly
# and only the between-block variance is taken out.

# Centred block effects, or NULL when the design cannot support the fit.
block_effects <- function(value, trt, blk) {
  ok <- !is.na(value) & !is.na(trt) & !is.na(blk)
  if (sum(ok) < 3L) {
    return(NULL)
  }
  v <- value[ok]
  f_trt <- factor(as.character(trt[ok]))
  f_blk <- factor(as.character(blk[ok]))
  if (nlevels(f_blk) < 2L) {
    return(NULL)
  }

  fit <- if (nlevels(f_trt) < 2L) {
    try(stats::lm(v ~ f_blk), silent = TRUE)
  } else {
    try(stats::lm(v ~ f_trt + f_blk), silent = TRUE)
  }
  if (inherits(fit, "try-error")) {
    return(NULL)
  }
  df_res <- stats::df.residual(fit)
  if (is.na(df_res) || df_res < 1L) {
    return(NULL)
  }

  cf <- stats::coef(fit)
  idx <- grep("^f_blk", names(cf))
  eff <- stats::setNames(rep(0, nlevels(f_blk)), levels(f_blk))
  if (length(idx) > 0L) {
    eff[sub("^f_blk", "", names(cf)[idx])] <- cf[idx]
  }
  eff[!is.finite(eff)] <- 0
  eff - mean(eff)
}

# Subtract the centred block effect from each observation. Values whose block
# could not be estimated are returned untouched.
block_adjust <- function(value, trt, blk) {
  eff <- block_effects(value, trt, blk)
  if (is.null(eff)) {
    return(value)
  }
  shift <- eff[as.character(blk)]
  shift[is.na(shift)] <- 0
  value - unname(shift)
}

# Apply the adjustment trait by trait over a whole table.
block_adjust_traits <- function(data, traits, treatment, block) {
  trt <- if (is.null(treatment)) rep("all", nrow(data)) else as.character(data[[treatment]])
  blk <- as.character(data[[block]])
  for (tr in traits) {
    data[[tr]] <- block_adjust(data[[tr]], trt, blk)
  }
  data
}

# Residual standard deviation of one treatment group once its block effects
# have been accounted for; falls back to the plain sd when blocks are unusable.
block_residual_sd <- function(value, blk) {
  ok <- !is.na(value) & !is.na(blk)
  v <- value[ok]
  b <- factor(as.character(blk[ok]))
  if (length(v) < 2L) {
    return(NA_real_)
  }
  if (nlevels(b) < 2L || nlevels(b) >= length(v)) {
    return(stats::sd(v))
  }
  fit <- try(stats::lm(v ~ b), silent = TRUE)
  if (inherits(fit, "try-error") || stats::df.residual(fit) < 1L) {
    return(stats::sd(v))
  }
  stats::sigma(fit)
}

# Temporarily set the RNG seed and restore the caller's stream on exit, so that
# reproducible graph layouts never leak into the user's random number state.
ps_local_seed <- function(seed, env = parent.frame()) {
  if (is.null(seed)) {
    return(invisible(NULL))
  }
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) {
    ps_abort("`seed` must be a numeric scalar or `NULL`.")
  }
  has_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if (has_seed) get(".Random.seed", envir = globalenv()) else NULL
  do.call(
    base::on.exit,
    list(
      substitute(
        if (is.null(old_seed)) {
          suppressWarnings(rm(".Random.seed", envir = globalenv()))
        } else {
          assign(".Random.seed", old_seed, envir = globalenv())
        },
        list(old_seed = old_seed)
      ),
      add = TRUE
    ),
    envir = env
  )
  set.seed(as.integer(seed))
  invisible(NULL)
}
