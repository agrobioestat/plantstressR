# Traits whose *increase* indicates damage rather than performance. Names are
# matched case-insensitively after removing separators, so "DIo/RC", "DI0_RC"
# and "di0rc" all resolve to the same entry.
.ps_lower_is_better <- c(
  "fo", "f0",
  "absrc", "di0rc", "diorc",
  "ci",
  "mda", "h2o2", "electrolyteleakage", "el",
  "canopytemperature", "cwsi",
  "na", "cl"
)

normalize_trait_name <- function(x) {
  x <- tolower(x)
  gsub("[^a-z0-9]", "", x)
}

#' Declare the Stress Direction of Each Trait
#'
#' @description
#' The Stress Response Index produced by [calculate_sri()] is *signed*: positive
#' values must always mean "the stress treatment impaired this trait". Whether a
#' numeric increase is good or bad depends on the trait, so every trait needs a
#' direction. `trait_directions()` builds and validates that named vector.
#'
#' Traits declared as `"higher_is_better"` (net photosynthesis, `Fv/Fm`,
#' relative water content, biomass) have their raw effect multiplied by `-1`;
#' traits declared as `"lower_is_better"` (`F0`, `DI0/RC`, `Ci`,
#' malondialdehyde, electrolyte leakage, canopy temperature) keep the raw sign.
#'
#' @param traits Character vector of trait names.
#' @param lower_is_better Character vector of trait names whose increase denotes
#'   stress damage. Overrides the built-in dictionary.
#' @param higher_is_better Character vector of trait names whose increase denotes
#'   better performance. Overrides the built-in dictionary.
#' @param default Direction used for traits that are not matched by the
#'   dictionary nor by the explicit arguments. One of `"higher_is_better"` or
#'   `"lower_is_better"`.
#' @param use_dictionary Logical. When `TRUE` (default) a small dictionary of
#'   well-established plant-physiology traits is consulted before `default` is
#'   applied. Set to `FALSE` for full manual control.
#'
#' @return A named character vector, one element per trait, with values
#'   `"higher_is_better"` or `"lower_is_better"`.
#'
#' @references
#' Strasser R.J., Tsimilli-Michael M., Srivastava A. (2004). Analysis of the
#' chlorophyll a fluorescence transient. \doi{10.1007/978-1-4020-3218-9_12}
#'
#' @seealso [calculate_sri()]
#'
#' @examples
#' trait_directions(c("A", "Fv_Fm", "DI0_RC", "Ci", "shoot_biomass"))
#'
#' # Override the dictionary for a trait measured in a non-standard way
#' trait_directions(c("A", "WUE"), lower_is_better = "WUE")
#' @export
trait_directions <- function(traits,
                             lower_is_better = NULL,
                             higher_is_better = NULL,
                             default = c("higher_is_better", "lower_is_better"),
                             use_dictionary = TRUE) {
  default <- match.arg(default)

  if (!is.character(traits) || length(traits) == 0L || anyNA(traits)) {
    ps_abort("`traits` must be a non-empty character vector without missing values.")
  }
  if (!is.logical(use_dictionary) || length(use_dictionary) != 1L || is.na(use_dictionary)) {
    ps_abort("`use_dictionary` must be `TRUE` or `FALSE`.")
  }

  overlap <- intersect(lower_is_better, higher_is_better)
  if (length(overlap) > 0L) {
    ps_abort(paste0(
      "Trait(s) declared in both `lower_is_better` and `higher_is_better`: ",
      paste(overlap, collapse = ", "), "."
    ))
  }

  out <- rep(default, length(traits))
  names(out) <- traits

  if (use_dictionary) {
    is_dict <- normalize_trait_name(traits) %in% normalize_trait_name(.ps_lower_is_better)
    out[is_dict] <- "lower_is_better"
  }

  unknown <- setdiff(c(lower_is_better, higher_is_better), traits)
  if (length(unknown) > 0L) {
    ps_warn(paste0(
      "Ignoring direction override for trait(s) absent from `traits`: ",
      paste(unknown, collapse = ", "), "."
    ))
  }

  out[intersect(traits, lower_is_better)] <- "lower_is_better"
  out[intersect(traits, higher_is_better)] <- "higher_is_better"
  out
}

# Accepts what users may realistically pass to `calculate_sri(direction = )`:
# "auto", a named character vector, or a named numeric vector of +1 / -1.
resolve_direction <- function(direction, traits) {
  if (is.null(direction) || identical(direction, "auto")) {
    return(trait_directions(traits))
  }

  if (is.numeric(direction)) {
    if (is.null(names(direction))) {
      ps_abort("A numeric `direction` must be named after the trait columns.")
    }
    if (!all(direction %in% c(-1, 1))) {
      ps_abort("A numeric `direction` must contain only -1 (higher is better) or 1 (lower is better).")
    }
    direction <- ifelse(direction == 1, "lower_is_better", "higher_is_better")
  }

  if (!is.character(direction)) {
    ps_abort('`direction` must be "auto", a named character vector, or a named numeric vector.')
  }

  if (length(direction) == 1L && is.null(names(direction))) {
    direction <- stats::setNames(rep(direction, length(traits)), traits)
  }
  if (is.null(names(direction))) {
    ps_abort("`direction` must be named after the trait columns.")
  }

  valid <- c("higher_is_better", "lower_is_better")
  if (!all(direction %in% valid)) {
    ps_abort(paste0(
      '`direction` values must be "higher_is_better" or "lower_is_better"; got: ',
      paste(unique(setdiff(direction, valid)), collapse = ", "), "."
    ))
  }

  missing_traits <- setdiff(traits, names(direction))
  if (length(missing_traits) > 0L) {
    filled <- trait_directions(missing_traits)
    direction <- c(direction, filled)
  }
  direction[traits]
}

# +1 keeps the raw sign, -1 flips it, so that SRI > 0 always means damage.
direction_sign <- function(direction) {
  ifelse(direction == "lower_is_better", 1, -1)
}
