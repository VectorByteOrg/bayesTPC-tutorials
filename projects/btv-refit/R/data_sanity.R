# data_sanity.R
# Data validation and transformation functions

assert_prob <- function(x, nm) {
  stopifnot(!any(is.na(x)))
  if (any(x < 0 | x > 1, na.rm = TRUE))
    stop(sprintf("%s must be in [0,1]; range = [%g,%g]", nm, min(x), max(x)))
}

assert_positive <- function(x, nm) {
  if (any(x <= 0, na.rm = TRUE))
    stop(sprintf("%s must be > 0; min = %g", nm, min(x, na.rm = TRUE)))
}

transform_btv <- function(df) {
  df <- df
  
  # probabilities: %. If >1, assume percent -> divide by 100
  prob_traits <- c("p", "pE", "pL", "pP", "b")
  df <- df |>
    dplyr::mutate(
      y = dplyr::case_when(
        trait %in% prob_traits & y > 1 ~ y/100,
        TRUE ~ y
      )
    )
  
  # ν = 1/EIP (already handled in the data processing)
  # μ: if your sheet had lifespan (days), convert upstream to mu = 1/lifespan
  # (If it's already a rate 1/day, leave as is.)
  
  df
}

validate_rowwise <- function(df) {
  for (tr in unique(df$trait)) {
    x <- df[df$trait == tr, "y", drop = TRUE]
    if (tr %in% c("pE", "pL", "pP", "b")) assert_prob(x, tr)
    if (tr %in% c("a", "F", "mu", "rhoE", "rhoL", "rhoP")) assert_positive(x, tr)
    # Special case for nu: can be 0 (no development at low temps)
    if (tr == "nu") {
      if (any(x < 0, na.rm = TRUE)) stop(sprintf("%s must be >= 0; min = %g", tr, min(x, na.rm = TRUE)))
    }
  }
  invisible(TRUE)
}
