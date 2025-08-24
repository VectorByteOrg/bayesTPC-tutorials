# param_utils.R
# Parameter naming and validation utilities

# Single source of truth for parameter names
param_map <- list(
  briere_normal    = c("T_min", "T_max", "q", "sigma.sq"),
  quadratic_normal = c("inter", "n.slope", "qd", "sigma.sq"),
  briere_binom     = c("T_min", "T_max", "q")  # no sigma.sq
)

# Validate posterior columns before plotting
require_params <- function(draws, model_key) {
  want <- param_map[[model_key]]
  if (is.null(want)) stop("Unknown model_key: ", model_key)
  
  miss <- setdiff(want, names(draws))
  if (length(miss)) stop("Missing params in posterior: ", paste(miss, collapse = ", "))
  
  draws[, want, drop = FALSE]
}

# Determine model key from model and likelihood
get_model_key <- function(model, likelihood) {
  if (model == "briere" && likelihood == "normal") return("briere_normal")
  if (model == "quadratic" && likelihood == "normal") return("quadratic_normal")
  if (model == "briere" && likelihood == "binomial") return("briere_binom")
  if (model == "briere" && likelihood == "truncnorm") return("briere_normal")  # Same as normal
  
  stop("Unknown model+likelihood combination: ", model, "+", likelihood)
}
