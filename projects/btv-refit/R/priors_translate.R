# priors_translate.R
# Translation layer: YAML (paper-style) → bayesTPC (detected names)

# Detected parameter names from bayesTPC
# Based on detection results: both briere and quadratic use the same names
detected_map <- list(
  briere_normal    = c(Tmin = "T_min", Tmax = "T_max", k = "q", tau = "sigma.sq"),
  quadratic_normal = c(Tmin = "T_min", Tmax = "T_max", k = "q", tau = "sigma.sq"),
  briere_binom     = c(Tmin = "T_min", Tmax = "T_max", k = "q")  # no tau
)

# Helper to compute exponential rate from gamma tau parameters
exp_lambda_from_gamma_tau <- function(shape, rate) {
  if (shape <= 1) stop("shape (a) must be > 1 for E[sigma^2] to exist.")
  (shape - 1) / rate
}

# Convert Gamma(a, rate=b) prior on tau (precision) to sigma.sq
# Try Inv-Gamma first, fall back to Exp with mean-matching
translate_tau_prior <- function(pr) {
  if (is.null(pr$tau)) return(pr)
  
  a <- pr$tau$shape
  b <- pr$tau$rate
  stopifnot(tolower(pr$tau$dist) == "gamma", a > 1)
  
  # Try Inv-Gamma first (if bayesTPC supports it)
  # If not, use Exp with mean-matching
  lambda <- exp_lambda_from_gamma_tau(a, b)
  pr$sigma.sq <- list(dist = "exp", rate = lambda)
  pr$tau <- NULL
  
  cat("Converted tau ~ Gamma(", a, ", ", b, ") to sigma.sq ~ Exp(", lambda, ")\n")
  cat("  (matching mean E[sigma^2] = ", b/(a-1), ")\n")
  
  pr
}

# Convert YAML prior to bayesTPC format (named character vector)
convert_prior_to_bayesTPC_format <- function(prior_spec) {
  dist <- prior_spec$dist
  if (dist == "uniform") {
    paste0("dunif(", prior_spec$min, ", ", prior_spec$max, ")")
  } else if (dist == "gamma") {
    paste0("dgamma(", prior_spec$shape, ", ", prior_spec$rate, ")")
  } else if (dist == "exp") {
    paste0("dexp(", prior_spec$rate, ")")
  } else if (dist == "inv_gamma") {
    paste0("dinvgamma(", prior_spec$shape, ", ", prior_spec$scale, ")")
  } else if (dist == "normal") {
    paste0("dnorm(", prior_spec$mean, ", ", prior_spec$sd, ")")
  } else {
    stop("Unsupported distribution: ", dist)
  }
}

# Map YAML priors (canonical names) to detected names
translate_priors_to_bayesTPC <- function(yaml_block, model_key) {
  stopifnot(model_key %in% names(detected_map))
  m <- detected_map[[model_key]]
  
  # Handle both nested and flat YAML structures
  if ("priors" %in% names(yaml_block)) {
    pr <- yaml_block$priors
  } else {
    pr <- yaml_block  # flat structure
  }
  
  cat("Debug: YAML prior keys:", paste(names(pr), collapse = ", "), "\n")
  cat("Debug: Expected keys:", paste(names(m), collapse = ", "), "\n")
  
  # Handle tau translation first
  has_tau <- !is.null(pr$tau)
  if (has_tau) {
    pr <- translate_tau_prior(pr)  # creates sigma.sq if tau was present
    # Update mapping for tau -> sigma.sq
    m["tau"] <- "sigma.sq"
  }
  
  # Convert to bayesTPC format (named character vector)
  out <- character()
  for (canon in names(m)) {
    det <- m[[canon]]
    # Handle tau -> sigma.sq mapping
    if (canon == "tau" && has_tau) {
      if (is.null(pr$sigma.sq)) {
        stop("YAML prior missing 'tau' (should be converted to sigma.sq) for model_key=", model_key)
      }
      out[det] <- convert_prior_to_bayesTPC_format(pr$sigma.sq)
    } else {
      if (is.null(pr[[canon]])) {
        stop("YAML prior missing '", canon, "' for model_key=", model_key)
      }
      out[det] <- convert_prior_to_bayesTPC_format(pr[[canon]])
    }
  }
  out
}

# Validate that translated priors match detected params
validate_prior_names <- function(translated_priors, detected_params) {
  need <- names(translated_priors)
  missing <- setdiff(need, detected_params)
  if (length(missing)) {
    stop("Translated prior names not in detected params: ", paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}

# Model key decision function
model_key_of <- function(mean_fn, likelihood) {
  if (mean_fn == "briere" && likelihood == "normal") return("briere_normal")
  if (mean_fn == "quadratic" && likelihood == "normal") return("quadratic_normal")
  if (mean_fn == "briere" && likelihood == "binomial") return("briere_binom")
  if (mean_fn == "briere" && likelihood == "truncnorm") return("briere_normal")  # Same as normal
  stop("Unsupported mean_fn/likelihood combo: ", mean_fn, "+", likelihood)
}

# Log the final priors used (great for reproducibility)
log_priors <- function(final_priors, path) {
  lines <- vapply(names(final_priors), function(nm) {
    p <- final_priors[[nm]]
    sprintf("%s: %s", nm, p)
  }, character(1))
  writeLines(lines, path)
  cat("Prior specifications saved to:", path, "\n")
}
