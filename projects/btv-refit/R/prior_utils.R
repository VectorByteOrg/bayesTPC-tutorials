# prior_utils.R
# Prior translation and utility functions

# Convert tau (precision) to sigma.sq (variance) priors
# If YAML contains tau ~ Gamma(a, b), then sigma.sq ~ Inv-Gamma(a, b)
translate_tau_to_sigma2 <- function(priors_in) {
  pr <- priors_in
  if (!is.null(pr$tau)) {
    if (tolower(pr$tau$dist) != "gamma") {
      stop("tau prior must be Gamma(shape=a, rate=b)")
    }
    a <- pr$tau$shape
    b <- pr$tau$rate
    pr$sigma.sq <- list(dist = "inv_gamma", shape = a, scale = b)
    pr$tau <- NULL
  }
  pr
}

# Inverse gamma density function
dinvgamma <- function(x, shape, scale) {
  # shape=a, scale=b; density ∝ b^a x^{-(a+1)} exp(-b/x)
  ifelse(x > 0, (scale^shape / gamma(shape)) * x^(-(shape + 1)) * exp(-scale/x), 0)
}

# Validate prior structure matches expected parameters
validate_priors <- function(priors, model_key) {
  source("R/param_utils.R")
  
  need <- param_map[[model_key]]
  if (is.null(need)) stop("Unknown model_key: ", model_key)
  
  missing <- setdiff(need, names(priors))
  if (length(missing)) {
    stop("YAML prior missing keys: ", paste(missing, collapse = ", "))
  }
  
  invisible(TRUE)
}
