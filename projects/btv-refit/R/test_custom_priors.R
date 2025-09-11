# test_custom_priors.R
# Test custom priors with tau → sigma.sq translation

library(nimble)
library(bayesTPC)
library(dplyr)
library(readr)
library(yaml)

# Load utilities
source("R/param_utils.R")
source("R/prior_utils.R")
source("R/qc_utils.R")

# Load the tidy data
btv <- read_csv("outputs/tidy_btv.csv")

# Load priors
priors_lib <- yaml::read_yaml("data/priors_tableA2.yaml")

# Test with pL trait (should have custom priors)
trait <- "pL"
species <- "sonorensis"

# Get data
one <- subset(btv, species == species & trait == trait)
cat("Testing custom priors for", trait, "in", species, "\n")
cat("Data points:", nrow(one), "\n")

# Get and translate priors
pri <- priors_lib[[trait]]
cat("Original priors:\n")
print(pri)

# Translate tau → sigma.sq (if needed)
pri <- translate_tau_to_sigma2(pri)
cat("Translated priors:\n")
print(pri)

# Validate priors
model_key <- "briere_normal"  # pL should be Brière with normal likelihood
validate_priors(pri, model_key)
cat("Prior validation passed!\n")

# Test fit with custom priors
cat("Running fit with custom priors...\n")
fit <- b_TPC(
  data = list(Temp = one$T, Trait = one$y),
  model = "briere",
  priors = pri,
  niter = 1000,
  burn = 100,
  nchains = 1
)

cat("Fit completed successfully!\n")

# Extract draws and validate parameters
extract_draws_df <- function(fit) {
  S <- fit$samples
  if (inherits(S, "mcmc.list")) {
    as.data.frame(do.call(rbind, lapply(S, as.matrix)))
  } else if (is.list(S)) {
    as.data.frame(do.call(rbind, S))
  } else if (is.matrix(S)) {
    as.data.frame(S)
  } else {
    stop("Unknown samples structure in fit$samples")
  }
}

draws <- extract_draws_df(fit)
pd <- require_params(draws, model_key)

cat("Parameter names in fit:", names(draws), "\n")
cat("Validated parameters:", names(pd), "\n")
cat("Number of posterior samples:", nrow(draws), "\n")

# Check convergence
rhats <- check_convergence(fit)
if (!is.null(rhats)) {
  cat("R-hat values:\n")
  print(rhats)
}

cat("Custom priors test completed successfully!\n")
