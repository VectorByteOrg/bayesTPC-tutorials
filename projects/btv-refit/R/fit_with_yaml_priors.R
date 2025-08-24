# fit_with_yaml_priors.R
# Fitting function with YAML prior translation

library(nimble)
library(bayesTPC)
library(dplyr)
library(readr)
library(yaml)

# Source utilities
source("R/priors_translate.R")
source("R/qc_utils.R")

# Function to detect parameter names from a fit
detect_param_names <- function(fit) {
  S <- fit$samples
  M <- if (inherits(S, "mcmc.list")) do.call(rbind, lapply(S, as.matrix))
       else if (is.list(S)) do.call(rbind, S)
       else if (is.matrix(S)) S
       else stop("Unknown fit$samples structure")
  
  unique(colnames(M))
}

# Main fitting function with YAML prior translation
fit_with_yaml_priors <- function(T, y, mean_fn, likelihood, yaml_priors_block, 
                                niter = 10000, burn = 2000, nchains = 3, thin = 2) {
  model_key <- model_key_of(mean_fn, likelihood)
  
  cat("Fitting", mean_fn, "+", likelihood, "model\n")
  cat("Model key:", model_key, "\n")
  
  # 1) Do a tiny default fit to detect names reliably
  cat("Detecting parameter names...\n")
  tmp <- b_TPC(
    data = list(Temp = T, Trait = y),
    model = mean_fn,
    priors = NULL,  # default priors
    niter = 100,
    burn = 50,
    nchains = 1
  )
  detected <- detect_param_names(tmp)
  cat("Detected parameters:", paste(detected, collapse = ", "), "\n")
  
  # 2) Translate YAML → bayesTPC names
  cat("Translating priors...\n")
  pri_trans <- translate_priors_to_bayesTPC(yaml_priors_block, model_key)
  cat("Translated prior names:", paste(names(pri_trans), collapse = ", "), "\n")
  
  # 3) Validate keys
  validate_prior_names(pri_trans, detected)
  cat("Prior validation passed!\n")
  
  # 4) Final fit with custom priors
  cat("Running final fit with custom priors...\n")
  fit <- b_TPC(
    data = list(Temp = T, Trait = y),
    model = mean_fn,
    priors = pri_trans,
    niter = niter,
    burn = burn,
    nchains = nchains,
    thin = thin
  )
  
  cat("Fit completed successfully!\n")
  fit
}

# Test function
test_yaml_priors <- function() {
  # Load data and priors
  btv <- read_csv("outputs/tidy_btv.csv")
  priors_lib <- yaml::read_yaml("data/priors_tableA2.yaml")
  
  # Test with pL trait
  trait <- "pL"
  species <- "sonorensis"
  
  # Get data
  one <- subset(btv, species == species & trait == trait)
  cat("Testing YAML priors for", trait, "in", species, "\n")
  cat("Data points:", nrow(one), "\n")
  
  # Get YAML block
  pl_block <- priors_lib[[trait]]
  cat("YAML block structure:\n")
  print(str(pl_block))
  
  # Run fit with translation
  fit <- fit_with_yaml_priors(
    T = one$T, 
    y = one$y,
    mean_fn = "briere", 
    likelihood = "normal",
    yaml_priors_block = pl_block,
    niter = 1000,  # small for testing
    burn = 200,
    nchains = 1
  )
  
  # Check results
  draws <- detect_param_names(fit)
  cat("Final fit parameter names:", paste(draws, collapse = ", "), "\n")
  cat("Number of posterior samples:", nrow(as.data.frame(fit$samples)), "\n")
  
  # Check convergence
  rhats <- check_convergence(fit)
  if (!is.null(rhats)) {
    cat("R-hat values:\n")
    print(rhats)
  }
  
  cat("YAML priors test completed successfully!\n")
  fit
}
