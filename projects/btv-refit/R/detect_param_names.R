# detect_param_names.R
# Detect parameter names that bayesTPC actually uses

library(nimble)
library(bayesTPC)
library(dplyr)
library(readr)
library(coda)

# Function to detect parameter names from a fit
detect_param_names <- function(fit) {
  S <- fit$samples
  M <- if (inherits(S, "mcmc.list")) do.call(rbind, lapply(S, as.matrix))
       else if (is.list(S)) do.call(rbind, S)
       else if (is.matrix(S)) S
       else stop("Unknown fit$samples structure")
  
  # keep unique column names in original order
  unique(colnames(M))
}

# Load data
btv <- read_csv("outputs/tidy_btv.csv")

# Test 1: Brière + Normal
cat("=== Testing Brière + Normal ===\n")
one_briere <- subset(btv, species == "sonorensis" & trait == "p")
cat("Data points:", nrow(one_briere), "\n")

fit_briere_normal <- b_TPC(
  data = list(Temp = one_briere$T, Trait = one_briere$y),
  model = "briere",
  priors = NULL,  # default priors
  niter = 500,
  burn = 100,
  nchains = 1
)

briere_normal_params <- detect_param_names(fit_briere_normal)
cat("Brière + Normal parameters:", paste(briere_normal_params, collapse = ", "), "\n\n")

# Test 2: Quadratic + Normal
cat("=== Testing Quadratic + Normal ===\n")
one_quad <- subset(btv, species == "sonorensis" & trait == "rho")
cat("Data points:", nrow(one_quad), "\n")

fit_quad_normal <- b_TPC(
  data = list(Temp = one_quad$T, Trait = one_quad$y),
  model = "quadratic",
  priors = NULL,  # default priors
  niter = 500,
  burn = 100,
  nchains = 1
)

quad_normal_params <- detect_param_names(fit_quad_normal)
cat("Quadratic + Normal parameters:", paste(quad_normal_params, collapse = ", "), "\n\n")

# Test 3: Brière + Binomial (if we have success/trials data)
cat("=== Testing Brière + Binomial ===\n")
# Check if we have binomial data
one_binom <- subset(btv, species == "sonorensis" & trait == "b")
cat("Data points:", nrow(one_binom), "\n")
cat("Has replicates column:", !all(is.na(one_binom$replicates)), "\n")

if (!all(is.na(one_binom$replicates))) {
  # Create binomial data structure
  binom_data <- list(
    Temp = one_binom$T,
    Trait = one_binom$y,
    n_succ = round(one_binom$y * one_binom$replicates),
    n_trials = one_binom$replicates
  )
  
  fit_briere_binom <- b_TPC(
    data = binom_data,
    model = "briere",
    priors = NULL,  # default priors
    niter = 500,
    burn = 100,
    nchains = 1
  )
  
  briere_binom_params <- detect_param_names(fit_briere_binom)
  cat("Brière + Binomial parameters:", paste(briere_binom_params, collapse = ", "), "\n\n")
} else {
  cat("No binomial data available for testing\n\n")
}

# Summary
cat("=== SUMMARY ===\n")
cat("Brière + Normal:", paste(briere_normal_params, collapse = ", "), "\n")
cat("Quadratic + Normal:", paste(quad_normal_params, collapse = ", "), "\n")
if (!all(is.na(one_binom$replicates))) {
  cat("Brière + Binomial:", paste(briere_binom_params, collapse = ", "), "\n")
}
