# 02_fit_trait.R
# Fit thermal performance curves using bayesTPC

library(nimble)  # Load nimble first
library(bayesTPC)
library(dplyr)
library(purrr)
library(coda)
library(yaml)

# MCMC control flag
USE_MCMC <- TRUE

# Mean functions
briere_mean <- function(T, q, T_min, T_max) {
  out <- q * T * (T - T_min) * sqrt(pmax(T_max - T, 0))
  out[T < T_min | T > T_max] <- 0
  out
}

quad_mean <- function(T, inter, n.slope, qd) {
  inter - n.slope * T + qd * T^2
}

# HPD band helper
hpd_band <- function(curves, prob = 0.95) {
  t(apply(curves, 1, function(x) {
    h <- HPDinterval(as.mcmc(x), prob = prob)
    c(lower = h[1], upper = h[2])
  }))
}

# Load priors
priors_lib <- yaml::read_yaml("data/priors_tableA2.yaml")

# Model selection
choose_model <- function(trait) {
  if (trait %in% c("rhoE", "rhoL", "rhoP", "mu")) "quadratic" else "briere"
}

# Likelihood selection
choose_lik <- function(trait, df) {
  if (trait %in% c("pE", "pL", "pP", "b") && 
      "replicates" %in% names(df) &&
      !any(is.na(df$replicates))) {
    "binomial"
  } else if (trait %in% c("pE", "pL", "pP", "b")) {
    "truncnorm"
  } else {
    "normal"
  }
}

# Main fitting function
fit_trait <- function(df, trait, species) {
  # Filter data for species and trait
  df <- df %>% filter(species == species, trait == trait)
  
  if (nrow(df) == 0) {
    cat("No data found for", trait, "in", species, "\n")
    return(NULL)
  }
  
  # Check if priors exist for this trait
  if (is.null(priors_lib[[trait]])) {
    cat("No priors found for trait:", trait, "\n")
    return(NULL)
  }
  
  # Select model and likelihood
  like <- choose_lik(trait, df)
  model <- choose_model(trait)
  pri <- priors_lib[[trait]]
  
  cat("Fitting", trait, "for", species, "with", model, "model and", like, "likelihood\n")
  
  # Temperature grid for predictions
  Tgrid <- seq(min(df$T), max(df$T), length.out = 400)
  
  # Fit model with error handling
  fit <- tryCatch({
    if (!USE_MCMC) stop("MCMC disabled by flag")
    
    cat("  Data summary: T range [", min(df$T), ",", max(df$T), "], y range [", min(df$y), ",", max(df$y), "]\n")
    cat("  Model:", model, "| Likelihood:", like, "\n")
    
    # Prepare data in the format expected by bayesTPC
    data_list <- list(Temp = df$T, Trait = df$y)
    
    # Use translation system for custom priors
    if (!is.null(pri)) {
      cat("  Using custom priors from YAML\n")
      source("R/priors_translate.R")
      model_key <- model_key_of(model, like)
      
      # Detect parameter names with a quick fit
      tmp_fit <- b_TPC(
        data = data_list,
        model = model,
        priors = NULL,  # default priors for detection
        niter = 100,
        burn = 50,
        nchains = 1
      )
      
      # Extract parameter names
      S <- tmp_fit$samples
      if (inherits(S, "mcmc.list")) {
        detected_params <- unique(colnames(do.call(rbind, lapply(S, as.matrix))))
      } else if (inherits(S, "mcmc")) {
        detected_params <- colnames(S)
      } else {
        stop("Unknown samples structure")
      }
      
      # Translate priors
      pri_trans <- translate_priors_to_bayesTPC(pri, model_key)
      
      # Validate
      validate_prior_names(pri_trans, detected_params)
      
      cat("  Translated priors:", paste(names(pri_trans), collapse = ", "), "\n")
      
      # Final fit with translated priors
      b_TPC(
        data = data_list,
        model = model,
        priors = pri_trans,
        niter = 25000,
        burn = 5000,
        nchains = 5
      )
    } else {
      cat("  Using default priors\n")
      b_TPC(
        data = data_list,
        model = model,
        priors = NULL,
        niter = 25000,
        burn = 5000,
        nchains = 5
      )
    }
  }, error = function(e) {
    writeLines(
      paste(Sys.time(), trait, species, "ERROR:", e$message),
      con = "outputs/logs/fitting_errors.log", 
      sep = "\n"
    )
    cat("Error fitting", trait, "for", species, ":", e$message, "\n")
    cat("  Full error details:", toString(e), "\n")
    NULL
  })
  
  if (is.null(fit)) return(NULL)
  
  # Extract posterior draws
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
  
  # Source parameter utilities
  source("R/param_utils.R")
  
  # Determine model key and validate parameters
  model_key <- get_model_key(model, like)
  pd <- require_params(draws, model_key)
  
  # Generate curves
  curves <- switch(model_key,
    briere_normal = apply(pd, 1, function(p) {
      briere_mean(Tgrid, p["q"], p["T_min"], p["T_max"])
    }),
    briere_binom = apply(pd, 1, function(p) {
      briere_mean(Tgrid, p["q"], p["T_min"], p["T_max"])
    }),
    quadratic_normal = apply(pd, 1, function(p) {
      p["inter"] - p["n.slope"] * Tgrid + p["qd"] * Tgrid^2
    })
  )
  curves <- matrix(curves, nrow = length(Tgrid))
  
  # Calculate summaries
  result <- list(
    fit = fit,
    Tgrid = Tgrid,
    median = apply(curves, 1, median),
    hpd = hpd_band(curves),
    draws = draws,
    data = df,
    model = model,
    likelihood = like
  )
  
  # Print diagnostics
  cat("Fit complete for", trait, "in", species, "\n")
  cat("  Model:", model, "| Likelihood:", like, "\n")
  cat("  Data points:", nrow(df), "\n")
  cat("  Posterior samples:", nrow(draws), "\n")
  
  return(result)
}
