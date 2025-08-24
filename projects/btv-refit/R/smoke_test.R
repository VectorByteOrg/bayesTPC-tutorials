# smoke_test.R
# Simple smoke test to verify bayesTPC pipeline

library(nimble)
library(bayesTPC)
library(dplyr)
library(readr)

# Load the tidy data
btv <- read_csv("outputs/tidy_btv.csv")

# Use a small subset for testing
one <- subset(btv, species == "sonorensis" & trait == "p")
stopifnot(nrow(one) > 3)

# Sanity check
cat("Data summary:\n")
cat("  Temperature range:", range(one$T), "\n")
cat("  y range:", range(one$y), "\n")
cat("  n observations:", nrow(one), "\n")

# Simple smoke test fit
USE_MCMC <- TRUE

fit <- tryCatch({
  if (!USE_MCMC) stop("MCMC disabled by flag")
  
  cat("Running smoke test fit...\n")
  
  b_TPC(
    data = list(Temp = one$T, Trait = one$y),
    model = "briere",
    priors = NULL,  # Use default priors to learn parameter names
    niter = 1000,
    burn = 100,
    nchains = 1
  )
}, error = function(e) {
  message("Fit failed: ", e$message)
  NULL
})

if (is.null(fit)) {
  stop("Smoke test fit failed; see message above.")
}

cat("Smoke test fit successful!\n")
cat("Fit object structure:\n")
print(names(fit))

# Extract draws to learn parameter names
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
cat("Parameter names in fit:\n")
print(names(draws))
cat("Number of posterior samples:", nrow(draws), "\n")

cat("Smoke test completed successfully!\n")
