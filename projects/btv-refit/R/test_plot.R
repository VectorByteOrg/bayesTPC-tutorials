# test_plot.R
# Test plotting with fitted data

library(nimble)
library(bayesTPC)
library(dplyr)
library(readr)
library(coda)

# Load the tidy data
btv <- read_csv("outputs/tidy_btv.csv")

# Use a small subset for testing
one <- subset(btv, species == "sonorensis" & trait == "p")
cat("Data points:", nrow(one), "\n")
cat("Temperature range:", range(one$T), "\n")
cat("y range:", range(one$y), "\n")

# Simple fit
fit <- b_TPC(
  data = list(Temp = one$T, Trait = one$y),
  model = "briere",
  priors = NULL,
  niter = 1000,
  burn = 100,
  nchains = 1
)

# Extract draws
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
cat("Parameter names:", names(draws), "\n")
cat("Number of samples:", nrow(draws), "\n")

# Generate curves
briere_mean <- function(T, q, T_min, T_max) {
  q * T * (T - T_min) * sqrt(pmax(T_max - T, 0))
}

Tgrid <- seq(min(one$T), max(one$T), length.out = 300)
curves <- apply(draws, 1, function(p) {
  briere_mean(Tgrid, p["q"], p["T_min"], p["T_max"])
})
curves <- matrix(curves, nrow = length(Tgrid))

# Calculate summaries
med <- apply(curves, 1, median)
hpd <- t(apply(curves, 1, function(x) HPDinterval(as.mcmc(x), prob = 0.95)))

# Simple plot
plot(Tgrid, med, type = "l", lwd = 2, col = "blue",
     xlab = "Temperature (°C)", ylab = "p(T)",
     main = "Posterior median + 95% HPD for p (sonorensis)")
lines(Tgrid, hpd[, 1], lty = 2, col = "blue")
lines(Tgrid, hpd[, 2], lty = 2, col = "blue")
points(one$T, one$y, pch = 16, col = "red", cex = 1.2)

cat("Test plot completed successfully!\n")
