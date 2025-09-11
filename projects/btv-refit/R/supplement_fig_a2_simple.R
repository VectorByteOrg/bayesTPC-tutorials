# supplement_fig_a2_simple.R
# Simplified Supplement Fig. A.2: Relative width analysis for available traits

# Load helpers and packages
library(coda)
library(dplyr)
library(purrr)
library(readr)

# Setup
sp <- "sonorensis"
root <- file.path("outputs", "traits", sp)

# 1) Utilities

# Mean functions
briere_mean <- function(T, q, T_min, T_max) {
  out <- q * T * (T - T_min) * sqrt(pmax(T_max - T, 0))
  out[T < T_min | T > T_max] <- 0
  out
}

# Build posterior curve matrix (nT x nDraws) from a fit
posterior_curves_from_fit <- function(fit, Tgrid) {
  draws <- {
    S <- fit$samples
    if (inherits(S, "mcmc.list")) do.call(rbind, lapply(S, as.matrix))
    else if (is.list(S)) do.call(rbind, S)
    else if (is.matrix(S)) S
    else stop("Unknown sample structure")
  } %>% as.data.frame()
  
  # All our traits use briere model
  need <- c("q", "T_min", "T_max")
  stopifnot(all(need %in% names(draws)))
  mat <- apply(draws[need], 1, \(p) briere_mean(Tgrid, p["q"], p["T_min"], p["T_max"]))
  
  matrix(mat, nrow = length(Tgrid))
}

# 95% relative width helper
rel_width <- function(samples_matrix) {
  qs <- apply(samples_matrix, 1, quantile, probs = c(.025, .5, .975), na.rm = TRUE)
  (qs[3, ] - qs[1, ]) / pmax(qs[2, ], .Machine$double.eps)
}

# 2) Load fits and build posterior curves

# Helper to load a fit and return posterior curves
load_trait_curves <- function(trait, Tgrid) {
  fit <- readRDS(file.path(root, trait, "fit.rds"))
  posterior_curves_from_fit(fit, Tgrid)
}

# Temperature grid based on observed data range
Tgrid <- seq(12, 35, length.out = 300)

# Available traits
available_traits <- c("p", "b", "nu", "rho")

# Build posterior curve matrices for each available trait
curves <- list()
for (trait in available_traits) {
  cat("Loading", trait, "\n")
  tryCatch({
    curves[[trait]] <- load_trait_curves(trait, Tgrid)
  }, error = function(e) {
    cat("Error loading", trait, ":", e$message, "\n")
  })
}

# 3) Compute Relative Widths for each trait

cat("Computing relative widths...\n")

# Compute relative widths for each trait
trait_RW <- list()
for (trait in names(curves)) {
  trait_RW[[trait]] <- rel_width(curves[[trait]])
}

# 4) Create derived quantities using available traits

# For demonstration, let's create some derived quantities using the traits we have
# We'll create simplified versions of V, f, and g using available traits

# Derived quantity 1: p * nu (survival * development rate)
derive_p_nu <- function(p_curves, nu_curves) {
  nT <- nrow(p_curves)
  nD <- ncol(p_curves)
  out <- matrix(NA_real_, nT, nD)
  for (j in 1:nD) {
    out[, j] <- p_curves[, j] * nu_curves[, j]
  }
  out
}

# Derived quantity 2: b * rho (vector competence * development time)
derive_b_rho <- function(b_curves, rho_curves) {
  nT <- nrow(b_curves)
  nD <- ncol(b_curves)
  out <- matrix(NA_real_, nT, nD)
  for (j in 1:nD) {
    out[, j] <- b_curves[, j] * rho_curves[, j]
  }
  out
}

# Compute derived quantity relative widths
derived_RW <- list()

# p * nu when p varies, nu fixed at median
if ("p" %in% names(curves) && "nu" %in% names(curves)) {
  nu_med <- apply(curves$nu, 1, median)
  p_nu_p_var <- matrix(NA_real_, nrow(curves$p), ncol(curves$p))
  for (j in 1:ncol(curves$p)) {
    p_nu_p_var[, j] <- curves$p[, j] * nu_med
  }
  derived_RW$p_nu_p <- rel_width(p_nu_p_var)
}

# p * nu when nu varies, p fixed at median
if ("p" %in% names(curves) && "nu" %in% names(curves)) {
  p_med <- apply(curves$p, 1, median)
  p_nu_nu_var <- matrix(NA_real_, nrow(curves$nu), ncol(curves$nu))
  for (j in 1:ncol(curves$nu)) {
    p_nu_nu_var[, j] <- p_med * curves$nu[, j]
  }
  derived_RW$p_nu_nu <- rel_width(p_nu_nu_var)
}

# b * rho when b varies, rho fixed at median
if ("b" %in% names(curves) && "rho" %in% names(curves)) {
  rho_med <- apply(curves$rho, 1, median)
  b_rho_b_var <- matrix(NA_real_, nrow(curves$b), ncol(curves$b))
  for (j in 1:ncol(curves$b)) {
    b_rho_b_var[, j] <- curves$b[, j] * rho_med
  }
  derived_RW$b_rho_b <- rel_width(b_rho_b_var)
}

# b * rho when rho varies, b fixed at median
if ("b" %in% names(curves) && "rho" %in% names(curves)) {
  b_med <- apply(curves$b, 1, median)
  b_rho_rho_var <- matrix(NA_real_, nrow(curves$rho), ncol(curves$rho))
  for (j in 1:ncol(curves$rho)) {
    b_rho_rho_var[, j] <- b_med * curves$rho[, j]
  }
  derived_RW$b_rho_rho <- rel_width(b_rho_rho_var)
}

# 5) Plot

# Create output directory
outdir <- file.path("outputs", "supplement_fig_a2")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Save the plot
png(file.path(outdir, "supplement_fig_a2_simple.png"), width = 1200, height = 900, res = 150)

op <- par(no.readonly = TRUE)
on.exit(par(op))
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

# Top-left: Individual trait relative widths
plot(Tgrid, trait_RW$p, type = "l", lwd = 2, col = "blue",
     ylim = c(0, max(unlist(trait_RW), na.rm = TRUE)), 
     xlab = "Temperature (°C)", ylab = "Relative width of quantiles",
     main = "Individual Trait Uncertainty")
for (trait in names(trait_RW)) {
  col <- switch(trait, p = "blue", b = "red", nu = "green", rho = "purple")
  lines(Tgrid, trait_RW[[trait]], lwd = 2, col = col)
}
legend("topright", legend = names(trait_RW), 
       col = c("blue", "red", "green", "purple"), lwd = 2, bty = "n")

# Top-right: p * nu uncertainty attribution
if (length(derived_RW) >= 2) {
  plot(Tgrid, derived_RW$p_nu_p, type = "l", lwd = 2, col = "blue",
       ylim = c(0, max(c(derived_RW$p_nu_p, derived_RW$p_nu_nu), na.rm = TRUE)),
       xlab = "Temperature (°C)", ylab = "Relative width of quantiles",
       main = "p × ν Uncertainty Attribution")
  lines(Tgrid, derived_RW$p_nu_nu, lwd = 2, col = "green")
  legend("topright", legend = c("p varies", "ν varies"), 
         col = c("blue", "green"), lwd = 2, bty = "n")
}

# Bottom-left: b * rho uncertainty attribution
if (length(derived_RW) >= 4) {
  plot(Tgrid, derived_RW$b_rho_b, type = "l", lwd = 2, col = "red",
       ylim = c(0, max(c(derived_RW$b_rho_b, derived_RW$b_rho_rho), na.rm = TRUE)),
       xlab = "Temperature (°C)", ylab = "Relative width of quantiles",
       main = "b × ρ Uncertainty Attribution")
  lines(Tgrid, derived_RW$b_rho_rho, lwd = 2, col = "purple")
  legend("topright", legend = c("b varies", "ρ varies"), 
         col = c("red", "purple"), lwd = 2, bty = "n")
}

# Bottom-right: caption
plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "")
text(1, 1, "Supplement Fig. A.2 (Simplified)\nRelative width of 95% credible intervals\nfor available traits and derived quantities", 
     cex = 1.2, font = 2)

dev.off()

cat("Supplement Fig. A.2 (simplified) saved to:", file.path(outdir, "supplement_fig_a2_simple.png"), "\n")

# Save the computed relative widths for further analysis
saveRDS(list(
  Tgrid = Tgrid,
  trait_RW = trait_RW,
  derived_RW = derived_RW
), file.path(outdir, "relative_widths_simple.rds"))

cat("Relative widths data saved to:", file.path(outdir, "relative_widths_simple.rds"), "\n")

# Print summary statistics
cat("\nSummary of relative widths:\n")
cat("Individual traits - Mean relative widths:\n")
for (nm in names(trait_RW)) {
  cat("  ", nm, ":", round(mean(trait_RW[[nm]], na.rm = TRUE), 3), "\n")
}

if (length(derived_RW) > 0) {
  cat("\nDerived quantities - Mean relative widths:\n")
  for (nm in names(derived_RW)) {
    cat("  ", nm, ":", round(mean(derived_RW[[nm]], na.rm = TRUE), 3), "\n")
  }
}
