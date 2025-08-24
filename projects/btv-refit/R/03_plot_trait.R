# 03_plot_trait.R
# Plot thermal performance curves in Appendix A.4 style

library(dplyr)
library(coda)

# Main plotting function for trait panels
plot_trait_panels <- function(res, trait, pri) {
  old <- par(no.readonly = TRUE)
  on.exit(par(old))
  
  # Layout: top panel (curve), bottom panels (priors vs posteriors)
  layout(matrix(c(1, 1, 2, 3, 4, 5), nrow = 3, byrow = TRUE), heights = c(2, 1, 1))
  par(mar = c(4, 4, 2, 1))
  
  # Top panel: curve + HPD + data
  plot(res$Tgrid, res$median, type = "l", lwd = 2, col = "blue",
       xlab = "Temperature (°C)", ylab = paste0(trait, "(T)"),
       main = paste("Posterior median + 95% HPD for", trait))
  
  # Add HPD bands
  lines(res$Tgrid, res$hpd[, "lower"], lty = 2, col = "blue")
  lines(res$Tgrid, res$hpd[, "upper"], lty = 2, col = "blue")
  
  # Add data points
  points(res$data$T, res$data$y, pch = 16, col = "red", cex = 1.2)
  
  # Add legend
  legend("topright", 
         legend = c("Posterior median", "95% HPD", "Data"),
         col = c("blue", "blue", "red"),
         lty = c(1, 2, NA),
         pch = c(NA, NA, 16),
         lwd = c(2, 1, NA),
         bty = "n")
  
  # Source utilities
  source("R/param_utils.R")
  source("R/prior_utils.R")
  
  # Determine model key and get validated parameters
  model_key <- get_model_key(res$model, res$likelihood)
  pd <- require_params(res$draws, model_key)
  
  # Bottom panels: prior vs posterior for parameters
  if (model_key == "briere_normal" || model_key == "briere_binom") {
    
    # T_min
    hist(pd$T_min, breaks = 30, main = "T_min", xlab = "", prob = TRUE,
         col = "lightblue", border = "blue")
    # Note: priors are now in bayesTPC format (character strings), not YAML format
    
    # T_max  
    hist(pd$T_max, breaks = 30, main = "T_max", xlab = "", prob = TRUE,
         col = "lightblue", border = "blue")
    
    # q
    hist(pd$q, breaks = 30, main = "q", xlab = "", prob = TRUE,
         col = "lightblue", border = "blue")
    
    # sigma.sq (only for normal likelihood)
    if (model_key == "briere_normal") {
      hist(pd$sigma.sq, breaks = 30, main = "sigma.sq (log scale)", xlab = "", prob = TRUE,
           col = "lightblue", border = "blue", log = "y")
      # Note: priors are now in bayesTPC format (character strings), not YAML format
    }
  } else if (model_key == "quadratic_normal") {
    # Quadratic parameters
    hist(pd$inter, breaks = 30, main = "intercept", xlab = "", prob = TRUE,
         col = "lightblue", border = "blue")
    
    hist(pd$n.slope, breaks = 30, main = "linear slope", xlab = "", prob = TRUE,
         col = "lightblue", border = "blue")
    
    hist(pd$qd, breaks = 30, main = "quadratic", xlab = "", prob = TRUE,
         col = "lightblue", border = "blue")
    
    hist(pd$sigma.sq, breaks = 30, main = "sigma.sq (log scale)", xlab = "", prob = TRUE,
         col = "lightblue", border = "blue", log = "y")
  } else {
    # Fallback: just show posterior histograms
    param_names <- names(pd)
    for (i in 1:min(4, length(param_names))) {
      param <- param_names[i]
      hist(pd[[param]], breaks = 30, main = param, xlab = "", prob = TRUE,
           col = "lightblue", border = "blue")
    }
  }
}

# Save plot function
save_trait_plot <- function(res, trait, species, outdir) {
  # Load priors if available
  priors_lib <- NULL
  tryCatch({
    priors_lib <- yaml::read_yaml("data/priors_tableA2.yaml")
  }, error = function(e) {
    cat("Warning: Could not load priors for plotting\n")
  })
  
  png(file.path(outdir, "trait_panels.png"), width = 1200, height = 900, res = 150)
  if (!is.null(priors_lib) && trait %in% names(priors_lib)) {
    # Translate priors to bayesTPC format for plotting
    source("R/priors_translate.R")
    model_key <- model_key_of(res$model, res$likelihood)
    pri_trans <- translate_priors_to_bayesTPC(priors_lib[[trait]], model_key)
    plot_trait_panels(res, trait, pri_trans)
  } else {
    # Fallback: just plot the curve without priors
    plot_trait_curve(res, trait)
  }
  dev.off()
  
  cat("Plot saved to", file.path(outdir, "trait_panels.png"), "\n")
}

# Summary statistics function
get_fit_summary <- function(res) {
  # Parameter summaries
  param_summary <- summary(res$fit)$statistics
  
  # Convergence diagnostics (if available)
  convergence <- NULL
  if ("Rhat" %in% names(res$fit)) {
    convergence <- res$fit$Rhat
  }
  
  # Residuals at observed temperatures
  if (res$model == "briere") {
    pred_obs <- apply(res$draws, 1, function(th) {
      briere_mean(res$data$T, th["k"], th["Tmin"], th["Tmax"])
    })
  } else {
    pred_obs <- apply(res$draws, 1, function(th) {
      quad_mean(res$data$T, th["inter"], th["n.slope"], th["qd"])
    })
  }
  
  pred_median <- apply(pred_obs, 1, median)
  residuals <- res$data$y - pred_median
  
  list(
    parameters = param_summary,
    convergence = convergence,
    residuals = residuals,
    rmse = sqrt(mean(residuals^2)),
    n_data = nrow(res$data)
  )
}

# Fallback plotting function (just the curve, no priors)
plot_trait_curve <- function(res, trait) {
  # Top panel: curve with data
  plot(res$Tgrid, res$median, type = "l", lwd = 2, col = "blue",
       xlab = "Temperature (°C)", ylab = trait,
       main = paste("Thermal Performance Curve -", trait))
  lines(res$Tgrid, res$hpd[, "lower"], lty = 2, col = "blue")
  lines(res$Tgrid, res$hpd[, "upper"], lty = 2, col = "blue")
  points(res$data$T, res$data$y, pch = 16, col = "red")
  legend("topright", legend = c("Posterior median", "95% HPD", "Data"),
         lty = c(1, 2, NA), pch = c(NA, NA, 16), col = c("blue", "blue", "red"),
         bty = "n")
}
