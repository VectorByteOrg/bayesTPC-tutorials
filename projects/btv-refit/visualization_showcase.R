# visualization_showcase.R
# Comprehensive visualization showcase for BTV refitting results

# Load libraries
library(dplyr)
library(readr)
library(jsonlite)

# Set up plotting
par(ask = FALSE)  # Don't ask between plots

cat("=== BTV Refitting Visualization Showcase ===\n\n")

# 1. Individual Trait Fits - Show one example
cat("1. Individual Trait Fit: Larval Survival (p) - C. sonorensis\n")
cat("   Loading trait panels plot...\n")

# Load the trait panels plot
trait_plot_path <- "outputs/traits/sonorensis/p/trait_panels.png"
if (file.exists(trait_plot_path)) {
  cat("   ✓ Trait panels plot available at:", trait_plot_path, "\n")
} else {
  cat("   ✗ Trait panels plot not found\n")
}

# 2. Supplement Fig. A.2 - Relative Width Analysis
cat("\n2. Supplement Fig. A.2: Relative Width Analysis\n")
cat("   Loading supplement figure...\n")

supplement_plot_path <- "outputs/supplement_fig_a2/supplement_fig_a2_robust.png"
if (file.exists(supplement_plot_path)) {
  cat("   ✓ Supplement figure available at:", supplement_plot_path, "\n")
} else {
  cat("   ✗ Supplement figure not found\n")
}

# 3. Create additional visualizations from the data
cat("\n3. Creating additional visualizations from fitted data...\n")

# Load curve summaries for comparison
curve_files <- list.files("outputs/traits/sonorensis", pattern = "curve_summary.csv", 
                          recursive = TRUE, full.names = TRUE)

curves_data <- list()
for (file in curve_files) {
  trait <- basename(dirname(file))
  curves_data[[trait]] <- read_csv(file)
  cat("   Loaded curve data for trait:", trait, "\n")
}

# Create a comparison plot of all traits
if (length(curves_data) > 0) {
  cat("\n   Creating trait comparison plot...\n")
  
  # Set up the plot
  png("outputs/trait_comparison.png", width = 1200, height = 800, res = 150)
  
  par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
  
  colors <- c("blue", "red", "green", "purple")
  trait_names <- c("p" = "Survival", "b" = "Vector Competence", 
                   "nu" = "Development Rate", "rho" = "Development Time")
  
  for (i in seq_along(curves_data)) {
    trait <- names(curves_data)[i]
    data <- curves_data[[trait]]
    
    plot(data$T, data$median, type = "l", lwd = 2, col = colors[i],
         xlab = "Temperature (°C)", ylab = trait_names[trait],
         main = paste(trait_names[trait], "- C. sonorensis"),
         ylim = c(0, max(data$upper, na.rm = TRUE) * 1.1))
    
    # Add HPD bands
    polygon(c(data$T, rev(data$T)), 
            c(data$lower, rev(data$upper)), 
            col = adjustcolor(colors[i], alpha = 0.2), border = NA)
    
    # Add median line
    lines(data$T, data$median, lwd = 2, col = colors[i])
  }
  
  dev.off()
  cat("   ✓ Trait comparison plot saved to: outputs/trait_comparison.png\n")
}

# 4. Load and display relative width analysis results
cat("\n4. Relative Width Analysis Results\n")

rw_file <- "outputs/supplement_fig_a2/relative_widths_robust.rds"
if (file.exists(rw_file)) {
  rw_data <- readRDS(rw_file)
  
  cat("   Creating relative width summary plot...\n")
  
  png("outputs/relative_width_summary.png", width = 1000, height = 600, res = 150)
  
  par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
  
  # Individual trait relative widths
  plot(rw_data$Tgrid, rw_data$trait_RW$p, type = "l", lwd = 2, col = "blue",
       ylim = c(0, 1), xlab = "Temperature (°C)", ylab = "Relative Width",
       main = "Individual Trait Uncertainty")
  
  lines(rw_data$Tgrid, rw_data$trait_RW$b, lwd = 2, col = "red")
  lines(rw_data$Tgrid, rw_data$trait_RW$nu, lwd = 2, col = "green")
  lines(rw_data$Tgrid, rw_data$trait_RW$rho, lwd = 2, col = "purple")
  
  legend("topright", legend = c("p (survival)", "b (competence)", "ν (rate)", "ρ (time)"),
         col = c("blue", "red", "green", "purple"), lwd = 2, bty = "n")
  
  # Derived quantity relative widths
  if (length(rw_data$derived_RW) >= 2) {
    plot(rw_data$Tgrid, rw_data$derived_RW$p_nu_p, type = "l", lwd = 2, col = "blue",
         ylim = c(0, 1), xlab = "Temperature (°C)", ylab = "Relative Width",
         main = "Derived Quantity Uncertainty")
    
    lines(rw_data$Tgrid, rw_data$derived_RW$p_nu_nu, lwd = 2, col = "green")
    lines(rw_data$Tgrid, rw_data$derived_RW$b_rho_b, lwd = 2, col = "red")
    lines(rw_data$Tgrid, rw_data$derived_RW$b_rho_rho, lwd = 2, col = "purple")
    
    legend("topright", legend = c("p×ν (p varies)", "p×ν (ν varies)", 
                                 "b×ρ (b varies)", "b×ρ (ρ varies)"),
           col = c("blue", "green", "red", "purple"), lwd = 2, bty = "n")
  }
  
  dev.off()
  cat("   ✓ Relative width summary plot saved to: outputs/relative_width_summary.png\n")
}

# 5. Parameter posterior distributions
cat("\n5. Parameter Posterior Distributions\n")

# Load fit summary for parameter analysis
fit_summary_file <- "outputs/traits/sonorensis/p/fit_summary.json"
if (file.exists(fit_summary_file)) {
  fit_summary <- fromJSON(fit_summary_file)
  
  cat("   Creating parameter summary plot...\n")
  
  png("outputs/parameter_summary.png", width = 1000, height = 600, res = 150)
  
  par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
  
  # Extract parameter estimates
  params <- fit_summary$parameters
  
  # Create bar plot of parameter estimates
  param_names <- c("T_max", "T_min", "q", "sigma.sq")
  param_means <- sapply(param_names, function(p) params[[p]]$mean)
  param_sds <- sapply(param_names, function(p) params[[p]]$sd)
  
  barplot(param_means, names.arg = param_names, 
          main = "Parameter Estimates - Larval Survival (p)",
          ylab = "Value", col = "lightblue")
  
  # Add error bars
  arrows(1:4, param_means - param_sds, 1:4, param_means + param_sds, 
         angle = 90, code = 3, length = 0.1)
  
  # Temperature range plot
  temp_range <- c(min(param_means[c("T_min", "T_max")]), 
                  max(param_means[c("T_min", "T_max")]))
  plot(temp_range, c(0, 1), type = "n", 
       xlab = "Temperature (°C)", ylab = "",
       main = "Temperature Range")
  rect(param_means["T_min"], 0, param_means["T_max"], 1, 
       col = "lightblue", border = "blue")
  text(mean(temp_range), 0.5, "Optimal Range", cex = 1.2)
  
  # Convergence diagnostics
  if ("convergence" %in% names(fit_summary)) {
    rhat_values <- fit_summary$convergence$rhat
    plot(1:length(rhat_values), rhat_values, type = "b", pch = 16,
         xlab = "Parameter", ylab = "R-hat", main = "Convergence (R-hat)",
         ylim = c(0.9, max(rhat_values) * 1.1))
    abline(h = 1.1, col = "red", lty = 2)
    text(1:length(rhat_values), rhat_values, names(rhat_values), 
         pos = 3, cex = 0.8)
  }
  
  # Model fit summary
  if ("fit_quality" %in% names(fit_summary)) {
    quality <- fit_summary$fit_quality
    plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "",
         main = "Model Fit Quality")
    text(1, 1, paste("RMSE:", round(quality$rmse, 4), "\n",
                     "R²:", round(quality$r_squared, 4), "\n",
                     "Data points:", quality$n_data), 
         cex = 1.2)
  }
  
  dev.off()
  cat("   ✓ Parameter summary plot saved to: outputs/parameter_summary.png\n")
}

# 6. Data overview
cat("\n6. Data Overview\n")

# Load tidied data
data_file <- "outputs/tidy_btv.csv"
if (file.exists(data_file)) {
  btv_data <- read_csv(data_file)
  
  cat("   Creating data overview plot...\n")
  
  png("outputs/data_overview.png", width = 1200, height = 800, res = 150)
  
  par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
  
  # Trait distribution
  trait_counts <- table(btv_data$trait)
  barplot(trait_counts, main = "Data Distribution by Trait",
          xlab = "Trait", ylab = "Number of Observations", col = "lightblue")
  
  # Temperature range by trait
  temp_ranges <- btv_data %>% 
    group_by(trait) %>% 
    summarise(min_T = min(T), max_T = max(T), .groups = "drop")
  
  plot(temp_ranges$min_T, 1:nrow(temp_ranges), type = "n",
       xlim = c(min(temp_ranges$min_T), max(temp_ranges$max_T)),
       ylim = c(0.5, nrow(temp_ranges) + 0.5),
       xlab = "Temperature (°C)", ylab = "", yaxt = "n",
       main = "Temperature Range by Trait")
  
  for (i in 1:nrow(temp_ranges)) {
    lines(c(temp_ranges$min_T[i], temp_ranges$max_T[i]), c(i, i), lwd = 3, col = "blue")
    points(temp_ranges$min_T[i], i, pch = 16, col = "red")
    points(temp_ranges$max_T[i], i, pch = 16, col = "red")
  }
  axis(2, at = 1:nrow(temp_ranges), labels = temp_ranges$trait, las = 2)
  
  # Value ranges by trait
  value_ranges <- btv_data %>% 
    group_by(trait) %>% 
    summarise(min_y = min(y), max_y = max(y), .groups = "drop")
  
  plot(value_ranges$min_y, 1:nrow(value_ranges), type = "n",
       xlim = c(min(value_ranges$min_y), max(value_ranges$max_y)),
       ylim = c(0.5, nrow(value_ranges) + 0.5),
       xlab = "Trait Value", ylab = "", yaxt = "n",
       main = "Value Range by Trait")
  
  for (i in 1:nrow(value_ranges)) {
    lines(c(value_ranges$min_y[i], value_ranges$max_y[i]), c(i, i), lwd = 3, col = "green")
    points(value_ranges$min_y[i], i, pch = 16, col = "red")
    points(value_ranges$max_y[i], i, pch = 16, col = "red")
  }
  axis(2, at = 1:nrow(value_ranges), labels = value_ranges$trait, las = 2)
  
  # Species comparison
  species_counts <- table(btv_data$species)
  barplot(species_counts, main = "Data Distribution by Species",
          xlab = "Species", ylab = "Number of Observations", col = "lightgreen")
  
  dev.off()
  cat("   ✓ Data overview plot saved to: outputs/data_overview.png\n")
}

cat("\n=== Visualization Summary ===\n")
cat("Generated plots:\n")
cat("1. Trait comparison: outputs/trait_comparison.png\n")
cat("2. Relative width summary: outputs/relative_width_summary.png\n")
cat("3. Parameter summary: outputs/parameter_summary.png\n")
cat("4. Data overview: outputs/data_overview.png\n")
cat("5. Individual trait panels: outputs/traits/sonorensis/*/trait_panels.png\n")
cat("6. Supplement Fig. A.2: outputs/supplement_fig_a2/supplement_fig_a2_robust.png\n")

cat("\nAll visualizations have been generated successfully!\n")
