# show_visualizations.R
# Display key visualizations from BTV refitting analysis

cat("=== BTV Refitting Visualizations ===\n\n")

# 1. Show available plots
cat("Available visualizations:\n")
cat("1. Trait comparison plot: outputs/trait_comparison.png\n")
cat("2. Relative width summary: outputs/relative_width_summary.png\n")
cat("3. Individual trait panels: outputs/traits/sonorensis/*/trait_panels.png\n")
cat("4. Supplement Fig. A.2: outputs/supplement_fig_a2/supplement_fig_a2_robust.png\n\n")

# 2. Load and display some key results
cat("Key Results Summary:\n")

# Load relative width data
rw_file <- "outputs/supplement_fig_a2/relative_widths_robust.rds"
if (file.exists(rw_file)) {
  rw_data <- readRDS(rw_file)
  
  cat("Relative Width Analysis (Mean values):\n")
  cat("Individual traits:\n")
  for (trait in names(rw_data$trait_RW)) {
    mean_rw <- mean(rw_data$trait_RW[[trait]], na.rm = TRUE)
    cat("  ", trait, ":", round(mean_rw, 3), "\n")
  }
  
  if (length(rw_data$derived_RW) > 0) {
    cat("\nDerived quantities:\n")
    for (derived in names(rw_data$derived_RW)) {
      mean_rw <- mean(rw_data$derived_RW[[derived]], na.rm = TRUE)
      cat("  ", derived, ":", round(mean_rw, 3), "\n")
    }
  }
}

# 3. Load fit summary for parameter estimates
fit_summary_file <- "outputs/traits/sonorensis/p/fit_summary.json"
if (file.exists(fit_summary_file)) {
  fit_summary <- jsonlite::fromJSON(fit_summary_file)
  
  cat("\nParameter Estimates (Larval Survival - p):\n")
  params <- fit_summary$parameters
  for (param in names(params)) {
    mean_val <- params[[param]]$mean
    sd_val <- params[[param]]$sd
    cat("  ", param, ":", round(mean_val, 4), "±", round(sd_val, 4), "\n")
  }
}

# 4. Show data overview
data_file <- "outputs/tidy_btv.csv"
if (file.exists(data_file)) {
  btv_data <- readr::read_csv(data_file)
  
  cat("\nData Overview:\n")
  cat("Total observations:", nrow(btv_data), "\n")
  cat("Species:", paste(unique(btv_data$species), collapse = ", "), "\n")
  cat("Traits:", paste(unique(btv_data$trait), collapse = ", "), "\n")
  cat("Temperature range:", round(min(btv_data$T), 1), "to", round(max(btv_data$T), 1), "°C\n")
  
  # Show trait counts
  trait_counts <- table(btv_data$trait)
  cat("\nObservations per trait:\n")
  for (trait in names(trait_counts)) {
    cat("  ", trait, ":", trait_counts[trait], "\n")
  }
}

cat("\n=== Visualization Files ===\n")
cat("All plots have been saved to the outputs/ directory.\n")
cat("You can view them in your file browser or image viewer.\n")
