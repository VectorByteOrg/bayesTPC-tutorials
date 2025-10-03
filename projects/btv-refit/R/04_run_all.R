# 04_run_all.R
# Main script to run all trait fits

# Source all scripts
source("R/01_load_tidy.R")
source("R/02_fit_trait.R") 
source("R/03_plot_trait.R")
source("R/qc_utils.R")

# Load tidied data
btv <- readr::read_csv("outputs/tidy_btv.csv")

# Define species and traits to fit
species <- c("sonorensis", "variipennis")
traits <- c("p", "b", "nu", "rho")  # Only traits available in the data

# Create output directories
dir.create("outputs/traits", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/logs", showWarnings = FALSE, recursive = TRUE)

# Save session info
save_session_info("outputs")

# Initialize log file
writeLines(paste("BTV Refit Log - Started at", Sys.time()), 
           con = "outputs/logs/refit_log.txt")

# Loop through species and traits
for (sp in species) {
  for (tr in traits) {
    cat("\n", paste(rep("=", 50), collapse = ""), "\n")
    cat("Processing:", tr, "for", sp, "\n")
    cat(paste(rep("=", 50), collapse = ""), "\n")
    
    # Check if data exists for this combination
    df <- btv %>% filter(species == sp, trait == tr)
    
    if (nrow(df) == 0) {
      cat("No data found for", tr, "in", sp, "- skipping\n")
      next
    }
    
    # Fit the trait
    res <- fit_trait(btv, tr, sp)
    
    if (is.null(res)) {
      cat("Fit failed for", tr, "in", sp, "- skipping\n")
      next
    }
    
    # Create output directory
    outdir <- file.path("outputs/traits", sp, tr)
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    
    # Save fit object
    saveRDS(res$fit, file.path(outdir, "fit.rds"))
    
    # Generate QC summary
    qc_summary <- generate_qc_summary(res, tr, sp)
    save_qc_summary(qc_summary, outdir)
    
    # Save curve summary
    curve_summary <- tibble(
      T = res$Tgrid,
      median = res$median,
      lower = res$hpd[, "lower"],
      upper = res$hpd[, "upper"]
    )
    readr::write_csv(curve_summary, file.path(outdir, "curve_summary.csv"))
    
    # Generate and save plot
    save_trait_plot(res, tr, sp, outdir)
    
    # Get and save fit summary
    summary_stats <- get_fit_summary(res)
    
    # Save summary as JSON
    jsonlite::write_json(summary_stats, file.path(outdir, "fit_summary.json"), 
                        pretty = TRUE, auto_unbox = TRUE)
    
    # Log success
    log_entry <- paste(Sys.time(), "SUCCESS:", tr, sp, 
                      "Data points:", nrow(df),
                      "RMSE:", round(summary_stats$rmse, 4))
    cat(log_entry, "\n", file = "outputs/logs/refit_log.txt", append = TRUE)
    
    cat("✓ Completed", tr, "for", sp, "\n")
    cat("  Outputs saved to:", outdir, "\n")
  }
}

# Final summary
cat("\n", paste(rep("=", 50), collapse = ""), "\n")
cat("REFIT COMPLETE\n")
cat(paste(rep("=", 50), collapse = ""), "\n")

# Count successful fits
successful_fits <- list.dirs("outputs/traits", recursive = TRUE, full.names = FALSE)
successful_fits <- successful_fits[successful_fits != ""]
cat("Successful fits:", length(successful_fits), "\n")

# Show output structure
cat("\nOutput structure:\n")
system("find outputs/traits -type f -name '*.png' | head -10")

cat("\nCheck outputs/logs/ for detailed logs\n")
