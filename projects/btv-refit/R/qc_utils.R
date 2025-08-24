# qc_utils.R
# Quality control and convergence diagnostics

# R-hat helper
get_rhats <- function(fit) {
  tryCatch({
    summary(fit)$rhat
  }, error = function(e) NULL)
}

# Check convergence
check_convergence <- function(fit) {
  rhats <- get_rhats(fit)
  if (!is.null(rhats) && any(rhats > 1.1, na.rm = TRUE)) {
    warning(sprintf("Convergence issue: %s", 
                   paste(names(rhats)[rhats > 1.1], collapse = ", ")))
  }
  rhats
}

# Generate quality control summary
generate_qc_summary <- function(res, trait, species) {
  # Convergence diagnostics
  rhats <- check_convergence(res$fit)
  
  # Data summary
  data_summary <- list(
    n_obs = nrow(res$data),
    T_range = range(res$data$T),
    y_range = range(res$data$y),
    trait = trait,
    species = species
  )
  
  # Parameter constraints (if applicable)
  constraints <- list()
  if (res$model == "briere") {
    # Check T_min < T_max
    T_min_vals <- res$draws$T_min
    T_max_vals <- res$draws$T_max
    constraints$T_min_lt_T_max <- all(T_min_vals < T_max_vals)
    
    # Check q > 0
    q_vals <- res$draws$q
    constraints$q_positive <- all(q_vals > 0)
  }
  
  # Package versions
  session_info <- list(
    R_version = R.version.string,
    bayesTPC_version = as.character(packageVersion("bayesTPC")),
    nimble_version = as.character(packageVersion("nimble"))
  )
  
  list(
    convergence = rhats,
    data = data_summary,
    constraints = constraints,
    session_info = session_info,
    timestamp = Sys.time()
  )
}

# Save QC summary to JSON
save_qc_summary <- function(qc_summary, outdir) {
  jsonlite::write_json(qc_summary, file.path(outdir, "qc.json"), 
                       pretty = TRUE, auto_unbox = TRUE)
}

# Save session info
save_session_info <- function(outdir) {
  session_info <- list(
    R_version = R.version.string,
    bayesTPC_version = packageVersion("bayesTPC"),
    nimble_version = packageVersion("nimble"),
    timestamp = Sys.time()
  )
  
  writeLines(
    c("Session Information:",
      paste("R version:", session_info$R_version),
      paste("bayesTPC version:", session_info$bayesTPC_version),
      paste("nimble version:", session_info$nimble_version),
      paste("Timestamp:", session_info$timestamp)
    ),
    file.path(outdir, "session_info.txt")
  )
}
