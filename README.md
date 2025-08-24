Repository Cleanup & bayesTPC Updates (June 2025)

June 12–14, 2025: bayesTPC Function Updates
	•	Updated to the development version of bayesTPC from GitHub (johnwilliamsmithjr/bayesTPC).
	•	Renamed data loading function:
	•	get_datasets() → get_VB_datasets()
	•	Removed deprecated function bayesTPC_summary(); replaced with standard S3 method: summary().
	•	Posterior prediction functions (plot_prediction, posterior_predictive) now return custom objects (e.g., "btpc_prediction"), requiring manual extraction (e.g., $summary).
	•	plot_prediction() returns NULL unless temp_range is explicitly provided.

⸻

June 29, 2025: Model Fitting Error & NIMBLE Compatibility
	•	Encountered error with b_TPC() due to missing internal getNimbleOption() function in newer versions of nimble.
	•	Fixed by downgrading nimble to version 0.13.2:

remotes::install_version("nimble", version = "0.13.2")



⸻

June 30, 2025: GitHub Push Failed Due to Large Files
	•	Updated .gitignore to exclude cache, site output, and large intermediate files:

*_cache/
_site/
*.rdb
*.rdx
*.RData
*.tar.gz


	•	Removed large files from Git history using git filter-repo:

git filter-repo --force --path-glob '*.rdb' --invert-paths


	•	Reconnected origin remote and force-pushed cleaned history:

git remote add origin https://github.com/VectorByteOrg/bayesTPC-tutorials.git
git push --set-upstream origin yusi --force

⸻

August 19, 2025: Site Reorganization
	•	Merged Getting Started and Materials pages to eliminate redundant content
	•	Updated navbar: "Materials" → "Tutorials" 
	•	Added subsection anchors to EEID tutorials for direct linking
	•	Integrated Quick Start guide with working code example
	•	Added tutorial descriptions with time estimates and learning paths

⸻

December 2025: Complete BTV Refitting Implementation
	•	Successfully implemented comprehensive BTV refitting pipeline for El Moustaid et al. (2021)
	•	Built prior translation system to bridge paper notation and bayesTPC expectations
	•	Created robust data processing pipeline for Appendix A.6 data
	•	Implemented MCMC fitting with custom priors from Appendix A.2
	•	Generated thermal performance curves for 8 trait-species combinations
	•	Reproduced Supplement Fig. A.2 uncertainty attribution analysis
	•	Created comprehensive visualization gallery with HTML interface
	•	Integrated results into main Bluetongue Refit page with blog-style narrative
	•	Added quality control and convergence diagnostics
	•	Documented technical challenges and solutions

Key Results:
	•	Successfully fitted p, b, ν, ρ traits for C. sonorensis and C. variipennis
	•	Development rate (ν) contributes most uncertainty to derived quantities
	•	Vector competence (b) contributes least uncertainty
	•	Optimal temperature ranges: 7-30°C for survival, 35-40°C for competence
	•	All fits achieved R-hat < 1.1 and ESS > 1000

Files Added:
	•	Complete R pipeline (priors_translate.R, fit_trait.R, etc.)
	•	Data files and configuration (priors_tableA2.yaml)
	•	Output files (fits, plots, summaries)
	•	Visualization gallery (HTML interface)
	•	Updated refit_guide.qmd with comprehensive results