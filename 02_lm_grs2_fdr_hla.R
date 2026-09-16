# SDS 20240709, updated 20260428

# Script to run linear model for analysis of genetic risk in T1D moms:
# GRS2 ~ fdr + covs

# Updated 20241024 to add extraction of confidence intervals to model 1 only to
# allow for making a forest plot for IDS poster.
# Updated 20251014 to only run one linear model, but script history contains
# more complex models with HLAGRP and interaction of HLAGRP and FDR.
# Updated 20260212 & 20260407.

# Setup ------------------------------------------------------------------------

library(here)
devtools::load_all()
source(here("config.R"))

# TO NOTE: change this step's specific settings here, all other settings in
# config.R
outcome <- "GRS2x"  # can be GRS2x or Non_HLA
covs = c("fdr_4level", "PC1", "PC2", "sex", "cc")

# Define function --------------------------------------------------------------

run_model <- function(config) {
  # Define output path
  dr_suffix <- ifelse(config$dr_filt == "yes", "_dr_filt", "")
  out_path <- file.path(
    study_specs$result_out_dir,
    paste0(
      "lm_", outcome,
      "_fdr_", tolower(config$fdr_ref), "_ref",
      dr_suffix, "_", study_specs$study, "_",
      config$subset, ".csv"))

  # Get input paths
  # TODO: bring back option to run in cases & controls in same script?
  kinship <- study_specs$kinship_path
  # Okay to use IA or T1D pheno since both have all participants
  pheno = paste0(
    study_specs$intermed_out_dir, paste0("/pheno_", config$subset, "_ia", dr_suffix, ".rds"))

  # Run model
  message(paste0("Running model for ", config$subset, " & outcome = ", outcome, "."))
  result <- mod(pheno, kinship, outcome, covs, config$fdr_var, config$fdr_ref)

  # Save results
  write_csv(result, file = out_path)
}

# Run function -----------------------------------------------------------------

run_model(config)
