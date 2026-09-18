# SDS 20240709

# Script to run linear model for analysis of genetic risk in T1D moms:
# GRS2 ~ fdr + covs

# Setup ------------------------------------------------------------------------

library(here)
library(tidyverse)
library(coxme)
devtools::load_all()
source(here("config.R"))

# TO NOTE: change this step's specific settings here. Also applies the dr_filt,
# fdr_var, and fdr_ref settings from config.R.
outcome <- "GRS2x" # "GRS2x" or "Non_HLA"
covs <- c("fdr_4level", "PC1", "PC2", "sex", "cc")
subset <- "all" # all, ctrls, female, or male
wald_test <- "fdr_4level" # "fdr_4level" or NULL

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
      subset, ".csv"
    )
  )

  # Get input paths
  kinship <- study_specs$kinship_path
  # Okay to use IA or T1D pheno since both have all participants
  pheno <- paste0(
    study_specs$intermed_out_dir, paste0("/pheno_", subset, "_ia", dr_suffix, ".rds")
  )

  # Run model
  message(paste0("Running model for ", subset, " & outcome = ", outcome, "."))
  result <- mod(
    pheno, kinship, outcome, covs, config$fdr_var, config$fdr_ref,
    wald_test = wald_test
  )

  # Save results
  write_csv(result, file = out_path)
}

# Run function -----------------------------------------------------------------

run_model(config)
