# SDS 20241023, updated 20260429

# Hypothesis is that if maternal protection is not explained by survival bias,
# adding GRS2 to the model will NOT cause maternal protective effect to go away.

# Running IA/T1D/T1Dprog models for TEDDY (survival):
  # Base: IA/T1D ~ FDR + PCs + clinical center (TEDDY) + 1|kinship
  # Base + GRS2x as covariate

#TODO: if keep interaction p-value in results, need to update file name to 
# write out differently

# Setup ------------------------------------------------------------------------

library(here)
library(tidyverse)
library(GMMAT)
library(coxme)
library(survival)
library(splines)
library(car)

# devtools::install()
# library(famhistPRS)
devtools::load_all()
source(here("config.R"))

# TO NOTE: change this step's specific settings here, all other settings in
# config.R and study_specs.R
study <- config$studies[1]
models <- c("IA")  # "IA", "T1D", "T1D_prog"
terms <- c("GRS2x")  # c("Non_HLA", "GRS2x", "dr34")
engine <- "coxph_tt"  # coxme, coxph, or coxph_tt
wald_test <- "fdr_4level"  # optional term to run Wald test one
  # TODO: needs to align with covs

# Required if engine = coxph_tt, can only include either "knots" or "df"
tt_spec <- list(
  type = "linear",  # "log", "linear", "spline"
  var = "fdr_4level",
  knots = NA,
  df = NA
)

# Derived pheno paths, not directly edited
pheno_surv_path <- list(
  IA = paste0(
    study_specs[[study]]$intermed_out_dir, "/pheno_", config$subset, "_ia", config$dr_suffix, ".rds"),
  T1D = paste0(
    study_specs[[study]]$intermed_out_dir, "/pheno_", config$subset, "_t1d", config$dr_suffix, ".rds"),
  T1D_prog = paste0(
    study_specs[[study]]$intermed_out_dir, "/pheno_", config$subset, "_prog", config$dr_suffix, ".rds"))

# Derived labels, not directly edited
dr_suffix <- ifelse(config$dr_filt == "yes", "_dr_filt", "")
engine_label <- if (engine == "coxph_tt") {
  ifelse(
    tt_spec$type == "spline",
    paste(engine, tt_spec$type,
          ifelse(!is.na(tt_spec$knots),
                 paste0(tt_spec$knots, "knots"),
                 paste0(tt_spec$df, "df")),
          sep = "_"),
    paste(engine, tt_spec$type, sep = "_"))
} else {
  engine
}


# Run models -------------------------------------------------------------------

results <- map_dfr(models, function(model) {
    # Get current study specs
    specs <- study_specs[[study]]
    
    # Load options that are same for all models
    kinship <- specs$kinship_path
    fdr_var <- config$fdr_var
    fdr_ref <- config$fdr_ref
    subset <- config$subset
    
    # Get survival variables for current model
    surv_def <- specs$surv_def[[model]]
    pheno <- pheno_surv_path[[model]]
    event <- surv_def$event
    time  <- surv_def$time
    covs <- specs$surv_covs[[model]]
    
    # If subset == female/male, remove sex from covs
    if (subset == "female" | subset == "male") {
      covs <- grep("sex", covs, value = T, invert = T)
    }
    
    result_base <- run_model(
      engine, model, study, pheno, event, time, covs, fdr_var, fdr_ref,
      wald_test = wald_test, kinship = kinship, tt_spec = tt_spec)
    base_row <- result_base$result_df
    base_row$term = NA
    
    # Save model object
    base_mod <- result_base$model_obj
    covs_mod <- attr(terms(base_mod), "term.labels")
    out_path_base_mod <- paste0(
      specs$result_out_dir, "/model_objects/",
      engine_label,
      "_", model,
      "_fdr_", tolower(fdr_ref), "_ref_",
      paste(covs_mod, collapse = "-"),
      dr_suffix, "_",
      study, "_",
      subset, ".rds")
    saveRDS(base_mod, file = out_path_base_mod)
    
    term_rows <- map_dfr(terms, function(term) {
      covs_term <- c(covs, term)
      
      result_term_tmp <- run_model(
        engine, model, study, pheno, event, time, covs_term, fdr_var, fdr_ref,
        wald_test = wald_test, kinship = kinship,tt_spec = tt_spec)
      
      result_term <- result_term_tmp$result_df
      result_term$term <- term
        
      # Save term model
      term_mod <- result_term_tmp$model_obj
      covs_mod <- attr(terms(term_mod), "term.labels")
      out_path_term_mod <- paste0(
        specs$result_out_dir, "/model_objects/",
        engine_label,
        "_", model,
        "_fdr_", tolower(fdr_ref), "_ref_",
        paste(covs_mod, collapse = "-"),
        dr_suffix, "_",
        study, "_",
        subset, ".rds")
      saveRDS(term_mod, file = out_path_term_mod)
      
      result_term
    })
    
    bind_rows(base_row, term_rows)
})

# Reformat results
results <- results %>%
  dplyr::relocate(term, .after = outcome)

# Save results
out_path <- paste0(
  study_specs[[study]]$result_out_dir, "/",
  engine_label,
  "_", paste(models, collapse = "_"),
  "_fdr_", tolower(config$fdr_ref), "_ref_",
  paste0(terms, collapse = "_"),
  dr_suffix, "_",
  study, "_",
  config$subset, ".csv")
write_csv(results, file = out_path)
