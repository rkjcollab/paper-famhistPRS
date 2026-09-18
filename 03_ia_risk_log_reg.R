# SDS 20241023

# Hypothesis is that if maternal protection is not explained by survival bias,
# adding GRS2 to the model will NOT cause maternal protective effect to go away.

# Running IA/T1D/T1Dprog models for TEDDY (survival):
# Base: IA/T1D ~ FDR + PCs + clinical center (TEDDY) + 1|kinship
# Base + GRS2x as covariate

# Setup ------------------------------------------------------------------------

library(here)
library(tidyverse)
library(coxme)
library(survival)
library(splines)
library(car)
devtools::load_all()
source(here("config.R"))

# TO NOTE: change this step's specific settings here. Also applies the dr_filt,
# fdr_var, and fdr_ref settings from config.R.
# TO NOTE: all engines account for relatedness in the model code, either with
# cluster(FID) (coxph, coxph_tt) or with kinship as a random effect (coxme).

# List of covariates to run models with and without
terms <- c("GRS2x") # c("Non_HLA", "GRS2x", "dr34")

# Each list in run_specs will generate one output file with:
# all models run by the specified engine
# all model X subset combinations
# for each model X subset combination, all intxn/wald_test pairs (except for
# sex intxn, which is blocked for when subset = female or male)
run_specs <- list(
  # Primary IA models with time-varying FDR effect
  list(
    # Single engine to use for all models in this set
    engine = "coxph_tt",
    # Every model X subset combination is run
    models = c("IA"),
    subsets = c("all", "female", "male"),
    # For each model X subset combination, run each intxn/wald_test pair
    intxn_wald_pairs = list(
      list(intxn = NA, wald_test = "fdr_4level"),
      list(intxn = "fdr_4level*sex", wald_test = "sex:fdr_4level")
    )
  ),
  # T1D & progression models
  list(
    engine = "coxph",
    models = c("T1D", "T1D_prog"),
    subsets = c("all", "female", "male"),
    intxn_wald_pairs = list(
      list(intxn = NA, wald_test = "fdr_4level"),
      list(intxn = "fdr_4level*sex", wald_test = "sex:fdr_4level")
    )
  ),
  # Plain (non-time-varying) IA models, for the PH-violation check only
  list(
    engine = "coxph",
    models = c("IA"),
    subsets = "all",
    intxn_wald_pairs = list(
      list(intxn = NA, wald_test = "fdr_4level")
    )
  )
)

# Required if any spec uses engine = coxph_tt; can only include either "knots"
# or "df"
tt_spec <- list(
  type = "linear", # "log", "linear", "spline"
  var = "fdr_4level",
  knots = NA,
  df = NA
)

# Run models -------------------------------------------------------------------

# Derived label
dr_suffix <- ifelse(config$dr_filt == "yes", "_dr_filt", "")

# Runs all options from one list entry in run_specs: every model X subset
# combination with each intxn_wald_pairs pair
run_one_spec <- function(engine, models, subsets, intxn_wald_pairs) {
  engine_label <- if (engine == "coxph_tt") {
    ifelse(
      tt_spec$type == "spline",
      paste(engine, tt_spec$type,
        ifelse(!is.na(tt_spec$knots),
          paste0(tt_spec$knots, "knots"),
          paste0(tt_spec$df, "df")
        ),
        sep = "_"
      ),
      paste(engine, tt_spec$type, sep = "_")
    )
  } else {
    engine
  }

  # Run once for each listed subset group
  for (subset in subsets) {
    pheno_surv_path <- list(
      IA = paste0(
        study_specs$intermed_out_dir, "/pheno_", subset, "_ia", dr_suffix, ".rds"
      ),
      T1D = paste0(
        study_specs$intermed_out_dir, "/pheno_", subset, "_t1d", dr_suffix, ".rds"
      ),
      T1D_prog = paste0(
        study_specs$intermed_out_dir, "/pheno_", subset, "_prog", dr_suffix, ".rds"
      )
    )

    # Run once in each subset for each intxn/wald_test pair
    for (pair in intxn_wald_pairs) {
      intxn_term <- pair$intxn
      wald_test <- pair$wald_test

      # If subset == female/male, skip interaction terms that include sex
      if ((subset == "female" | subset == "male") &&
        !is.na(intxn_term) &&
        grepl("sex", intxn_term, ignore.case = TRUE)) {
        message("Skipping interaction term ", intxn_term, " for subset = ", subset, ".")
        next
      }

      intxn_suffix <- ifelse(is.na(intxn_term), "", "_intxn")

      model_results <- vector("list", length(models))

      for (m in seq_along(models)) {
        model <- models[m]

        message(
          "Running: engine=", engine, " subset=", subset, " model=", model,
          " wald_test=", wald_test, " intxn=", intxn_term, "."
        )

        # Load options that are same for all models
        kinship <- study_specs$kinship_path
        fdr_var <- config$fdr_var
        fdr_ref <- config$fdr_ref

        # Get survival variables for current model
        surv_def <- study_specs$surv_def[[model]]
        pheno <- pheno_surv_path[[model]]
        event <- surv_def$event
        time <- surv_def$time
        covs <- study_specs$surv_covs[[model]]

        # If subset == female/male, remove sex from covs
        if (subset == "female" | subset == "male") {
          covs <- grep("sex", covs, value = T, invert = T)
        }

        if (!is.na(intxn_term)) {
          # Replace individual terms with interaction term
          intxn_terms <- paste(str_split(intxn_term, "\\*")[[1]], collapse = "|")
          covs <- grep(intxn_terms, covs, value = T, invert = T)
          covs <- c(covs, intxn_term)
        }

        result_base <- run_model(
          engine, model, pheno, event, time, covs, fdr_var, fdr_ref,
          wald_test = wald_test, kinship = kinship, tt_spec = tt_spec
        )
        base_row <- result_base$result_df
        base_row$model <- model
        base_row$term <- NA

        # Save model object
        base_mod <- result_base$model_obj
        covs_mod <- attr(terms(base_mod), "term.labels")
        out_path_base_mod <- paste0(
          study_specs$result_out_dir, "/model_objects/",
          engine_label,
          "_", model,
          "_fdr_", tolower(fdr_ref), "_ref_",
          paste(covs_mod, collapse = "-"),
          dr_suffix, "_", study_specs$study, "_",
          subset, ".rds"
        )
        saveRDS(base_mod, file = out_path_base_mod)

        term_result_rows <- vector("list", length(terms))

        for (t in seq_along(terms)) {
          term <- terms[t]
          covs_term <- c(covs, term)

          result_term_tmp <- run_model(
            engine, model, pheno, event, time, covs_term, fdr_var, fdr_ref,
            wald_test = wald_test, kinship = kinship, tt_spec = tt_spec
          )

          result_term <- result_term_tmp$result_df
          result_term$model <- model
          result_term$term <- term

          # Save term model
          term_mod <- result_term_tmp$model_obj
          covs_mod <- attr(terms(term_mod), "term.labels")
          out_path_term_mod <- paste0(
            study_specs$result_out_dir, "/model_objects/",
            engine_label,
            "_", model,
            "_fdr_", tolower(fdr_ref), "_ref_",
            paste(covs_mod, collapse = "-"),
            dr_suffix, "_", study_specs$study, "_",
            subset, ".rds"
          )
          saveRDS(term_mod, file = out_path_term_mod)

          term_result_rows[[t]] <- result_term
        }

        model_results[[m]] <- bind_rows(base_row, bind_rows(term_result_rows))
      }

      # Reformat results
      results <- bind_rows(model_results) %>%
        dplyr::relocate(model, term, .after = outcome)

      # Save results
      out_path <- paste0(
        study_specs$result_out_dir, "/",
        engine_label,
        "_", paste(models, collapse = "_"),
        "_fdr_", tolower(config$fdr_ref), "_ref_",
        paste0(terms, collapse = "_"),
        dr_suffix, "_", study_specs$study, "_",
        subset, intxn_suffix, ".csv"
      )
      write_csv(results, file = out_path)
    }
  }
}

for (spec in run_specs) {
  run_one_spec(spec$engine, spec$models, spec$subsets, spec$intxn_wald_pairs)
}
