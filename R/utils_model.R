#' @export
run_model <- function(
    engine, study, model, pheno, event, time, covs, fdr_var, fdr_ref,
    wald_test = NULL, kinship = NULL, tt_spec = NULL) {
  
  message(paste0("Running ", engine, " model for ", study, " & ", model, "."))
  
  if (engine == "coxph") {
    mod_coxph(study, pheno, event, time, covs, fdr_var, fdr_ref, wald_test)
    
  } else if (engine == "coxph_tt") {
    mod_coxph_tt(study, pheno, event, time, covs, fdr_var, fdr_ref, tt_spec, wald_test)
    
  } else if (engine == "coxme") {
    mod_coxme(study, pheno, kinship, event, time, covs, fdr_var, fdr_ref, wald_test)
    
  } else {
    stop("Unknown engine")
  }
}