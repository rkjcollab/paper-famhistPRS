#' @export
run_model <- function(
    engine, model, pheno, event, time, covs, fdr_var, fdr_ref,
    wald_test = NULL, kinship = NULL, tt_spec = NULL) {

  message(paste0("Running ", engine, " model for ", model, "."))

  if (engine == "coxph") {
    mod_coxph(pheno, event, time, covs, fdr_var, fdr_ref, wald_test)

  } else if (engine == "coxph_tt") {
    mod_coxph_tt(pheno, event, time, covs, fdr_var, fdr_ref, tt_spec, wald_test)

  } else if (engine == "coxme") {
    mod_coxme(pheno, kinship, event, time, covs, fdr_var, fdr_ref, wald_test)

  } else {
    stop("Unknown engine")
  }
}