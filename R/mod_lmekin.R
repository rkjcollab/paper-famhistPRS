library(coxme)

# fdr_var <- config$fdr_var
# fdr_ref <- config$fdr_ref

# Note: pheno files must be set to have id = "ID"
#' @export
mod <- function(
    study, pheno, kinship, outcome, covs, fdr_var, fdr_ref, wald_test = NULL) {
  
  # Read in pheno
  pheno <- readRDS(pheno)
  ids <- pheno$ID  # the ID column used in (1|ID)
  
  # Read in kinship, subset to analytical data
  kinship <- readRDS(kinship)
  kinship_mat <- as.matrix(kinship)[ids, ids] # reorder rows/cols
  stopifnot(all(rownames(kinship_mat) == ids))
  
  # Set FDR reference level for given FDR (3 or 4-level)
  pheno[[fdr_var]] <- relevel(as.factor(pheno[[fdr_var]]), ref = fdr_ref)
  
  ### Model: GRS2 ~ fdr + covs
  form = as.formula(paste0(outcome, "~", paste(covs, collapse = "+"), "+ (1|ID)"))
  
  # Run linear mixed model
  mod = lmekin(form, data = pheno, varlist = kinship_mat)
  
  # Wald test
  if (!is.null(wald_test)) {
    beta <- fixef(mod)
    v <- vcov(mod)
    
    if (grepl(":", wald_test)) {
      stop("Cannot compute a Wald test on an interaction term.")
    } else {
      # Term
      term_pattern <- paste0("^", wald_test)
    }
    term_idx <- which(grepl(term_pattern, names(beta)))
    if (length(term_idx) == 0) {
      stop("Wald test term not found in coefficients: ", wald_test)
    }
    
    beta_sub <- beta[term_idx]
    v_sub <- v[term_idx, term_idx, drop = FALSE]
    
    # Wald statistic
    chisq <- as.numeric(t(beta_sub) %*% solve(v_sub) %*% beta_sub)
    df <- length(beta_sub)
    pval <- pchisq(chisq, df = df, lower.tail = FALSE)
    
    wald <- data.frame(
      chisq = chisq,
      df = df,
      p = pval)
  }
  
  # Get overall model values
  results <- data.frame(study = study)
  results$nobs = mod$n
  results$form = paste0(
    as.character(form)[2], as.character(form)[1], as.character(form)[3])
  results$which_fdr <- fdr_var
  
  # Based on Erika's lmekin code
  # Extract all fixed effects
  beta <- as.numeric(mod$coefficients$fixed)
  results_est <- as.data.frame(t(beta))
  colnames(results_est) <- paste0("estimate_", names(mod$coefficients$fixed))
  
  # Get standard errors for fixed effects
  nvar <- length(beta)
  nfrail <- nrow(mod$var) - nvar  # number random effects
  se <- as.numeric(sqrt(diag(mod$var)[nfrail + 1:nvar]))
  results_se <- as.data.frame(t(se))
  colnames(results_se) <- paste0("std.error_", names(mod$coefficients$fixed))
  
  # Get z-statistics and p-values
  z <- as.numeric(round(beta/se, 5))
  results_z <- as.data.frame(t(z))
  colnames(results_z) <- paste0("z.value_", names(mod$coefficients$fixed))
  p <- as.numeric(signif(1-pchisq((beta/se)^2,1), 5))
  results_p <- as.data.frame(t(p))
  colnames(results_p) <- paste0("p.value_", names(mod$coefficients$fixed))
  
  # Get CIs
  ci_low <- beta - 1.96 * se
  ci_high <- beta + 1.96 * se
  results_ci_low <- as.data.frame(t(ci_low))
  colnames(results_ci_low) <- paste0("conf.low_", names(mod$coefficients$fixed))
  results_ci_high <- as.data.frame(t(ci_high))
  colnames(results_ci_high) <- paste0("conf.high_", names(mod$coefficients$fixed))
  
  # Extract Wald test on given variable
  if (!is.null(wald_test)) {
    results_w <- data.frame(
      wald$p,
      wald$chisq,
      wald$df)
    colnames(results_w) <- c(
      paste0(wald_test, "_global_p"),
      paste0(wald_test, "_global_chisq"),
      paste0(wald_test, "_global_df"))
  } else {
    results_w <- data.frame()
  }
  
  # Combine, standardize names
  results_all_tmp <- cbind(
    results_est, results_se, results_z, results_p,
    results_ci_low, results_ci_high, results_w)
  
  colnames(results_all_tmp) <- colnames(results_all_tmp) %>%
    gsub(fdr_var, "fdr_", .) %>%
    gsub("__+", "_", .) %>%
    gsub("Mom", "mom", .) %>%
    gsub("None", "none", .) %>%
    gsub("Sib", "sib", .) %>%
    gsub("Dad", "dad", .)
  results_all <- cbind(results, results_all_tmp)
  
  return(as.data.frame(results_all))
}