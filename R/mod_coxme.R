
#' @export
mod_coxme <- function(
    study, pheno, kinship, outcome, time, covs, fdr_var, fdr_ref,
    wald_test = NULL) {
  
  ### Prep
  # Read in pheno
  pheno <- readRDS(pheno)
  ids <- pheno$ID  # the ID column used in (1|ID)
  
  # Read in kinship, subset to analytical data
  kinship <- readRDS(kinship)
  kinship_mat <- as.matrix(kinship)[ids, ids] # reorder rows/cols
  stopifnot(all(rownames(kinship_mat) == ids))
  
  # Set FDR reference level for given FDR (3 or 4-level)
  pheno[[fdr_var]] <- relevel(as.factor(pheno[[fdr_var]]), ref = fdr_ref)
  
  ### Model:
  # Formula automatically includes (1|ID) for kinship
  form <- as.formula(
    paste(paste0(
      "Surv(", time, ",", outcome, ") ~ (1|ID) + "),
      paste(covs, collapse = "+")))
  
  mod <- coxme(form, data = pheno, varlist = kinship_mat)
  
  # Wald test
  if (!is.null(wald_test)) {
    beta <- coef(mod)
    v <- vcov(mod)
    
    if (grepl(":", wald_test)) {
      # Interaction term
      parts <- strsplit(wald_test, ":")[[1]]
      term_pattern <- paste0(
        "(", parts[1], ".*:", parts[2], ")|(",
        parts[2], ".*:", parts[1], ")")
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
  results$outcome <- outcome
  results$nobs = mod$n[2]  # get total N, not just cases
  results$form = paste0(
    as.character(form)[2], as.character(form)[1], as.character(form)[3])
  results$which_fdr <- fdr_var
  
  # Extract all fixed effects
  beta <- as.numeric(mod$coefficients)
  results_est <- as.data.frame(t(beta))
  colnames(results_est) <- paste0("estimate_", names(mod$coefficients))
  
  # Get standard errors for fixed effects
  se <- summary(mod)$coefficients[, "se(coef)"]
  results_se <- as.data.frame(t(se))
  colnames(results_se) <- paste0("std.error_", names(mod$coefficients))
  
  # Get z-statistics and p-values
  z <- summary(mod)$coefficients[, "z"]
  results_z <- as.data.frame(t(z))
  colnames(results_z) <- paste0("z.value_", names(mod$coefficients))
  p <- summary(mod)$coefficients[, "p"]
  results_p <- as.data.frame(t(p))
  colnames(results_p) <- paste0("p.value_", names(mod$coefficients))
  
  # Get CIs
  ci_low <- beta - 1.96 * se
  ci_high <- beta + 1.96 * se
  results_ci_low <- as.data.frame(t(ci_low))
  colnames(results_ci_low) <- paste0("conf.low_", names(mod$coefficients))
  results_ci_high <- as.data.frame(t(ci_high))
  colnames(results_ci_high) <- paste0("conf.high_", names(mod$coefficients))
  
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
  
  # Combine
  results_all_tmp <- cbind(
    results_est, results_se, results_z, results_p,
    results_ci_low, results_ci_high, results_w)
  colnames(results_all_tmp) <- colnames(results_all_tmp) %>%
    gsub(fdr_var, "fdr_", .) %>%
    gsub("__+", "_", .) %>%
    gsub("Mom", "mom", .) %>%
    gsub("None", "none", .) %>%
    gsub("Sib", "sib", .) %>%
    gsub("Dad", "dad", .) %>%
    gsub("_:", ":", .)  # in case fdr_var is also wald_test
  
  results_all <- as.data.frame(cbind(results, results_all_tmp))
  
  return(list(
    result_df = results_all,
    model_obj = mod))
}
