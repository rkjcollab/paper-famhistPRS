# Test
# model <- models[1]
# specs <- study_specs[[study]]
# pheno <- specs$pheno_surv_path[[model]]
# surv_def <- specs$surv_def[[paste0(model, "_surv_def")]]
# outcome <- surv_def$event
# time  <- surv_def$time
# covs <- specs$surv_covs[[model]]
# fdr_var <- config$fdr_var
# fdr_ref <- config$fdr_ref
# tt_spec <- tt_spec
# type <- tt_spec$type
# n_knots <- tt_spec$knots
# df <- tt_spec$df

mod_coxph <- function(
    study, pheno, outcome, time, covs, fdr_var, fdr_ref) {

  ### Prep
  # Read in pheno
  pheno <- readRDS(pheno)
  ids <- pheno$ID  # the ID column used in (1|ID)

  # Set FDR reference level for given FDR (3 or 4-level)
  pheno[[fdr_var]] <- relevel(as.factor(pheno[[fdr_var]]), ref = fdr_ref)

  ### Model:
  form <- as.formula(
    paste(paste0(
      "Surv(", time, ",", outcome, ") ~ "),
      paste(covs, collapse = "+")))

  mod <- coxph(form, data = pheno)

  # Get overall model values
  results <- data.frame(study = study)
  results$outcome <- outcome
  results$nobs = mod$n
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
  p <- summary(mod)$coefficients[, "Pr(>|z|)"]
  results_p <- as.data.frame(t(p))
  colnames(results_p) <- paste0("p.value_", names(mod$coefficients))

  # Get CIs
  ci_low <- beta - 1.96 * se
  ci_high <- beta + 1.96 * se
  results_ci_low <- as.data.frame(t(ci_low))
  colnames(results_ci_low) <- paste0("conf.low_", names(mod$coefficients))
  results_ci_high <- as.data.frame(t(ci_high))
  colnames(results_ci_high) <- paste0("conf.high_", names(mod$coefficients))

  # Combine
  results_all_tmp <- cbind(
    results_est, results_se, results_z, results_p,
    results_ci_low, results_ci_high)
  colnames(results_all_tmp) <- colnames(results_all_tmp) %>%
    gsub(fdr_var, "fdr_", .) %>%
    gsub("Mom", "mom", .) %>%
    gsub("None", "none", .) %>%
    gsub("Sib", "sib", .) %>%
    gsub("Dad", "dad", .)
  
  results_all <- as.data.frame(cbind(results, results_all_tmp))
  
  # Return result dataframe and model object
  return(list(
    result_df = results_all,
    model_obj = mod))
}


build_tt <- function(type, pheno, time, outcome, n_knots = NA, df = NA) {

  knots <- NULL
  boundary <- NULL

  if (type == "spline" && !is.na(n_knots) && !is.na(df)) {
    exit("For spline, can only specify either tt_spec$knots or tt_spec$df.")
  }

  if (type == "spline" && !is.na(n_knots)) {
    # Set spline knots based on quantiles in cases
    pheno_case <- pheno[pheno[[outcome]] == 1, ]

    # Set quantile options based on:
      # PMID 25883970 for knots = 3 & 5
      # knowledge of data/prior literature for knots = 1
    if (n_knots == 1) {
      # Single knot at age 5
      knots <- as.numeric(5)
    } else if (n_knots == 3) {
      probs <- c(0.05, 0.5, 0.95)
      knots <- as.numeric(quantile(pheno_case[[time]], probs = probs))
    } else if (n_knots == 5) {
      probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
      knots <- as.numeric(quantile(pheno_case[[time]], probs = probs))
    } else {
      stop("Current model only supports spline knots = 1, 3, or 5.")
    }
    boundary <- range(pheno[[time]])
  }

  if (type == "spline" && !is.na(df)) {
    t_ref <- pheno[[time]]
    s_ref <- splines::ns(t_ref, df = df)
    
    knots <- attr(s_ref, "knots")
    boundary <- attr(s_ref, "Boundary.knots")
  }

  tt_fun <- function(x, t, ...) {

    # Convert factor to dummy matrix
    if (is.factor(x)) {
      mm <- model.matrix(~ x)[, -1, drop = FALSE]
      colnames(mm) <- gsub("^x", "", colnames(mm))
    } else {
      mm <- cbind(x = x)
    }

    if (type == "log") {
      return(mm * log(t))
    } else if (type == "linear") {
      return(mm * t)
    } else if (type == "spline") {
      s <- splines::ns(t, knots = knots, Boundary.knots = boundary)
      out <- do.call(cbind, lapply(seq_len(ncol(mm)), function(j) {
        mm[, j] * s
      }))
      colnames(out) <- as.vector(
        outer(colnames(mm), colnames(s), paste0)
      )
      return(out)
    }
    stop("Invalid tt type.")
  }
  list(fun = tt_fun, knots = knots, boundary = boundary)
}

mod_coxph_tt <- function(
    study, pheno, outcome, time, covs, fdr_var, fdr_ref, tt_spec) {
  
  ### Prep
  # Read in pheno
  pheno <- readRDS(pheno)
  ids <- pheno$ID  # the ID column used in (1|ID)
  
  # Set FDR reference level for given FDR (3 or 4-level)
  pheno[[fdr_var]] <- relevel(as.factor(pheno[[fdr_var]]), ref = fdr_ref)
  
  ### Model
  form <- as.formula(
    paste(paste0(
      "Surv(", time, ",", outcome, ") ~ "),
      paste(covs, collapse = "+"),
      "+ tt(", tt_spec$var, ")"))
  
  # Get time-transform object
  knots <- if (is.null(tt_spec$knots)) NA else tt_spec$knots
  df <- if (is.null(tt_spec$df)) NA else tt_spec$df
  if (!is.na(knots) && !is.na(df)) {
    exit("Can only specify either tt_spec$knots or tt_spec$df.")
  }
  tt_obj <- build_tt(
    type = tt_spec$type,
    pheno = pheno,
    time = time,
    outcome = outcome,
    n_knots = knots,
    df = df)
  if (tt_spec$type == "spline") {
    if (!is.na(knots)) {
      message("Using spline with ", knots, " knots (quantile-based).")
    } else if (!is.na(df)) {
      message("Using spline with df = ", df, ".")
    }
  } else {
    message("Using ", tt_spec$type, " time transform.")
  }
  
  # Run model
  mod <- coxph(
    form,
    data = pheno,
    tt = tt_obj$fun)
  mod$tt_meta <- list(
    type = tt_spec$type,
    var = tt_spec$var,
    knots = tt_obj$knots,
    boundary = tt_obj$boundary)
  
  # Get overall model values
  results <- data.frame(study = study)
  results$outcome <- outcome
  results$nobs = mod$n
  results$form = paste0(
    as.character(form)[2], as.character(form)[1], as.character(form)[3])
  results$which_fdr <- fdr_var
  results$which_tt <- tt_spec$type
  results$which_tt_var <- tt_spec$var
  
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
  p <- summary(mod)$coefficients[, "Pr(>|z|)"]
  results_p <- as.data.frame(t(p))
  colnames(results_p) <- paste0("p.value_", names(mod$coefficients))
  
  # Get CIs
  ci_low <- beta - 1.96 * se
  ci_high <- beta + 1.96 * se
  results_ci_low <- as.data.frame(t(ci_low))
  colnames(results_ci_low) <- paste0("conf.low_", names(mod$coefficients))
  results_ci_high <- as.data.frame(t(ci_high))
  colnames(results_ci_high) <- paste0("conf.high_", names(mod$coefficients))
  
  # Combine
  results_all_tmp <- cbind(
    results_est, results_se, results_z, results_p,
    results_ci_low, results_ci_high)
  colnames(results_all_tmp) <- colnames(results_all_tmp) %>%
    gsub(fdr_var, "fdr_", .) %>%
    gsub("fdr_)", "fdr)", .) %>%
    gsub("Mom", "mom", .) %>%
    gsub("None", "none", .) %>%
    gsub("Sib", "sib", .) %>%
    gsub("Dad", "dad", .)
  
  results_all <- as.data.frame(cbind(results, results_all_tmp))
  
  # Return result dataframe and model object
  return(list(
    result_df = results_all,
    model_obj = mod))
}
