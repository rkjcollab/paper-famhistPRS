# SDS 20260511
# Originally based on code from Kirk Hohsfield.

# Use for 'log' or 'identity' (linear) time-transform, must match the transform
# used when model was run
# mod = m_ia_lin_g
# var = "fdr_4level"
# levels = c("None", "Mom", "Sib")
# ref_label = "Dad"
# time_transform = identity
# lvl <- "Mom"
# times = seq(0, 16, 1)
get_tt_hr <- function(
    mod,
    var,
    levels,
    ref_label = "Reference",
    times = seq(0, 16, 0.5),
    time_transform) {
  
  coefs <- coef(mod)
  vcov_mat <- vcov(mod)
  
  # Loop over levels
  results <- lapply(levels, function(lvl) {
    
    name_var_main <- paste0(var, lvl)
    name_var_tt   <- paste0("tt(", var, ")", lvl)
    
    if (!(name_var_main %in% names(coefs)) || !(name_var_tt %in% names(coefs))) {
      stop(paste("Missing coefficients for level:", lvl))
    }
    
    df <- calc_hr_ci(
      beta_main = coefs[name_var_main],
      beta_tt   = coefs[name_var_tt],
      var_main  = vcov_mat[name_var_main, name_var_main],
      var_tt    = vcov_mat[name_var_tt, name_var_tt],
      cov_main_tt = vcov_mat[name_var_main, name_var_tt],
      times = times,
      time_transform = time_transform
    )
    
    df$cohort <- paste(lvl, "vs", ref_label)
    df
  })
  
  do.call(rbind, results)
}


# This function is for spline models, same knots/boundaries should be used as
# when model was run
# mod = m_ia_spline_3
# var = mod$tt_meta$var
# levels = c("None", "Mom", "Sib")
# ref_label = "Dad"
# times = seq(0, 16, 0.5)
# knots = mod$tt_meta$knots
# boundary = mod$tt_meta$boundary
get_tt_hr_spline <- function(
    mod,
    var,
    levels,
    ref_label,
    times = seq(0, 16, 0.5),
    knots,
    boundary) {
  coefs <- coef(mod)
  vcov_mat <- vcov(mod)
  
  # Build spline basis (must match model!)
  S <- splines::ns(
    times,
    knots = knots,
    Boundary.knots = boundary)
  
  results <- lapply(levels, function(lvl) {
    
    var_main <- paste0(var, lvl)
    
    # Find all spline coefficients for this level
    pattern <- paste0("tt\\(", var, "\\)", lvl)
    coef_names <- names(coefs)
    tt_names <- coef_names[grepl(pattern, coef_names)]
    spline_names <- colnames(S)
    expected_names <- paste0("tt(", var, ")", lvl, spline_names)
    beta_vec <- coefs[expected_names]
    
    if (any(is.na(beta_vec))) {
      stop(paste("Mismatch between spline basis and coefficients for", lvl))
    }
    beta_main <- coefs[var_main]
    
    # Compute log HR
    log_hr <- as.numeric(beta_main + S %*% beta_vec)
    
    # Variance
    var_main_val <- vcov_mat[var_main, var_main]
    vcov_tt      <- vcov_mat[expected_names, expected_names]
    cov_main_tt  <- vcov_mat[var_main, expected_names]
    
    # Var = Var(main) + S Σ S^T + 2 * cov
    se_log_hr <- sapply(1:nrow(S), function(i) {
      s <- S[i, ]
      
      var_tt_part <- t(s) %*% vcov_tt %*% s
      cov_part <- 2 * sum(s * cov_main_tt)
      
      sqrt(var_main_val + var_tt_part + cov_part)
    })
    
    data.frame(
      time = times,
      hr = exp(log_hr),
      lower = exp(log_hr - 1.96 * se_log_hr),
      upper = exp(log_hr + 1.96 * se_log_hr),
      cohort = paste(lvl, "vs", ref_label)
    )
  })
  
  do.call(rbind, results)
}

# Function to calculate HR and CI at each time point
calc_hr_ci <- function(
    beta_main,
    beta_tt,
    var_main,
    var_tt,
    cov_main_tt,
    times,
    time_transform) {
  
  g_t <- time_transform(times)
  
  log_hr <- beta_main + beta_tt * g_t
  se_log_hr <- sqrt(var_main + (g_t^2) * var_tt + 2 * g_t * cov_main_tt)
  
  data.frame(
    time = times,
    hr = exp(log_hr),
    lower = exp(log_hr - 1.96 * se_log_hr),
    upper = exp(log_hr + 1.96 * se_log_hr))
}


plot_tt_hr <- function(
    data,
    var,
    ref_label,
    times = seq(0, 16, 0.5)) {
  
  p <- ggplot(data, aes(x = time, y = hr, color = cohort, fill = cohort)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
    scale_y_log10() +
    scale_x_continuous(breaks = times) +
    labs(
      x = "Follow-up Time",
      y = "Hazard Ratio",
      title = paste("Time-Varying HR for", var),
      subtitle = paste("Reference:", ref_label),
      color = "Comparison",
      fill = "Comparison"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank())
  
  return(p)
  
}


