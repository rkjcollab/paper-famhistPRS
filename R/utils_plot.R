
# data <- m1_dad
# ref <- "dad"
make_h1_forest_plot <- function(data, ref) {
  
  # Map reference to the other terms
  fdr_all <- c("dad", "mom", "sib", "none")
  fdr <- setdiff(fdr_all, ref)
  
  # Dynamically select columns
  cols_to_select <- c(
    "study",  # "group",
    paste0("estimate_fdr_", fdr),
    paste0("conf.low_fdr_", fdr),
    paste0("conf.high_fdr_", fdr))
  
  data_plot <- data %>%
    select(all_of(cols_to_select)) %>%
    pivot_longer(
      cols = matches("estimate_|conf\\.low_|conf\\.high_"),
      names_to = c(".value", "term"),
      names_pattern = "(estimate|conf\\.low|conf\\.high)_(.*)") %>%
    dplyr::mutate(term = recode(
      term,
      fdr_mom = "Mother",
      fdr_sib = "Sibling",
      fdr_dad = "Father",
      fdr_none = "None")) %>%
    dplyr::mutate(term = factor(term, levels = rev(fdr_order)))
  
  # Build plot
  plot <- ggplot(data_plot, aes(
    x = term, y = estimate, ymin = conf.low, ymax = conf.high, color = term, fill = term)) +
    geom_linerange(position = position_dodge(width = 0.5), show.legend = FALSE, size = 0.8) +
    geom_hline(yintercept = 0, lty = 2) +
    geom_point(size = 2.5, shape = 21, colour = "white", stroke = 0.5,
               position = position_dodge(width = 0.5)) +
    scale_fill_manual(values = c("Father" = color_fdr_dad,
                                 "Mother" = color_fdr_mom,
                                 "Sibling" = color_fdr_sib,
                                 "None" = color_fdr_none)) +
    scale_color_manual(values = c("Father" = color_fdr_dad,
                                  "Mother" = color_fdr_mom,
                                  "Sibling" = color_fdr_sib,
                                  "None" = color_fdr_none)) +
    coord_flip() +
    theme_bw() +
    labs(x = "", y = "Difference in GRS2x", fill = "FDR") +
    theme(legend.position = "none",
          strip.background = element_rect(fill = "white", color = "black"),
          axis.text.y = element_text(angle = 90, hjust = 0.5))
  
  return(plot)
  
}

#TODO: this version from time-varying report, need to reconcile with manuscript
# version!
# data <- m2_dad
# ref <- "dad"
# term <- "GRS2x"
make_h2_forest_plot <- function(data, ref, fdr_order, outcome_order) {
  
  fdr_all <- c("dad", "mom", "sib", "none")
  fdr <- setdiff(fdr_all, ref)
  
  cols_to_select <- c(
    "study", "model", "outcome", "model_group",
    paste0("estimate_fdr_", fdr),
    paste0("conf.low_fdr_", fdr),
    paste0("conf.high_fdr_", fdr))
  
  data_plot <- data %>%
    select(all_of(cols_to_select)) %>%
    rename(model_type = model) %>%
    pivot_longer(
      cols = matches("estimate_|conf\\.low_|conf\\.high_"),
      names_to = c(".value", "term"),
      names_pattern = "(estimate|conf\\.low|conf\\.high)_(.*)") %>%
    mutate(term = recode(
      term,
      fdr_mom = "Mother",
      fdr_sib = "Sibling",
      fdr_dad = "Father",
      fdr_none = "None")) %>%
    mutate(
      or = exp(estimate),
      or_conf_low = exp(conf.low),
      or_conf_high = exp(conf.high)) %>%
    mutate(term = factor(term, levels = rev(fdr_order))) %>%
    mutate(outcome = case_when(
      outcome == "case_t1dprog" ~ "Progression",
      outcome == "ia" ~ "IA",
      outcome == "persist_conf_ab" ~ "IA",
      outcome == "t1d" ~ "T1D")) %>%
    mutate(outcome = factor(outcome, levels = outcome_order)) %>%
    mutate(
      model_type_ind = ifelse(grepl("GRS2x", model_type), "With GRS2x", "Without GRS2x"),
      model_type_ind = factor(model_type_ind, levels = c("With GRS2x", "Without GRS2x"))
    )
  
  plot <- ggplot(
    data_plot,
    aes(x = term, y = or, ymin = or_conf_low, ymax = or_conf_high,
        color = term, fill = term, alpha = model_type_ind)
  ) + 
    geom_linerange(position = position_dodge(width = 0.5), show.legend = FALSE, size = 0.8) +
    geom_hline(yintercept = 1, lty = 2) +
    geom_point(size = 2.5, shape = 21, colour = "white", stroke = 0.5,
               position = position_dodge(width = 0.5)) +
    scale_color_manual(values = c(
      "Father" = color_fdr_dad,
      "Mother" = color_fdr_mom,
      "Sibling" = color_fdr_sib,
      "None" = color_fdr_none),
      guide = "none") +
    scale_fill_manual(values = c(
      "Father" = color_fdr_dad,
      "Mother" = color_fdr_mom,
      "Sibling" = color_fdr_sib,
      "None" = color_fdr_none),
      guide = "none") +
    scale_alpha_manual(
      values = c("Without GRS2x" = 0.4, "With GRS2x" = 1)) +
    guides(
      alpha = guide_legend(
        reverse = T,
        override.aes = list(
          shape = 21,
          size  = 3,
          fill  = "grey50",
          color = "grey50"))) +
    coord_flip() +
    theme_bw() +
    labs(
      x = "",
      y = "Hazard Ratio"
    ) +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "white", color = "black"),
      axis.text.y = element_text(angle = 90, hjust = 0.5)
    ) +
    facet_wrap(~ model_group, ncol = 1, strip.position = "left")
  
  return(plot)
}


#TODO: this is the manuscript version, need to de-duplicate!
#TODO: remove DAISY & TEDDY support?
# data <- m2_dad
# ref <- "dad"
# term <- "GRS2x"
make_h2_forest_plot_manuscript <- function(data, ref, term) {
  
  fdr_all <- c("dad", "mom", "sib", "none")
  fdr <- setdiff(fdr_all, ref)
  
  cols_to_select <- c(
    "term", "outcome",
    paste0("estimate_fdr_", fdr),
    paste0("conf.low_fdr_", fdr),
    paste0("conf.high_fdr_", fdr))
  
  data_plot <- data %>%
    select(all_of(cols_to_select)) %>%
    dplyr::rename(model_type = term) %>%
    pivot_longer(
      cols = matches("estimate_|conf\\.low_|conf\\.high_"),
      names_to = c(".value", "term"),
      names_pattern = "(estimate|conf\\.low|conf\\.high)_(.*)") %>%
    dplyr::mutate(term = recode(
      term,
      fdr_mom = "Mother",
      fdr_sib = "Sibling",
      fdr_dad = "Father",
      fdr_none = "None")) %>%
    dplyr::mutate(
      or = exp(estimate),
      or_conf_low = exp(conf.low),
      or_conf_high = exp(conf.high)) %>%  # so plotting OR, not beta
    dplyr::mutate(term = factor(term, levels = rev(fdr_order))) %>%
    dplyr::mutate(outcome = factor(outcome, levels = outcome_order))
  
  data_plot <- data_plot %>%
    dplyr::filter(is.na(model_type) | model_type == !!term) %>%
    dplyr::filter(!duplicated(.)) %>%  # remove extra identical rows from without term
    dplyr::mutate(model_type_ind = ifelse(is.na(model_type), "Without GRS2x", "With GRS2x")) %>%
    dplyr::mutate(model_type_ind = factor(model_type_ind, levels = c("With GRS2x", "Without GRS2x")))
  
  plot <- ggplot(
    data_plot,
    aes(x=term, y=or, ymin=or_conf_low, ymax=or_conf_high, color=term,
        fill=term, alpha=model_type_ind)) + 
    geom_linerange(position=position_dodge(width = 0.5), show.legend = F, size=0.8) +
    geom_hline(yintercept = 1, lty=2) +
    geom_point(size=2.5, shape=21, colour="white", stroke = 0.5,
               position=position_dodge(width = 0.5)) +
    scale_color_manual(values = c(
      "Father" = color_fdr_dad,
      "Mother" = color_fdr_mom,
      "Sibling" = color_fdr_sib,
      "None" = color_fdr_none),
      guide = "none") +
    scale_fill_manual(values = c(
      "Father" = color_fdr_dad,
      "Mother" = color_fdr_mom,
      "Sibling" = color_fdr_sib,
      "None" = color_fdr_none),
      guide = "none") +
    scale_alpha_manual(
      values = c("Without GRS2x" = 0.4, "With GRS2x" = 1)) +
    guides(
      alpha = guide_legend(
        reverse = T,
        override.aes = list(
          shape = 21,
          size  = 3,
          fill  = "grey50",
          color = "grey50"))) +
    coord_flip() +
    theme_bw() +
    labs(
      x = "Affected Relative",
      y = "Hazard Ratio",
      alpha = "",
    ) +
    scale_x_discrete(labels = c("teddy" = "TEDDY","daisy" = "DAISY")) +
    theme(legend.position = "right",
          legend.direction = "vertical",
          strip.background = element_rect(fill = "white", color = "black"),
          axis.text.y = element_text(angle = 90, hjust = 0.5)) +
    facet_wrap(~outcome)
  
  return(plot)
  
}

#TODO: this is related to version in utils_coxph_tt, need to de-duplicate!
plot_tt_hr_manuscript <- function(
    data,
    var,
    ref_label,
    times = seq(0, 16, 0.5),
    add_last_matprot = FALSE,
    data_n = NULL) {
  
  data$cohort <- case_when(
    data$cohort == "None vs Dad" ~ "None",
    data$cohort == "Mom vs Dad" ~ "Mother",
    data$cohort == "Sib vs Dad" ~ "Sibling")
  data$cohort = factor(data$cohort, levels = fdr_order)
  
  p <- ggplot(data, aes(x = time, y = hr, color = cohort, fill = cohort)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
    scale_y_log10() +
    scale_x_continuous(breaks = times) +
    scale_color_manual(values = c(
      "Father" = color_fdr_dad,
      "Mother" = color_fdr_mom,
      "Sibling" = color_fdr_sib,
      "None" = color_fdr_none)) +
    scale_fill_manual(values = c(
      "Father" = color_fdr_dad,
      "Mother" = color_fdr_mom,
      "Sibling" = color_fdr_sib,
      "None" = color_fdr_none)) +
    labs(
      x = "Follow-up Time",
      y = "Hazard Ratio",
      color = "",
      fill = ""
    ) +
    theme_bw() +
    theme(
      legend.position = "right",
      legend.direction = "vertical",
      panel.grid.minor = element_blank())
  
  
  if (!is.null(data_n)) {
    times_n <- c(0, 3, 6, 9, 12, 15)
    
    plot_n <- data_n %>%
      dplyr::filter(time %in% times_n) %>%
      dplyr::group_by(time) %>%
      dplyr::summarize(
        n_remain = sum(n_remain))
    
    p <- p +
      geom_text(
        data = plot_n,
        aes(
          x = time,
          y = 0.1,   # adjust vertically as needed
          label = n_remain),
        inherit.aes = F,
        show.legend = F)
  }
  
  if (add_last_matprot) {
    # Get latest time T1D mom was protective, add point and CI line
    matprot_t <- data %>%
      dplyr::filter(cohort == "Mother") %>%
      dplyr::filter(upper < 1) %>%
      slice_max(upper, n = 1) %>%
      dplyr::pull(time)
    
    p <- p +
      geom_point(
        data = subset(data, cohort == "Mother" & time == matprot_t),
        size = 2.5,
        color = "black",
        show.legend = F
      ) +
      geom_linerange(
        data = subset(data, cohort == "Mother" & time == matprot_t),
        aes(ymin = lower, ymax = upper),
        linewidth = 0.75,
        color = "black",
        show.legend = F
      )
    
    
  }
  
  return(p)
  
}

