# data <- m2_dad
# ref <- "dad"
# term <- "GRS2x"
make_h2_forest_plot <- function(data, ref) {
  
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
