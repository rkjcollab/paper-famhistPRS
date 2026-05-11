library(dplyr)

set_factor_levels <- function(data) {
  # Keep outcomes as 0/1 for coxme/surv
  data <- data %>%
    mutate(
      mAA_at_IA = relevel(factor(mAA_at_IA), ref = "0"),
      sex = relevel(factor(sex), ref = "0"),
      dr34 = relevel(factor(dr34), ref = "0"),
      fdr_3level = relevel(factor(fdr_3level), ref = "None"),
      HLAGRP = relevel(factor(HLAGRP), ref = "DR4/X")
    )
  
  if ("cc" %in% names(data)) {
    data$cc <- relevel(factor(data$cc), ref = 1)
  }
  
  if ("fdr_4level" %in% names(data)) {
    data$fdr_4level <- relevel(factor(data$fdr_4level), ref = "None")
  }
  
  data
}

# Apply HLD-DR filter to keep only HLA DR3/4, DR4/4, DR3/3 or DR4/X
apply_dr_filt <- function(data) {
  dr_list <- c("DR3/3", "DR3/4", "DR4/4", "DR4/X")
  dplyr::filter(data, HLAGRP %in% dr_list)
}

# Functions for subsetting exports based on outcome
get_all <- function(df, outcome) {
  switch(
    outcome,
    IA = df %>% filter(!is.na(IA)),
    T1D = df %>% filter(!is.na(T1D)),
    T1D_prog = df %>% filter(!is.na(fupT1D_prog))
  )
}

# Functions for subsetting exports to controls only
get_controls <- function(df, outcome) {
  switch(
    outcome,
    IA = df %>% filter(IA == 0 & (T1D == 0 | is.na(T1D))),
    T1D = df %>% filter(T1D == 0),
    T1D_prog = df %>% filter(T1D == 0)
  )
}

# Note that data should be coded so male = 1, female = 0
split_sex <- function(df) {
  list(
    female = df %>% filter(sex == 0),
    male = df %>% filter(sex == 1)
  )
}

make_exports <- function(df) {
  out <- list()

  for (outcome in c("IA", "T1D")) {
    df_all <- get_all(df, outcome)
    df_ctrl <- get_controls(df_all, outcome)
    sex_split <- split_sex(df_all)

    out[[paste0("all_", tolower(outcome))]] <- df_all
    out[[paste0("ctrls_", tolower(outcome))]] <- df_ctrl
    out[[paste0("female_", tolower(outcome))]] <- sex_split$female
    out[[paste0("male_", tolower(outcome))]] <- sex_split$male
  }

  df_prog <- get_all(df, "T1D_prog")
  df_prog_ctrl <- get_controls(df_prog, "T1D_prog")
  sex_prog <- split_sex(df_prog)

  out$all_prog <- df_prog
  out$ctrls_prog <- df_prog_ctrl
  out$female_prog <- sex_prog$female
  out$male_prog <- sex_prog$male

  out
}

write_phenos <- function(pheno_list, study, out_prefix, dr_suffix) {
  dir.create(out_prefix, recursive = TRUE, showWarnings = FALSE)

  for (nm in names(pheno_list)) {
    saveRDS(
      pheno_list[[nm]],
      file = file.path(out_prefix, paste0("pheno_", nm, dr_suffix, ".rds"))
    )
  }
}
