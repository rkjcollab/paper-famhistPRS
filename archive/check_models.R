# SDS 20251125

# Adding separate script to check specific questions about models run in
# 3_lm_grs2_fdr_hla.R and 4_ia_risk_log_reg.R. 

# Setup ------------------------------------------------------------------------

library(tidyverse)
setwd(Sys.getenv("RKJCOLLAB"))

# DAISY inputs
d_all_ia <- readRDS(
  "Maternal_Protection/data/daisy/famhist_prs/pheno_all_ia.rds")  # pheno for hypoth 1
d_all_t1d <- readRDS(
  "Maternal_Protection/data/daisy/famhist_prs/pheno_all_t1d.rds")
d_all_prog <- readRDS(
  "Maternal_Protection/data/daisy/famhist_prs/pheno_all_prog.rds")
d_ctrls_ia <- readRDS(
  "Maternal_Protection/data/daisy/famhist_prs/pheno_ctrls_ia.rds")  # pheno for hypoth 1
d_ctrls_t1d <- readRDS(
  "Maternal_Protection/data/daisy/famhist_prs/pheno_ctrls_t1d.rds")
d_ctrls_prog <- readRDS(
  "Maternal_Protection/data/daisy/famhist_prs/pheno_ctrls_prog.rds")
d_kinship <- readRDS(
  "DAISY/genetics/daisy_ask_genetics/genesis/study_nhw_GRM.rds")

# TEDDY inputs
t_all <- readRDS(
  "Maternal_Protection/data/teddy/famhist_prs/pheno_all.rds")  # pheno for hypoth 1
t_all_prog <- readRDS(
  "Maternal_Protection/data/teddy/famhist_prs/pheno_all_prog.rds")
t_ctrls_ia <- readRDS(
  "Maternal_Protection/data/teddy/famhist_prs/pheno_ctrls_ia.rds")
t_ctrls_t1d <- readRDS(
  "Maternal_Protection/data/teddy/famhist_prs/pheno_ctrls_t1d.rds")
t_ctrls_prog <- readRDS(
  "Maternal_Protection/data/teddy/famhist_prs/pheno_ctrls_prog.rds")
t_kinship <- readRDS(
  "Immunogenetics_T1D/genetics/teddy_r01/genesis_2/study_nhw_GRM.rds")

# Hypoth 2 results
m2 <- read_csv(
  "Maternal_Protection/data/results/famhist_prs/lm_outcome_fdr_none_ref_daisy_teddy.csv") 
m2_sex <- read_csv(
  "Maternal_Protection/data/results/famhist_prs/lm_outcome_fdr_none_ref_daisy_teddy_with_sex.csv")


# Check kinship in prog analysis -----------------------------------------------

# All analysis

### DAISY
# Check if need to include kinship in the subset
d_kinship_mat <- as.matrix(d_kinship)
d_kinship_mat <- d_kinship_mat[d_all_prog$ID, d_all_prog$ID]

# Check if any non-0 values not on diagonal
d_kinship_ct <- sum(d_kinship_mat[lower.tri(d_kinship_mat)] > 0)
# 18 non-0 pairs - so do include in analysis

### TEDDY
# Check if need to include kinship in the subset
t_kinship_mat <- as.matrix(t_kinship)
t_kinship_mat <- t_kinship_mat[t_all_prog$ID, t_all_prog$ID]

# Check if any non-0 values not on diagonal
t_kinship_ct <- sum(t_kinship_mat[lower.tri(t_kinship_mat)] > 0)
  # 27 non-0 pairs - so do include in analysis

# Controls analysis
### DAISY
# Check if need to include kinship in the subset
d_kinship_mat <- d_kinship_mat[d_ctrls_prog$ID, d_ctrls_prog$ID]

# Check if any non-0 values not on diagonal
d_kinship_ct <- sum(d_kinship_mat[lower.tri(d_kinship_mat)] > 0)
# 5 non-0 pairs - so do include in analysis

### TEDDY
# Check if need to include kinship in the subset
t_kinship_mat <- t_kinship_mat[t_ctrls_prog$ID, t_ctrls_prog$ID]

# Check if any non-0 values not on diagonal
t_kinship_ct <- sum(t_kinship_mat[lower.tri(t_kinship_mat)] > 0)
# 11 non-0 pairs - so do include in analysis


# Check if sex needed in hypoth 2 ----------------------------------------------

# Check if any sex term p-values significant
m2_sex %>% dplyr::filter(p.value_sex1 < 0.05) %>%
  dplyr::select(study, outcome, term, nobs, p.value_sex1, estimate_sex1)
#   study outcome term   nobs p.value_sex1 estimate_sex1
# 1 teddy IA      NA     6942       0.0228         0.172
# 2 teddy IA      GRS2x  6942       0.0115         0.194
# 3 teddy IA      NA     6942       0.0228         0.172
# 4 teddy IA      dr34   6942       0.0174         0.180

# Only significant in TEDDY IA models...
# TODO: ask Randi about including

