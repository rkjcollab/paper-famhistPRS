# SDS 20250925

# Script while building analytical plan for MatProt PRS paper.

# 1. Relatedness in TEDDY
  # Want to figure out how to handle the relatedness in TEDDY -
  # so want to understand number of siblings, what is their family history, what
  # is their outcome.

# 2. Check if TEDDY outcome varies by clinical center, by PCs, or by sex
  # Want to check this in analytical dataset.

# 3. Check if kinship is needed for TEDDY and DAISY in progression analysis
  # subset.

# Setup ------------------------------------------------------------------------

library(tidyverse)
setwd(paste0(Sys.getenv("RKJCOLLAB"), "/Immunogenetics_T1D"))

pheno <- sas7bdat::read.sas7bdat(
  "raw/teddy_r01/2025-09-17/rj_immuno_all_demo_masked_v2.sas7bdat")
  # No family ID included here

# Load genetic fam file for family ID
fam <- read_delim(
  "raw/teddy_r01/2025-03-05/OmicsDatasets/t1dexome_masked.fam",
  col_names = c("fid", "iid", "pid", "mid", "sex", "pheno"))

# Load analytical dataset in all
data_t <- readRDS("../Maternal_Protection/data/teddy/famhist_prs/pheno_all.rds")
data_d <- readRDS("../Maternal_Protection/data/daisy/famhist_prs/pheno_all.rds")

kin_t <- readRDS(
  "../Immunogenetics_T1D/genetics/teddy_r01/ancestry_estimation/TEDDY_GRM.rds")
kin_d <- readRDS(
  "../DAISY/genetics/daisy_ask_genetics/genesis/study_nhw_GRM.rds")

# 1. Relatedness in TEDDY ------------------------------------------------------

# Add FDR 4-level variable, and filter to remove individuals with more than
# one FDR
# family_mem: 0=none, 1=all, 2=mom & dad, 3=mom & sibling, 4=dad & sibling,
# 5=mom, 6=dad, 7=sibling
table(pheno$family_mem_screen)
# 0    2    3    4    5    6    7 
# 7747   11    7    9  335  444  123

# Don't have mom, dad, sibling values available, so need to recreate
pheno$mom = ifelse(
  pheno$family_mem_screen == 2 |
    pheno$family_mem_screen == 3 |
    pheno$family_mem_screen == 5, 1, 0)
pheno$dad = ifelse(
  pheno$family_mem_screen == 2 |
    pheno$family_mem_screen == 4 |
    pheno$family_mem_screen == 6, 1, 0)
pheno$sibling = ifelse(
  pheno$family_mem_screen == 3 |
    pheno$family_mem_screen == 4 |
    pheno$family_mem_screen == 7, 1, 0)

# Remove those with 2, 3, 4, add FDR 4-level
pheno_filt_4level <- pheno %>%
  dplyr::filter(
    family_mem_screen != 2 & family_mem_screen != 3 & family_mem_screen != 4) %>%
  dplyr::mutate(fdr_4level = case_when(
    mom == 1 ~ "Mom",
    dad == 1 ~ "Dad",
    sibling == 1 ~ "Sib",
    family_mem_screen == 0 ~ "None"))
table(pheno_filt_4level$family_mem_screen)
# 0    5    6    7 
# 7747  335  444  123 
table(pheno_filt_4level$fdr_4level)
# Dad  Mom None  Sib 
# 444  335 7747  123 

# Find siblings
#TODO: not sure whether to focus on fid, mid, or all?
n_distinct(fam$fid)  # 7683
n_distinct(fam$iid)  # 8215
n_distinct(fam$pid)  # 7700
n_distinct(fam$mid)  # 7702

# Merge
pheno_fam <- inner_join(
  pheno_filt_4level,
  fam,
  by = c("RJohnson_Immuno_anc_MaskID" = "iid"))

# Get fid dups
pheno_fam_fid_dup_list <- pheno_fam %>%
  dplyr::filter(duplicated(fid))
pheno_fam_fid_dups <- pheno_fam %>%
  dplyr::filter(fid %in% pheno_fam_fid_dup_list$fid)

# How many fids appear more than once?
n_distinct(pheno_fam_fid_dups$fid)  # 503
# How many kids are in this?
n_distinct(pheno_fam_fid_dups$RJohnson_Immuno_anc_MaskID)  # 1027

# Family history of these kids?
table(pheno_fam_fid_dups$fdr_4level)
# Dad  Mom None  Sib 
# 63   62  894    8 
table(pheno_fam$fdr_4level)
# Dad  Mom None  Sib 
# 424  319 7327  120 

# Outcomes of these kids?
table(pheno_fam_fid_dups$t1d)
# 0   1 
# 978  49 
table(pheno_fam_fid_dups$persist_conf_ab)
# 0   1 
# 898 129 

# Check if score is normal
score <- read_csv(
  "genetics/teddy_r01/risk_scores/grs2x/score_results_imputed.csv")

ggplot(score, aes(x = `t1dgrs2-luckett25_total`)) +
  geom_histogram(bins = 100)

# Both family history and outcome?
table(pheno_fam_fid_dups$fdr_4level, pheno_fam_fid_dups$t1d)
#       0   1
# Dad   57   6
# Mom   58   4
# None 855  39
# Sib    8   0

table(pheno_fam_fid_dups$fdr_4level, pheno_fam_fid_dups$persist_conf_ab)
#       0   1
# Dad   50  13
# Mom   55   7
# None 786 108
# Sib    7   1


# 2. Check TEDDY outcome -------------------------------------------------------

##### Check if outcome varies by clinical center --------
### IA
data_t_tab <- table(data_t$IA, data_t$cc)
data_t_tab
#     1    2    3    4    5    6    133  134
# 0  706  616  928 1493  457 1986    2    7
# 1   97   60   83  205   54  291    0    0
chisq.test(data_t_tab)
# Pearson's Chi-squared test
# data:  data_tab
# X-squared = 21.475, df = 7, p-value = 0.003128

data_t %>% 
  dplyr::filter(cc == 133 | cc == 134) %>%
  dplyr::select(fdr, country)
# All FDR, all US

# Re-run without these
data_t_filt <- data_t %>%
  dplyr::filter(cc != 133 & cc != 134) %>%
  droplevels()
data_t_filt_tab <- table(data_t_filt$IA, data_t_filt$cc)
data_t_filt_tab
chisq.test(data_t_filt_tab)
# Pearson's Chi-squared test
# data:  data_filt_tab
# X-squared = 20.303, df = 5, p-value = 0.001096

# Also check country in all
data_t_tab_2 <- table(data_t$IA, data_t$country)
data_t_tab_2
#     1    2    3    4
# 0 2259 1493  457 1986
# 1  240  205   54  291
chisq.test(data_t_tab_2)
# Pearson's Chi-squared test
# data:  data_tab_2
# X-squared = 13.424, df = 3, p-value = 0.003804

### T1D
data_t_tab <- table(data_t$T1D, data_t$cc)
data_t_tab
#     1    2    3    4    5    6    133  134
# 0  739  644  975 1600  475 2168    2    7
# 1   64   32   36   98   36  109    0    0
chisq.test(data_t_tab)
# Pearson's Chi-squared test
# data:  data_tab
# X-squared = 21.475, df = 7, p-value = 0.003128

# Re-run without these
data_t_filt_tab <- table(data_t_filt$T1D, data_t_filt$cc)
data_t_filt_tab
chisq.test(data_t_filt_tab)
# Pearson's Chi-squared test
# data:  data_filt_tab
# X-squared = 22.594, df = 5, p-value = 0.0004036

# Also check country in all
data_t_tab_2 <- table(data_t$T1D, data_t$country)
data_t_tab_2
#     1    2    3    4
# 0 2367 1600  475 2168
# 1  132   98   36  109
chisq.test(data_t_tab_2)
# Pearson's Chi-squared test
# data:  data_tab_2
# X-squared = 4.9221, df = 3, p-value = 0.1776

# Conclusion: with or without odd clinical center entries, both IA and T1D
# vary by clinical center. Only IA also varies by country.

##### Check if outcome varies by PCs -------------------------------------------
### IA
mod <- glm(IA ~ PC1, data = data_t, family = binomial)
summary(mod)
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -2.06339    0.03791 -54.422   <2e-16 ***
#   PC1         -5.37781    2.89732  -1.856   0.0634 . 

### T1D
mod <- glm(T1D ~ PC1, data = data_t, family = binomial)
summary(mod)
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -2.87187    0.05322 -53.959   <2e-16 ***
#   PC1         -3.71119    4.09311  -0.907    0.365

# Conclusion: neither IA or T1D significantly vary with PC1.


##### Check if outcome varies by sex -------------------------------------------

### IA
data_t_tab <- table(data_t$IA, data_t$sex)
data_t_tab
#   Female Male
# 0   3070 3125
# 1    356  434
chisq.test(data_t_tab)
# Pearson's Chi-squared test with Yates' continuity correction
# data:  data_tab
# X-squared = 5.4808, df = 1, p-value = 0.01923

### T1D
data_t_tab <- table(data_t$T1D, data_t$sex)
data_t_tab
#   Female Male
# 0   3256 3354
# 1    170  205
chisq.test(data_t_tab)
# Pearson's Chi-squared test with Yates' continuity correction
# data:  data_tab
# X-squared = 2.0337, df = 1, p-value = 0.1538


# 3. Check kinship in prog -----------------------------------------------------

# For now, just checking in IA cases
data_d_ia <- data_d %>%
  dplyr::filter(IA == 1)  # 167
data_t_ia <- data_t %>%
  dplyr::filter(IA == 1)  # 790

# Subset kinship, keeping order of rows/cols
ids_d <- data_d_ia$ID  # the ID column used in (1|ID)
kin_d <- kin_d[ids_d, ids_d]  # reorder rows/cols
ids_t <- data_t_ia$ID
kin_t <- kin_t[ids_t, ids_t]

# Only get upper triangle, not including diagonal
kin_d_tri <- kin_d[upper.tri(kin_d)]
kin_t_tri <- kin_t[upper.tri(kin_t)]

# Count how many pairs > 0.2 - siblings
sum(kin_d_tri > 0.2)  # 14
sum(kin_t_tri > 0.2)  # 11
