# Script gets the last age-specific HR that shows significant maternal
# protection from the final coxph() model with linear time-transform of
# fdr_4level and with GRS2x adjustment.

# Setup ------------------------------------------------------------------------

library(here)
source(here("R/study_specs.R"))
source(here("R/utils_coxph_tt_hr.R"))
library(tidyverse)

result_dir <- study_specs$teddy$result_out_dir

model <- readRDS(paste0(
  result_dir,
  "/model_objects/coxph_tt_linear_fdr_4level_IA_fdr_dad_ref_fdr_4level-PC1-PC2-sex-cc-cluster(FID)-GRS2x_teddy_all.rds"))
model <- readRDS(paste0(
  result_dir,
  "/model_objects/coxph_tt_linear_fdr_4level_IA_fdr_dad_ref_fdr_4level-PC1-PC2-sex-cc-cluster(FID)_teddy_all.rds"))

# Get PH over time -------------------------------------------------------------

tt_hr <- get_tt_hr(
  mod = model,
  var = "fdr_4level",
  levels = c("None", "Mom", "Sib"),
  ref_label = "Dad",
  times = seq(0, 16, 0.5),
  time_transform = identity)


#TODO: add write out of results

