# One-time preprocessing step to subset genotype data to GRS2x SNPs for
# manuscript analysis (not for GRS2x score generation).

# Setup ------------------------------------------------------------------------

library(here)
source(here("R/study_specs.R"))
library(tidyverse)
plink <- "plink2"

snp_list_path <- here("PRSedm/snplists/grs2_version_snplists.xlsx")
risk_score_dir <- dirname(study_specs$teddy$grs2_path)
plink_path <- study_specs$teddy$plink_path

# Write SNP list ---------------------------------------------------------------

snp_list <- read_xlsx(snp_list_path, sheet = "GRS2X_TOPMED_R3")

snp_id <- snp_list %>%
  dplyr::mutate(id = paste0(POSITION_HG38, ":", REF, ":", ALT)) %>%
  dplyr::select(id)

# Write out SNP list
readr::write_tsv(
  snp_id,
  paste0(risk_score_dir, "/grs2x_snp_ids.txt"),
  col_names = FALSE)

# Run PLINK --------------------------------------------------------------------

plink_args <- c(
  "--pfile", study_specs$teddy$plink_path,
  "--extract", paste0(risk_score_dir, "/grs2x_snp_ids.txt"),
  "--make-pgen",
  "--out", file.path(risk_score_dir, "/chr_all_concat_grs2x"))

system2(plink, plink_args)
