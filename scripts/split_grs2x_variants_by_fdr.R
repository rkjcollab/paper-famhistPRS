# One-time preprocessing step to subset GRS2x genotype data individuals with
# T1D mothers or T1D fathers (not for GRS2x score generation).
# Script "extract_grs2x_variants.R" should be run before this one.

# Setup ------------------------------------------------------------------------

library(here)
devtools::load_all()
source(here("config.R"))
library(tidyverse)
library(readxl)
plink <- "plink2"
plink1 <- "plink"

pheno_path <- paste0(study_specs$teddy$intermed_out_dir, "/pheno_all_ia.rds")
risk_score_dir <- dirname(study_specs$teddy$grs2_path)
# same as pheno_all_t1d.rds
snp_list_path <- here("PRSedm/snplists/grs2_version_snplists.xlsx")

# Write lists ------------------------------------------------------------------

pheno <- readRDS(pheno_path)

pheno_d <- pheno %>%
  dplyr::filter(fdr_4level == "Dad") %>%
  dplyr::select(FID, ID)
pheno_m <- pheno %>%
  dplyr::filter(fdr_4level == "Mom") %>%
  dplyr::select(FID, ID)

write_tsv(pheno_d, paste0(risk_score_dir, "/t1d_dad_ids.txt"), col_names = F)
write_tsv(pheno_m, paste0(risk_score_dir, "/t1d_mom_ids.txt"), col_names = F)

# Make mom vs. dad file for --glm
pheno_plink <- rbind(
  pheno_d %>% dplyr::mutate(PHENO = 1),
  pheno_m %>% dplyr::mutate(PHENO = 2))
write_tsv(pheno_plink, paste0(risk_score_dir, "/t1d_dad_mom_pheno_assoc.txt"), col_names = F)

# Run PLINK --------------------------------------------------------------------

# Get AFs for each group separately
plink_args <- c(
  "--pfile", paste0(risk_score_dir, "/chr_all_concat_grs2x"),
  "--keep", paste0(risk_score_dir, "/t1d_dad_ids.txt"),
  "--freq",
  "--out", paste0(risk_score_dir, "/chr_all_concat_grs2x_t1d_dad"))
system2(plink, plink_args)

plink_args <- c(
  "--pfile", paste0(risk_score_dir, "/chr_all_concat_grs2x"),
  "--keep", paste0(risk_score_dir, "/t1d_mom_ids.txt"),
  "--freq",
  "--out", paste0(risk_score_dir, "/chr_all_concat_grs2x_t1d_mom"))
system2(plink, plink_args)

# Make file with only mom and dads with T1D
plink_args <- c(
  "--pfile", paste0(risk_score_dir, "/chr_all_concat_grs2x"),
  "--keep", paste0(risk_score_dir, "/t1d_dad_ids.txt"),
    paste0(risk_score_dir, "/t1d_mom_ids.txt"),
  "--make-bed",
  "--out", paste0(risk_score_dir, "/chr_all_concat_grs2x_t1d_dad_mom"))
system2(plink, plink_args)

# Make table with AFs ----------------------------------------------------------

# Read in PLINK results
af_d <- read_delim(
  paste0(risk_score_dir, "/chr_all_concat_grs2x_t1d_dad.afreq"))
af_m <- read_delim(
  paste0(risk_score_dir, "/chr_all_concat_grs2x_t1d_mom.afreq"))

# Read in GRS2x SNP list
snp_list <- read_xlsx(snp_list_path, sheet = "GRS2X_TOPMED_R3")

# First, get AFs split by sex
af_fam <- inner_join(
  af_m %>% dplyr::select(-`#CHROM`, -OBS_CT),
  af_d %>% dplyr::select(-`#CHROM`, -OBS_CT),
  by = c("ID", "REF", "ALT"),
  suffix = c("_mother", "_father"))  # all 64, so same ALT

af_fam_grs2x <- left_join(
  snp_list %>% dplyr::mutate(ID = paste0(POSITION_HG38, ":", REF, ":", ALT)),
  af_fam,
  by = c("ID", "REF", "ALT"))

write_tsv(af_fam_grs2x, paste0(risk_score_dir, "/afreq_by_fdr_grs2x.txt"))

# Run assoc allelic chi-square test --------------------------------------------

plink_args <- c(
  "--bfile", paste0(risk_score_dir, "/chr_all_concat_grs2x_t1d_dad_mom"),
  "--pheno", paste0(risk_score_dir, "/t1d_dad_mom_pheno_assoc.txt"),
  "--assoc",
  "--keep-allele-order",
  "--out", paste0(risk_score_dir, "/t1d_dad_mom_pheno_assoc_result"))
system2(plink1, plink_args)

