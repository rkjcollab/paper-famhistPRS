# Script defines paths and settings for the TEDDY study analysis.

base_dir <- Sys.getenv("RKJCOLLAB")

# Study-specific paths and settings
study_specs <- list(
  # Label used in output file names
  study = "teddy",

  # Inputs
  pheno_raw_path =
    paste0(base_dir, "/Immunogenetics_T1D/pheno/pheno_teddy_r01.tsv"),
  fid_path =
    paste0(base_dir, "/Immunogenetics_T1D/raw/teddy_r01/2025-03-05/OmicsDatasets/t1dexome_masked.fam"),
  twin_list_path =
    paste0(base_dir, "/Immunogenetics_T1D/genetics/teddy_r01/genesis_1/study_full_twin_to_rm.txt"),
  ancestry_path =
    paste0(base_dir, "/Immunogenetics_T1D/genetics/teddy_r01/admixture/ancestry_estimation/study_all_pop.txt"),
  grs2_path =
    paste0(base_dir, "/Immunogenetics_T1D/genetics/teddy_r01/risk_scores/grs2x/tm_r3_imp/full_qc/score_results_with_imp.csv"),
  pcair_path =
    paste0(base_dir, "/Immunogenetics_T1D/genetics/teddy_r01/genesis_2/study_nhw_pcair.rds"),
  kinship_path =
    paste0(base_dir, "/Immunogenetics_T1D/genetics/teddy_r01/genesis_2/study_nhw_GRM.rds"),
  plink_path =
    paste0(base_dir, "/Immunogenetics_T1D/genetics/teddy_r01/imputation/tm_r3_imp/exome_chip/full_qc/imputed_clean_maf0_rsq0.3/chr_all_concat"),

  # Outputs
  intermed_out_dir = paste0(base_dir, "/Maternal_Protection/data/teddy/famhist_prs"),
  result_out_dir = paste0(base_dir, "/Maternal_Protection/data/results/famhist_prs"),
  report_out_dir = paste0(base_dir, "/Maternal_Protection/reports/famhist_prs"),
  manuscript_out_dir = paste0(base_dir, "/Maternal_Protection/dissemination/famhist_prs_paper"),

  # Settings for survival model in script 03
  # TO NOTE: all models are hardcoded to account for TEDDY relatedness, either
  # with cluster(FID) (coxph, coxph_tt) or with kinship as a random effect
  # (coxme).
  surv_def = list(
    IA = list(time = "fupIA", event = "IA"),
    T1D = list(time = "fupT1D", event = "T1D"),
    T1D_prog = list(time = "fupT1D_prog", event = "T1D")
  ),
  surv_covs = list(
    IA = c("fdr_4level", "sex", "PC1", "PC2", "cc"),
    T1D = c("fdr_4level", "sex", "PC1", "PC2", "cc"),
    T1D_prog = c("fdr_4level", "sex", "PC1", "PC2", "cc", "fupIA", "mAA_at_IA")
  )
)
