# Script defines study-specific paths

# Add lists here for each study

base_dir <- Sys.getenv("RKJCOLLAB")

study_specs <- list(
  teddy = list(
    # Inputs
    pheno_raw_path =
      paste0(base_dir, "/Immunogenetics_T1D/pheno/pheno_teddy_r01.tsv"),
    fid_path =
      paste0(base_dir, "/Immunogenetics_T1D/raw/teddy_r01/2025-03-05/OmicsDatasets/t1dexome_masked.fam"),
    twin_list_path =
      paste0(base_dir, "/Immunogenetics_T1D/genetics/teddy_r01/genesis_1/study_full_twin_to_rm.txt"),
    grs2_path =
      paste0(base_dir, "/Immunogenetics_T1D/genetics/teddy_r01/risk_scores/grs2x/tm_r3_imp/full_qc/score_results_with_imp.csv"),
    pcair_path =
      paste0(base_dir, "/Immunogenetics_T1D/genetics/teddy_r01/genesis_2/study_nhw_pcair.rds"),
    kinship_path =
      paste0(base_dir, "/Immunogenetics_T1D/genetics/teddy_r01/genesis_2/study_nhw_GRM.rds"),

    # Outputs
    intermed_out_dir = paste0(base_dir, "/Maternal_Protection/data/teddy/famhist_prs"),
    result_out_dir = paste0(base_dir, "/Maternal_Protection/data/results/famhist_prs"),
    report_out_dir = paste0(base_dir, "/Maternal_Protection/reports/famhist_prs")
  )
)

