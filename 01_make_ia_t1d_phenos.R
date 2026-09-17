# SDS 20240627

# Make TEDDY phenotype file for analysis of genetic risk in T1D moms. Start from
# Immunogenetics TEDDY R01 pheno with shared outcome definitions and exclusions:
# Immunogenetics_T1D/pheno/pheno_teddy_r01.tsv defined by script
# immuno_t1d/pheno/make-pheno-file-teddy-r01.R

# Shared criteria:
# Have eligible HLA
# From one of six primary clinical centers

# Additional analytical criteria
# Have exome chip data
# Have primarily European ancestry
# Have only one or no first-degree relatives with T1D
# One genetically-identical twin removed at random
# OPTIONAL: have only HLA DR3/4, DR4/4, DR3/3 or DR4/X

# Additional notes
# Subset to everyone or to controls (no IA or T1D) only
# Use 3- and 4-level FDR (dad/sib combined & separated)
# Can run sex-stratified analyses
# Writes out separate analytical pheno file for each outcome & subset

# Setup ------------------------------------------------------------------------

library(here)
devtools::load_all()
source(here("config.R"))

# Define function --------------------------------------------------------------

# Code automatically makes files with all individuals and with controls only,
# and applies critera as set in config.R.

make_pheno <- function(config) {
  # Do TEDDY-specific prep
  df_tmp <- prep_teddy_base(study_specs)
  df <- prep_teddy_final(study_specs, df_tmp)

  message("After prep: ", nrow(df))

  # Do general prep
  df <- set_factor_levels(df)

  if (config$dr_filt == "yes") {
    df <- apply_dr_filt(df)
    dr_suffix <- "_dr_filt"
    message("After DR filter: ", nrow(df))
  } else {
    dr_suffix <- ""
  }

  # Write exports
  pheno_list <- make_exports(df)

  write_phenos(
    pheno_list,
    out_prefix = study_specs$intermed_out_dir,
    dr_suffix = dr_suffix
  )
}


# Run function -----------------------------------------------------------------

make_pheno(config)
