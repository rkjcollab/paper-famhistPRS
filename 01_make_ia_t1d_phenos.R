# SDS 20240627, updated 20260428

# Make TEDDY phenotype file for analysis of genetic risk in T1D moms. Start from
# Immunogenetics TEDDY R01 pheno with shared outcome definitions and exclusions:
# Immunogenetics_T1D/pheno/pheno_teddy_r01.tsv defined by script
# immuno_t1d/pheno/make-pheno-file-teddy-r01.R

# Shared criteria:
  # Have eligible HLA
  # From one of six primary clinical centers

# Additional analytical criteria
  # Have exome chip data - based on anc_d & anc_t inputs
  # Have primarily European ancestry
  # Have only one or no first-degree relatives with T1D
  # One genetically-identical twin removed at random
  # OPTIONAL: have only HLA DR3/4, DR4/4, DR3/3 or DR4/X

# Additional notes
  # Run analysis two ways: in everyone and in controls (no IA or T1D) only
  # Include 3- and 4-level FDR (dad/sib combined & separated)
  # Write out separate analytical pheno file for each outcome & subset
  # OPTIONAL: sex-stratified analyses

# Setup ------------------------------------------------------------------------

library(here)

source(here("config.R"))
source(here("R/study_specs.R"))
source(here("R/utils_pheno.R"))
source(here("R/utils_teddy.R"))

# TO NOTE: change study list and other settings in config.R

# Define function --------------------------------------------------------------

# Code automatically makes files with all individuals and with controls only,
# and applies all criteria above except for optional HLA-DR filter.
# dr_filt is defined in config.R. Set to "yes" or "no" to filter
# DR3/4, DR4/4, DR3/3 or DR4/X.

make_pheno <- function(study, config) {
  if (!study %in% names(study_specs)) {
    stop("Unknown study: ", study)
  }
  specs <- study_specs[[study]]
  
  # Do study-specific prep
  df <- switch(
    study,
    teddy = prep_teddy(specs),
    stop("Unknown study")
  )
  message("Study: ", study)
  message("After prep: ", nrow(df))
  
  # Do shared prep
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
    study = study,
    out_prefix = specs$intermed_out_dir,
    dr_suffix = dr_suffix
  )
}


# Run function -----------------------------------------------------------------

#TODO: switch to only single study?

for (s in config$studies) {
  make_pheno(s, config)
}
