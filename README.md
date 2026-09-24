# MatProt PRS Family History Analysis

This repository was dervied from the internal project
[`rkjcollab/explore_matprot`](https://github.com/rkjcollab/explore_matprot),
specifically from the `genetics/famhistPRS` subdirectory at tag
`paper-release-famhistPRS` and commit `8765959`. The code was extracted using
git subtree split to preserve history.

This repository contains the analysis code for a manuscript that investigated
whether genetic susceptibility to T1D (measured by the GRS2x) explains that
relative maternal protection from T1D in offspring exposed to maternal T1D.

## Data

The repository is organized as a reproducible research compendium: the code is
here but the TEDDY data used cannot be made public and is stored outside the
repo at the RKJcollab data root `RKJCOLLAB`.

Before running the analysis, set `RKJCOLLAB` to the RKJcollab OneDrive data
root that contains top-level folders `Immunogenetics_T1D` and
`Maternal_Protection`.

## Setup

Restore the R environment:

```sh
Rscript -e 'renv::restore()'
```

PLINK1.9 and PLINK2 are also required but not managed by `renv`.

## Pipeline

Run the scripts in order from the repository root:

```sh
Rscript 01_make_ia_t1d_phenos.R
Rscript 02_lm_grs2_fdr_hla.R
Rscript 03_ia_risk_log_reg.R
Rscript 04_extract_grs2x_variants.R
Rscript 05_split_grs2x_variants_by_fdr.R
Rscript -e 'source("config.R"); rmarkdown::render("reports/manuscript.Rmd")'

```

- `01_make_ia_t1d_phenos.R`: builds TEDDY analytical phenotype files.
- `02_lm_grs2_fdr_hla.R`: runs the linear mixed model for GRS2x by FDR.
- `03_ia_risk_log_reg.R`: runs the survival models for IA, T1D, and
    progression from IA to T1D.
- `04_extract_grs2x_variants.R`: subsets TEDDY genotype data to GRS2x SNPs.
- `05_test_grs2x_variants_by_fdr.R`: tests for differences in GRS2x variant
    allele frequency by FDR.
- `reports/manuscript.Rmd`: renders the manuscript report.


## Outputs

Derived analysis files, reports, and results are written under `RKJCOLLAB`, in
the paths defined by `study_specs` in [R/study_specs.R](R/study_specs.R).

## Other Files

The input genetic data was imputed using the repo
[`rkjcollab/imputation`](https://github.com/rkjcollab/imputation) and the
config file [`imputation_config.yml`](imputation_config.yml) copied here but
run in the separate `rkjcollab/immuno_t1d` repo since used by multiple
projects.

The GRS2x was calculated using the script
[`reports/teddy_immunoT1D_tmr3_grs2x.Rmd`](reports/teddy_immunoT1D_tmr3_grs2x.Rmd),
again copied here but run in the separate
[`rkjcollab/immuno_t1d`](https://github.com/rkjcollab/immuno_t1d) repo.

The report [`reports/time_varying_analysis.Rmd`](reports/time_varying_analysis.Rmd)
is an archived exploratory report that documents a preliminary step in the
time-varying analysis. It is frozen at the point it was last rendered for
decision making and is kept for reference only.
