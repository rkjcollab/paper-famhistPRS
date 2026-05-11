# MatProt PRS Family History Analysis

This repository contains the analysis code for a manuscript evaluating whether type 1 diabetes polygenic risk scores explain first-degree relative patterns in TEDDY and related analyses.

The repository is organized as a reproducible research compendium: the code is version controlled here, while private/raw study data are expected to live outside the repository under the collaboration data root defined by `RKJCOLLAB`.

## Data Access

Raw and derived study data are not included in this repository. Before running the analysis, set `RKJCOLLAB` to the private collaboration data root that contains directories such as `Immunogenetics_T1D`, `Maternal_Protection`, and `DAISY`.

For example:

```sh
export RKJCOLLAB=/path/to/collaboration/root
```

The scripts use [config.R](config.R) to resolve repo-local helper files and private data paths.

## Main Workflow

The publication workflow is encoded in [Makefile](Makefile). From the repository root, run:

```sh
make help
make check-env
make all
```

The main targets are:

- `make phenotypes`: builds TEDDY analytical phenotype files with the manuscript exclusions and derived variables.
- `make data-summary`: renders the analytical data summary report.
- `make lm-grs`: runs the linear mixed model for `GRS2x ~ FDR + covariates`.
- `make outcome-models`: runs the survival outcome models for IA, T1D, and progression.
- `make manuscript`: renders the manuscript report.
- `make time-varying`: renders the time-varying analysis report.

## Key Files

- [1_make_ia_t1d_phenos.R](1_make_ia_t1d_phenos.R): creates analytical phenotype files.
- [2_analytical_data_summary.Rmd](2_analytical_data_summary.Rmd): summarizes analytical datasets.
- [3_lm_grs2_fdr_hla.R](3_lm_grs2_fdr_hla.R): models PRS differences by family history.
- [4_ia_risk_log_reg.R](4_ia_risk_log_reg.R): models IA, T1D, and progression outcomes.
- [Manuscript_FamhistPRS.Rmd](Manuscript_FamhistPRS.Rmd): final manuscript report source.
- [R/study_specs.R](R/study_specs.R): study-specific input paths, phenotype paths, covariates, and survival definitions.
- [mod_lmekin.R](mod_lmekin.R), [mod_coxme.R](mod_coxme.R), and [mod_coxph.R](mod_coxph.R): model helper functions.
- [set_paper_colors.R](set_paper_colors.R): manuscript plotting colors.
- [0_prs_famhist_analytical_plan.R](0_prs_famhist_analytical_plan.R): exploratory analytical planning notes.
- [archive](archive): historical/provenance code retained for reference.
- [PRSedm](PRSedm): polygenic risk score tooling and resources used by the project.

## Outputs

Derived analysis files and rendered reports are written under `RKJCOLLAB`, primarily:

- `Maternal_Protection/data/teddy/famhist_prs`
- `Maternal_Protection/data/results/famhist_prs`
- `Maternal_Protection/reports/risk_scores`

## Notes For Reuse

This repository is intended to document and reproduce the analyses supporting the manuscript. It is not yet packaged as a general-purpose R package. The immediate reproducibility contract is the Makefile workflow plus the shared path configuration in [config.R](config.R).
