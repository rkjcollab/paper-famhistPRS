.PHONY: help check-env all phenotypes data-summary lm-grs outcome-models manuscript time-varying

R ?= Rscript

help:
	@echo "MatProt PRS analysis workflow"
	@echo ""
	@echo "Required environment:"
	@echo "  RKJCOLLAB  Private collaboration data root"
	@echo ""
	@echo "Targets:"
	@echo "  make check-env       Check required environment variables"
	@echo "  make phenotypes      Build TEDDY analytical phenotype files"
	@echo "  make data-summary    Render analytical data summary"
	@echo "  make lm-grs          Run GRS2x/FDR linear mixed model"
	@echo "  make outcome-models  Run outcome survival models"
	@echo "  make manuscript      Render manuscript report"
	@echo "  make time-varying    Render time-varying analysis report"
	@echo "  make all             Run the main manuscript workflow"

check-env:
	@test -n "$$RKJCOLLAB" || (echo "RKJCOLLAB is not set"; exit 1)
	@test -d "$$RKJCOLLAB" || (echo "RKJCOLLAB does not point to an existing directory: $$RKJCOLLAB"; exit 1)

all: phenotypes data-summary lm-grs outcome-models manuscript

phenotypes: check-env
	$(R) 1_make_ia_t1d_phenos.R

data-summary: check-env
	$(R) -e 'source("config.R"); rmarkdown::render("2_analytical_data_summary.Rmd", output_dir = risk_scores_report_dir())'

lm-grs: check-env
	$(R) 3_lm_grs2_fdr_hla.R

outcome-models: check-env
	$(R) 4_ia_risk_log_reg.R

manuscript: check-env
	$(R) -e 'source("config.R"); rmarkdown::render("Manuscript_FamhistPRS.Rmd", output_dir = risk_scores_report_dir())'

time-varying: check-env
	$(R) -e 'source("config.R"); rmarkdown::render("time_varying_analysis.Rmd", output_dir = risk_scores_report_dir())'
