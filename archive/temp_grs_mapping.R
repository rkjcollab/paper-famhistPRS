

# SNP list TEDDY MatProt GRS2 with 61 SNPs

setwd(Sys.getenv("RKJCOLLAB"))

library(tidyverse)
library(readxl)

# From TEDDY MatProt scoring used in poster
matprot_list <- read_tsv(
  "Maternal_Protection/data/teddy/genetics/risk_scores/grs2/scoring_1000g/grs2_score_input_with_info.txt",
  col_types = "ccccccccc")

# From previous repo
grs2_1000g_list <- read_xlsx(
  "/Users/slacksa/repos/explore_matprot/hla-prs-toolkit/Snplists/T1D_GRS67/T1D_GRS67_1000G_nopalin_pos_hg19.xlsx")

# From new PRSedm (GRS2x) repo
# Think need to map from 1000G -> TM R2 -> TM R3 used "replaced" column
grs2x_map <- read_xlsx(
  "/Users/slacksa/repos/explore_matprot/genetics/risk_scores/PRSedm/snplists/grs2_version_snplists.xlsx",
  sheet = "GRS2X_TOPMED_R3")
grs2_tm_map <- read_xlsx(
  "/Users/slacksa/repos/explore_matprot/genetics/risk_scores/PRSedm/snplists/grs2_version_snplists.xlsx",
  sheet = "GRS2_TOPMED_R2")
grs2_1000g_map <- read_xlsx(
  "/Users/slacksa/repos/explore_matprot/genetics/risk_scores/PRSedm/snplists/grs2_version_snplists.xlsx",
  sheet = "GRS2_HRC_1000G")

# Check new & old repo 1000G map


# Check overlap
length(intersect(grs2_1000g_list$RSID, matprot_input$RSID))  # 61

# Get missing SNPs
grs2_not_used <- grs2_1000g_list %>%
  dplyr::filter(!RSID %in% matprot_input$RSID)
grs2_not_used_snps <- grs2_not_used$RSID
grs2_not_used_snps
# "rs17843689"  "rs111485156" "rs1281935"   "rs3129727"   "rs9271346"   "rs1281934" 

# See if missing SNPs in GRS2x as is?
grs2x_map %>%
  dplyr::filter(RSID %in% grs2_not_used_snps)
# HLA-DQ73  rs1281935 6:32616043    G     T       0.0373 +      NA    T             NA 
# One as is

# See if missing SNPs in "replaced" columns?
grs2_tm_map %>%
  dplyr::filter(REPLACED %in% grs2_not_used_snps)
# HLA-DQ62    rs9273342 6:32655666    T     G       0.864  +      NA    T             rs17843689
# HLA Class 2 rs9271347 6:32615766    A     G       0.870  +      1.69  G             rs9271346 
# HLA Class 2 rs1281943 6:32630313    T     C       0.0191 +      0.9   C             rs1281934 
# Three replaced in TM GRS2

grs2x_map %>%
  dplyr::filter(REPLACED %in% grs2_not_used_snps)
# 0

# 2 not tracked by replaced column are DQ-##, so know which are
# DQ75 & DQ53
grs2x_map %>%
  dplyr::filter(COMPONENT == "HLA-DQ75" | COMPONENT == "HLA-DQ53")
# HLA-DQ53  rs1794268 6:32706288    T     C       0.0230 +      NA    C             rs1794265
# HLA-DQ75  rs9469200 6:32635435    T     C       0.130  +      NA    C             NA   

# Note that both replacements are listed in the previous repo list for 1000G,
# but not the updated repo - weird?

