library(dplyr)
library(readr)
library(tibble)

prep_teddy <- function(specs) {
  
  # Load
  pheno_raw <- read_tsv(specs$pheno_raw_path)
  fid <- read_delim(specs$fid_path, col_names = F)
  twin_list <- read_tsv(specs$twin_list_path)
  grs2 <- read_delim(specs$grs2_path)
  pc <- readRDS(specs$pcair_path)
  
  # Add FID to pheno data
  # Assume input given is .bim, so first two cols are FID IID
  fid <- fid %>%
    dplyr::rename(
      FID = X1,
      IID = X2) %>%
    dplyr::select(FID, IID)
  pheno <- pheno_raw %>%
    dplyr::left_join(
      fid,
      by = c("IID"))
  
  # Study-specific prep
  grs2$ID <- gsub("^.+_", "", grs2$IID)
  
  # Because of new delivery structure, want to add in follow up time for
  # controls to outcome variables
  pheno <- pheno %>%
    mutate(
      fupIA = ifelse(is.na(fupIA), fup, fupIA),
      fupT1D = ifelse(is.na(fupT1D), fup, fupT1D),
      fupT1D_prog = ifelse(
        IA == 1, ifelse(
          T1D == 1, fupT1D_prog, fup - fupIA), NA))
  
  # Merge and filter
  df <- inner_join(
    grs2,
    pheno %>% mutate(IID = as.character(IID)),
    by = c("ID" = "IID"))
  
  # Remove tiwns (identified genetically in
  # immuno_t1d/genetics/ancestry_estimation/TwinFinder.qmd)
  df <- df %>% filter(!ID %in% twin_list$IID)
  
  # Add genetic PCs (in EUR only)
  pc_df <- as.data.frame(pc$vectors) %>%
    rownames_to_column("id")
  colnames(pc_df) <- gsub("^V", "PC", colnames(pc_df))
  
  df <- inner_join(df, pc_df[, 1:11], by = c("ID" = "id"))
  
  # Filter to remove individuals with more than one FDR
  df <- df %>%
    filter(!(family_mem_screen %in% 1:4))
  
  # HLA mapping
  df <- df %>%
    mutate(
      HLAGRP = case_when(
        hla_category == 1 ~ "DR3/4",
        hla_category %in% c(2,3,7) ~ "DR4/4",
        hla_category %in% c(4,5,6,8) ~ "DR4/X",
        hla_category == 9 ~ "DR3/3",
        hla_category == 10 ~ "DR3/X"))
  
  # Standardize column names across versions
  df <- df %>%
    rename(
      fdr_4level = any_of(c("fdr_4level", "fdr_4level_screen")),
      fdr_3level = any_of(c("fdr_3level", "fdr_3level_screen")),
      sex = Sex,
      GRS2x = `t1dgrs2-luckett25_total`
    ) %>%
    mutate(
      sex = ifelse(sex == "Male", 1, 0))
  
  # Final columns
  cols_keep <- c(
    "FID", "ID", "dr34", "sex", "IA", "fupIA", "T1D", "fupT1D", "fupT1D_prog",
    "mAA_at_IA", "GRS2x", "Non_HLA", "fdr_3level", "fdr_4level", "cc",
    "PC1", "PC2", "PC3", "PC4", "PC5", "HLAGRP"
  )
  
  df[, intersect(cols_keep, names(df))]
}
