# Shared settings

config <- list(
  studies = c("teddy"),
  dr_filt = "no",  # set to "no", or "yes" to filter to only DR3/4, DR4/4, DR3/3 or DR4/X
  fdr_var = "fdr_4level",  # fdr_3level or fdr_4level
  fdr_ref = "Dad",  # Mom, Dad, Sib, None
  subset = "all"  # all, ctrls, female, male
)