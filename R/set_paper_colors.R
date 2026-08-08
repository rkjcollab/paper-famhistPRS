# SDS 20251125

# Set colors for MatProt PRS paper.

# Setup ------------------------------------------------------------------------

library(colorspace)
library(khroma)

# Set colors ---------------------------------------------------------------

# Using "muted" by Paul Tol for colorblind friendly palette
# https://cran.r-project.org/web/packages/khroma/vignettes/tol.html

muted <- color("muted")

# FDR group colors
color_fdr_none <- unname(muted(9)[3])
color_fdr_dad <- unname(muted(9)[5])
color_fdr_mom <- unname(muted(9)[1])
color_fdr_sib <- unname(muted(9)[4])

colors_fdr <- c(
  "None" = as.character(color_fdr_none),
  "Dad" = as.character(color_fdr_dad),
  "Mom" = as.character(color_fdr_mom),
  "Sib" = as.character(color_fdr_sib))

# Darker FDR none color for graphical abstract only
color_fdr_none_abs <- darken(color_fdr_none, amount = 0.2)

# Case/control colors
color_case <- "#DDDDDD"  # grey included with muted
color_control <- "white"

# barplot(rep(1, length(colors_fdr)), col = colors_fdr, border = "white", space = 0,
#         names.arg = names(colors_fdr))
