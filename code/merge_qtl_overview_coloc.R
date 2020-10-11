library(data.table)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
setwd("~/work/qtl")

C = fread("tables/coloc_results.csv")
O=fread("qtl_pseudo_overview.csv")
unique_to_C_cols = colnames(C)[!colnames(C)%in%colnames(O)]
M = merge(O,
      C[,.SD,.SDcols=c(unique_to_C_cols, "gwas")],
      by="gwas",
      all.x=TRUE)
setorder(M,-PP.H4.abf)
M
fwrite(M,"tables/qtl_overview_with_coloc.csv")
