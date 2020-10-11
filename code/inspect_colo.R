library(data.table)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
setwd("~/work/qtl")

#Load coloc results
files = list.files("coloc",pattern = "\\.coloc\\.RDS$")
coloc_paths = paste0("coloc/",files)
C = lapply(coloc_paths,readRDS)

C = data.table(t(sapply(C,function(coloc_result){coloc_result$summary})))
C[,coloc:=files]

setorder(C,-PP.H4.abf)

#Format tissue and ensembl ID information
setnames(C,old="coloc",new="tissue")
C[,tissue:=gsub("\\.coloc\\.RDS$","",tissue)]
C[,gwas:=tissue]
C[,c("gene_ensembl","tissue"):=tstrsplit(tissue,"_TISSUE-")]
C[,tissue:=gsub("_slash_","/",tissue)]
C[,tissue:=gsub("_"," ",tissue)]
C[,tissue:=gsub("narrow-","[higher granularity clustering]  ",tissue)]
C[,tissue:=gsub("broad-","[lower granularity clustering]  ",tissue)]
C[,tissue:=gsub("bulk-","[bulk] ",tissue)]
C[,c("granularity","tissue_name"):=tstrsplit(tissue,"  ")]
C = C[,.SD,.SDcols=setdiff(colnames(C),"tissue")]


#Add gene names
gene_desc = c("XCR1", "FYCO1", "CXCR6", "CCR9", "LZTFL1","SLC6A20")
gene_ensembl = c("ENSG00000173578", "ENSG00000163820", "ENSG00000172215", "ENSG00000173585", "ENSG00000163818", "ENSG00000163817")
GI = data.table(gene_desc, gene_ensembl)
C = merge(C,
           GI,
           by="gene_ensembl",
           all.x=TRUE, sort=FALSE)


#Output
setcolorder(C,c("gwas","granularity","tissue_name","gene_desc","gene_ensembl"))
fwrite(C,"tables/coloc_results.csv")
