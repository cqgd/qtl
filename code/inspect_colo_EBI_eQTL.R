library(data.table)
library(ggplot2)
library(cowplot)
library(stringr)
library(EnsDb.Hsapiens.v79)
theme_set(theme_cowplot())
setwd("~/work/qtl")

#Load coloc results
files = list.files("coloc_EBI_eQTL",pattern = "\\.coloc\\.RDS$")
coloc_paths = paste0("coloc_EBI_eQTL/",files)
C = lapply(coloc_paths,readRDS)

C = data.table(t(sapply(C,function(coloc_result){coloc_result$summary})))
C[,coloc:=files]

#Formatting file names
C[,gwas:=coloc]
C[,gene_ensembl:=tstrsplit(coloc,"_")[[3]]]
C[,coloc:=gsub("\\.chr3p21\\.tsv\\.coloc\\.RDS","",coloc)]
C[,coloc:=tstrsplit(coloc,"_FROM-")[[2]]]
src_split = "(20..)|(V8)|(CEDAR)|(GENCORD)|(BLUEPRINT)|(TwinsUK)|(GEUVADIS)|(HipSci)|(ROSMAP)|(BrainSeq)"
src_split_short = "(CEDAR)|(GENCORD)|(BLUEPRINT)|(TwinsUK)|(GEUVADIS)|(HipSci)|(ROSMAP)|(BrainSeq)"
C[,c("source","tissue_name"):=tstrsplit(coloc,src_split)]
C[nchar(source)==0,source:=str_match(coloc,src_split_short)[,1]]
C[,tissue_name:=trimws(gsub("_"," ",tissue_name))]
C[,source:=trimws(gsub("_"," ",source))]
C[,coloc:=NULL]

#Annotate genes
GA <- ensembldb::select(EnsDb.Hsapiens.v79, keys= C$gene_ensembl, keytype = "GENEID", columns = c("SYMBOL"))
setnames(GA,c("gene_ensembl","gene_desc"))
C = merge(C,GA,by="gene_ensembl",all.x=TRUE)

#Output
setcolorder(C,c("gwas","source","tissue_name","gene_desc","gene_ensembl"))
setorder(C,-PP.H4.abf)
fwrite(C,"tables/coloc_EBI_eQTL_results.csv")

#Merge with our own results

C_own = fread("tables/qtl_overview_with_coloc.csv")
C_own

C_both = rbind(C_own,C,fill=TRUE)
setcolorder(C_both,c("gwas","source","granularity"))
setorder(C_both,-PP.H4.abf)
fwrite(C_both[,.SD,.SDcols=setdiff(colnames(C_both),"tissue")],"tables/coloc_own_and_EBI_eQTL_results.csv")
