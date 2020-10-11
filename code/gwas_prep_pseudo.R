setwd("~/work/qtl")


library(data.table)

#Load sample meta-data to see if sample names match
M = fread("data/Final_COMBAT_basic_clinical_data_freeze_170820.txt")

#Format pseudo single-cell expression data

##BROAD 
B = readRDS("data/pseudobulk_counts_broad.rds")
genes_interest = fread("resource/region_genes.txt")$SYMBOL
#genes_interest = c("ENSG00000173578", "ENSG00000163820", "ENSG00000172215", "ENSG00000173585", "ENSG00000163818", "ENSG00000163817")
B = B[dimnames(B)[[1]]%in%genes_interest,,]

slices = list()
tissues = dimnames(B)[[2]]
t=1
for(t in 1:length(tissues)){
  slice = data.table(t(B[,t,]),keep.rownames = FALSE)
  colnames(slice)<- paste0(colnames(slice),"_TISSUE-broad-",tissues[t])
  slices[[t]] = slice
}

z = do.call(cbind,slices)
z[,pseudo_sample:=dimnames(B)[[3]]]
setcolorder(z,"pseudo_sample")
z[,pseudo_sample:=gsub("\\.","-",pseudo_sample)]

all(z$pseudo_sample %in% M$scRNASeq_sample_ID)
setnames(z,old="pseudo_sample",new="scRNASeq_sample_ID")

fwrite(z,"data/pseudobulk_counts_broad.csv")


##NARROW
B = readRDS("data/pseudobulk_counts_narrow.rds")
genes_interest = fread("resource/region_genes.txt")$SYMBOL
#genes_interest = c("ENSG00000173578", "ENSG00000163820", "ENSG00000172215", "ENSG00000173585", "ENSG00000163818", "ENSG00000163817")
B = B[dimnames(B)[[1]]%in%genes_interest,,]

slices = list()
tissues = dimnames(B)[[2]]
t=1
for(t in 1:length(tissues)){
  slice = data.table(t(B[,t,]),keep.rownames = FALSE)
  colnames(slice)<- paste0(colnames(slice),"_TISSUE-narrow-",tissues[t])
  slices[[t]] = slice
}

z = do.call(cbind,slices)
z[,pseudo_sample:=dimnames(B)[[3]]]
setcolorder(z,"pseudo_sample")
z[,pseudo_sample:=gsub("\\.","-",pseudo_sample)]

all(z$pseudo_sample %in% M$scRNASeq_sample_ID)
setnames(z,old="pseudo_sample",new="scRNASeq_sample_ID")

fwrite(z,"data/pseudobulk_counts_narrow.csv")




#FWRITE OUT

#DO THE SAME FOR NARROW

#Merge this into the previous pheno_cov, then put the relevant columns in a file and load CC_ID from there when submitting job
#Make it two columns for convenience, number and name