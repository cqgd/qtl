setwd("~/work/qtl")

library(data.table)
S = fread("data/COMBAT_topmed_filtered_chr3p21.sample")
M = fread("data/Final_COMBAT_basic_clinical_data_freeze_170820.txt")
E = fread("data/Logcpm_143_23063.txt")
W = fread("data/COMBAT_100gpca_ancestryassignments.txt")
PC = fread("data/PCs.csv")
Epb = fread("data/pseudobulk_counts_broad.csv")
Epn = fread("data/pseudobulk_counts_narrow.csv")


#Transpose expression data to have 1 column per gene
Et = transpose(E,keep.names="names")
Et[1,1]="RNASeq_sample_ID"
colnames(Et) <- unlist(Et[1,])
Et = Et[2:nrow(Et)]
Et[1:10,1:10]
colnames(Et)[2:ncol(Et)] <- paste0(colnames(Et)[2:ncol(Et)],"_TISSUE-bulk")

#Take sample meta-data and add expression data

X = merge(M[,.(COMBAT_ID,RNASeq_sample_ID,scRNASeq_sample_ID,sample_priority,age=Age,sex=Sex)],
          Et,
          by="RNASeq_sample_ID",
          all.x = TRUE)

X = merge(X,
          Epb,
          by="scRNASeq_sample_ID",
          all.x = TRUE)

X = merge(X,
          Epn,
          by="scRNASeq_sample_ID",
          all.x = TRUE)

X[,.(RNASeq_sample_ID,sample_priority)]

#Normalize age and sex
X[,sex.raw:=sex] #works
X[,age.raw:=age]

X[,sex:=as.numeric(as.factor(sex))]
X[,sexage := sex*age]
X[,sexage := sexage - mean(sexage,na.rm=T)]
X[,age := age - mean(age,na.rm=T)]
X[,sex := sex - mean(sex,na.rm=T)]
X[,agesq := age^2]
X[,sexagesq := sexage^2]

#Add inferred ethnicity and PCs
X = merge(X,
          W,
          by="COMBAT_ID",
          all.x=TRUE)

#Subset by inferred ethnicity in EUR/SAS
X = X[pca_pop_gp%in%c("EUR","SAS")]

#Subset to most severe expression sample for each individual
X = X[sample_priority==1]

#Subset to only genes of interest
#XCR1, FYCO1, CXCR6, CCR9, LZTFL1 and SLC6A20

extra_cols = colnames(X)[!grepl("ENSG",colnames(X))]


genes_interest = fread("resource/region_genes.txt")$SYMBOL
#genes_interest = c("ENSG00000173578", "ENSG00000163820", "ENSG00000172215", "ENSG00000173585", "ENSG00000163818", "ENSG00000163817") ##replace with get_genes_from_region.R
expr_cols = grep(paste(genes_interest,collapse="|"),colnames(X),value=TRUE)
X = X[,.SD,.SDcols=c(extra_cols,expr_cols)]


#Convert expression values to numeric, you can do this by removing V1 from E before transposing also.
#Then quantile-normalize
quantile.normalize = function(v){
  qnorm((frank(v,ties.method="average",na.last="keep")-0.5)/length(v[!is.na(v)]))
}

# for(i in seq_along(expr_cols)){
#   cat(paste0("\n",i,"/",length(expr_cols)))
#   col = expr_cols[i]
#   X[,(col):=quantile.normalize(as.numeric(get(col)))]
# }

#much faster:
system.time({ 
  
  M = as.matrix(X[,.SD,.SDcols=expr_cols])
  Mqn =  apply(M,2,quantile.normalize)
  X[,(expr_cols):=as.data.table(Mqn)]       
  
})

#Rename columns as per bolt requirements
X[,FID:=COMBAT_ID]
X[,IID:=COMBAT_ID]
setcolorder(X,c("FID","IID"))

#Output
safe_expr_cols = gsub("\\/","_slash_",gsub(" ","_",expr_cols))
setnames(X,old=expr_cols,new=safe_expr_cols)

#figure out which genes have to run with which set of loci/bgen (determined by the hit they're near)
region_genes = fread("resource/region_genes.txt")
runs = data.table(PHENOTYPE=safe_expr_cols)
runs[,SYMBOL:=tstrsplit(safe_expr_cols,"_")[[1]]]
##some expression cols will have to be run multiple times with different bgens (because the genes are relevant to hit A but also hit B), hence cartesian merge.
runs = merge(runs,region_genes,by="SYMBOL",all.x = TRUE,allow.cartesian = TRUE) 

fwrite(X,"pheno_cov/pheno_cov.csv",sep="\t",col.names=T,row.names=F,quote=F,na="NA")
fwrite(runs[,.(HIT, PHENOTYPE)],"pheno_cov/runs.csv",col.names = FALSE)
# library(ggplot2)
# library(cowplot)
# theme_set(theme_cowplot())
# ggplot(X,aes(x=PC_1,y=PC_2)) + geom_point(aes(col=Ethnicity)) + background_grid() + ggtitle("Clustering of ethnicity by first two genome-wide PCs")
