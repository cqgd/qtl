setwd("~/work/qtl")

library(data.table)
S = fread("data/COMBAT_topmed_filtered_chr3p21.sample")
M = fread("data/Final_COMBAT_basic_clinical_data_freeze_170820.txt")
E = fread("data/Logcpm_143_23063.txt")
W = fread("data/COMBAT_100gpca_ancestryassignments.txt")
PC = fread("data/PCs.csv")


#Transpose expression data to have 1 column per gene
Et = transpose(E,keep.names="names")
Et[1,1]="RNASeq_sample_ID"
colnames(Et) <- unlist(Et[1,])
Et = Et[2:nrow(Et)]
Et[1:10,1:10]

#Add sample meta-data
X = merge(Et,
          M[,.(COMBAT_ID,RNASeq_sample_ID,sample_priority,age=Age,sex=Sex, Ethnicity)],
          by="RNASeq_sample_ID",
          all.x=TRUE)

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


#Add home-made PCs
# colnames(PC) = gsub("_","",colnames(PC))
# colnames(PC)[1] <-"COMBAT_ID"
# PC = PC[,-2]
# X = merge(X,
#           PC,
#           by="COMBAT_ID",
#           all.x=TRUE)


#Subset by inferred ethnicity in EUR/SAS
X = X[pca_pop_gp%in%c("EUR","SAS")]

#Subset to most severe expression sample for each individual
X = X[sample_priority==1]

#Subset to only genes of interest
#XCR1, FYCO1, CXCR6, CCR9, LZTFL1 and SLC6A20

extra_cols = colnames(X)[!grepl("^ENSG",colnames(X))]
expr_cols = c("ENSG00000173578", "ENSG00000163820", "ENSG00000172215", "ENSG00000173585", "ENSG00000163818", "ENSG00000163817")
X = X[,.SD,.SDcols=c(extra_cols,expr_cols)]

#Convert expression values to numeric, you can do this by removing V1 from E before transposing also.
#Then quantile-normalize
quantile.normalize = function(v){
  qnorm((frank(v,ties.method="average",na.last="keep")-0.5)/length(v[!is.na(v)]))
}

for(col in expr_cols){
  X[,(col):=quantile.normalize(as.numeric(get(col)))]
}

#Rename columns as per bolt requirements
X[,FID:=COMBAT_ID]
X[,IID:=COMBAT_ID]
setcolorder(X,c("FID","IID"))

#Output
fwrite(X,"pheno_cov.csv",sep="\t",col.names=T,row.names=F,quote=F,na="NA")


# library(ggplot2)
# library(cowplot)
# theme_set(theme_cowplot())
# ggplot(X,aes(x=PC_1,y=PC_2)) + geom_point(aes(col=Ethnicity)) + background_grid() + ggtitle("Clustering of ethnicity by first two genome-wide PCs")
