library(data.table)

setwd("~/work/qtl/")

#get promising loci from sumstats

prune_by_p = function(x,P_threshold=5e-8){
  x = copy(x[P<=P_threshold])
  if(!nrow(x)>=1){return(NULL)}
  
  x[,P_copy:=P][,picked:=as.numeric(NA)]

  i=1
  now = x[which.min(P)]

  while(nrow(now)==1 & now$P<=P_threshold){
    x[which.min(P),picked:=i]
    x[,dist:=abs(BP-now$BP)]
    x[dist<=250e3,P:=Inf]
    now = x[which.min(P)]
    i=i+1
  }

  x[,P:=P_copy][,P_copy:=NULL][,dist:=NULL]
  return(x[is.finite(picked)])
}

X = fread("hgi/COVID19_HGI_B2_ALL_20200930.txt.gz_1.0E-5.txt")
relevant_cols = c("#CHR","POS","all_inv_var_meta_p","rsid")
new_names=c("CHR","BP","P","rsid")
X = X[,.SD,.SDcols=relevant_cols]
setnames(X,old=relevant_cols,new=new_names)

pruned = rbindlist(by(X,X$CHR,prune_by_p))

#add our old top hit by hand
hits = rbind(pruned,X[rsid=="rs34668658"],fill=TRUE)

setorder(hits, CHR, BP)

max_dist_from_locus = 1e6

bed = hits[,.(CHR,
              BP.START=BP-max_dist_from_locus,
              BP.END=BP+max_dist_from_locus,
              NAME=paste0(CHR,":",BP))]
fwrite(bed,"resource/hits.bed",col.names = FALSE,sep="\t")


#now run get_genes_from_region.R




# setwd("~/work/qtl")
# x = fread("gwas/ENSG00000163817_TISSUE-bulk/blmm_output.bgen.stats.gz") #random GWAS to see te
# 1e6-diff(x[,range(BP)]) #confirming we have 1MBp of SNPs around the bait via Luke, so now let's see 1Mbp on either side of that
# bed = x[SNP=="rs34668658",.(CHR,BP.START=BP-4e6,BP.END=BP+4e6)]
# fwrite(bed,"resource/hits.bed",col.names = FALSE,sep="\t")


