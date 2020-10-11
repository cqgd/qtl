library(data.table)
library(ggplot2)
library(cowplot)
library(coloc)
theme_set(theme_cowplot())
setwd("~/work/qtl")

eqtl_dataset = function(path){
  SS = fread(path)
  
  N_hardcoded=median(SS$an,na.rm = TRUE)
  
  SS = SS[,.(snp=rsid, pvalues=pvalue, MAF=maf)]
  dataset = list(type="quant", #Use this with P_LINREG?
                   pvalues=SS$pvalues,  
                   MAF=SS$MAF,
                   snp=SS$snp,
                   sdY=1,
                   N=N_hardcoded)
  return(dataset)
}

gwas_dataset = function(path){
  SS = fread(path)
  SS = SS[nchar(rsid)>=1]
  s = 3199/median(SS$all_meta_sample_N)
  SS = SS[,.(snp=rsid, pvalues=all_inv_var_meta_p, beta=all_inv_var_meta_beta, varbeta=all_inv_var_meta_sebeta^2)] #Should MAF be min(MAF,1-MAF)?
  dataset = list(type="cc",
                 s=s,
                 #pvalues=SS$pvalues, 
                 beta=SS$beta, 
                 varbeta=SS$varbeta,
                 snp=SS$snp)
  return(dataset)
}


gwas = gwas_dataset("hgi_gwas/COVID19_HGI_ANA_B2_V2_20200701.b37.chr3.txt")

gwases = list.files("EBI_eQTL/by_gene/",pattern="\\.tsv$")
paths = paste0(getwd(),"/EBI_eQTL/by_gene/",gwases)
exists = file.exists(paths)
gwases = gwases[exists]
paths = paths[exists]

for(e in 1:length(paths)){
  print(paste0(e,"/",length(paths)))
  eqtl = eqtl_dataset(paths[e])
  coloc_result = coloc.abf(gwas,eqtl)
  print(coloc_result)
  saveRDS(coloc_result,file=paste0("coloc_EBI_eQTL/EBI_eQTL_",gwases[e],".coloc.RDS"))
}
