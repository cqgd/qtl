library(data.table)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
setwd("~/work/qtl")

eqtl_dataset = function(path){
  SS = fread(path)
  
  N=70 
  
  if("P_BOLT_LMM"%in%colnames(SS)){
    SS = SS[,.(snp=SNP, pvalues=P_BOLT_LMM, beta=BETA, varbeta=SE^2)] #If using P_LINEREG, maybe just use p-values only below?
    dataset = list(type="quant", #When run with -forceLMM and you have P_BOLT_LMM
                   beta=SS$beta, 
                   varbeta=SS$varbeta,
                   snp=SS$snp,
                   sdY=1)
  }else{
    SS = SS[,.(snp=SNP, pvalues=P_LINREG, MAF=A1FREQ)]
    dataset = list(type="quant", #Use this with P_LINREG?
                   pvalues=SS$pvalues,  
                   MAF=SS$MAF,
                   snp=SS$snp,
                   sdY=1,
                   N=70)
  }
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

gwases = list.files("gwas")
paths = paste0(getwd(),"/gwas/",gwases,"/blmm_output.bgen.stats.gz")
exists = file.exists(paths)
gwases = gwases[exists]
paths = paths[exists]

for(e in 1:length(paths)){
  print(paste0(e,"/",length(paths)))
  eqtl = eqtl_dataset(paths[e])
  coloc_result = coloc.abf(gwas,eqtl)
  print(coloc_result)
  saveRDS(coloc_result,file=paste0("coloc/",gwases[e],".coloc.RDS"))
}
