library(data.table)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
setwd("~/work/qtl")

gwases = list.files("gwas")
paths = paste0(getwd(),"/gwas/",gwases,"/blmm_output.bgen.stats.gz")
exists = file.exists(paths)
gwases = gwases[exists]
paths = paths[exists]

Ss = lapply(paths,fread)

add_gwas_name = function(S,gwas){
  S[,gwas:=gwas]
}

Ss = mapply(add_gwas_name,Ss,gwases,SIMPLIFY = FALSE)


Ss = rbindlist(Ss)

Ss = Ss[,c("gene_ensembl","tissue"):=tstrsplit(gwas,"_TISSUE-")]
Ss[,tissue:=gsub("_slash_","/",tissue)]
Ss[,tissue:=gsub("_"," ",tissue)]
Ss[,tissue:=gsub("narrow-","[higher granularity clustering]  ",tissue)]
Ss[,tissue:=gsub("broad-","[lower granularity clustering]  ",tissue)]
Ss[,tissue:=gsub("bulk-","[bulk] ",tissue)]

gene_desc = c("XCR1", "FYCO1", "CXCR6", "CCR9", "LZTFL1","SLC6A20")
gene_ensembl = c("ENSG00000173578", "ENSG00000163820", "ENSG00000172215", "ENSG00000173585", "ENSG00000163818", "ENSG00000163817")
GI = data.table(gene_desc, gene_ensembl)


Ss = merge(Ss,
      GI,
      by="gene_ensembl",
      all.x=TRUE, sort=FALSE)

fwrite(Ss,"qtls_v1.csv")

Ominp = Ss[,.(min_p=min(P_LINREG), min_SNP=SNP[which.min(P_LINREG)]),by=gwas]
Ophit = Ss[SNP=="rs34668658",.(rs34668658_p=P_LINREG),by=gwas]
O = merge(Ominp, Ophit, by="gwas",all=TRUE)
O = O[,c("gene_ensembl","tissue"):=tstrsplit(gwas,"_TISSUE-")]
O[,tissue:=gsub("_slash_","/",tissue)]
O[,tissue:=gsub("_"," ",tissue)]
O[,tissue:=gsub("narrow-","[higher granularity clustering]  ",tissue)]
O[,tissue:=gsub("broad-","[lower granularity clustering]  ",tissue)]
O[,tissue:=gsub("bulk-","[bulk] ",tissue)]
O[,c("granularity","tissue_name"):=tstrsplit(tissue,"  ")]

O = merge(O,
          GI,
          by="gene_ensembl",
          all.x=TRUE, sort=FALSE)

O = O[,.(gwas,tissue, granularity,tissue_name,gene_desc,gene_ensembl,rs34668658_p,min_p,min_SNP)]
setorder(O,min_p)
O
fwrite(O,"qtl_pseudo_overview.csv")

mhplot = function(S){
  hit_p = formatC(S[SNP=="rs34668658",P_LINREG], format = "e", digits = 2)
  
  thres = S[,.(gene_has_hit=any(P_LINREG<5e-8)),by=gene_desc]
  thres = thres[gene_has_hit==TRUE]
  
  p.mh = ggplot(S,aes(x=BP,y=-log10(P_LINREG))) + 
    geom_point() + 
    geom_point(data=S[SNP=="rs34668658"],color="red",alpha=0.4,cex=8) +
    geom_point(data=S[SNP=="rs34668658"],color="red",cex=2)
    
  if(nrow(thres)>=1){
    p.mh = p.mh + geom_hline(data=thres,aes(yintercept =-log10(5e-8)),col="orange",alpha=0.7,size=2)
  }
  
  p.mh = p.mh +  ggtitle(paste0("",unique(S$tissue))) + 
    facet_wrap(~gene_desc,ncol=3,scales = "free",nrow=2,drop=FALSE) +
    xlab("CHR3 BP") + ylab("-log10(p)") +
    background_grid() +
    scale_y_continuous(breaks=seq(0,15,by=1)) + 
    theme(legend.position = "none")
  
  print(p.mh)
  plot_path = paste0("plots/",gsub("\\/","_slash_",unique(S$tissue)),".pdf")
  print(plot_path)
  ggsave2(filename=plot_path,
          plot=p.mh,
          width=16,
          height=7,
          units="in")
  return(p.mh)
}



Ss[,gene_desc:=as.factor(gene_desc)]
tissues = unique(Ss$tissue)
for(t in 1:length(tissues)){
  print(t)
  S = Ss[tissue==tissues[t]]
  mhplot(S)
}


