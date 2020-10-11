library(data.table)
#Convert Luke's tabix eQTL datasets from EBI_eQTL to one file per gene/tissue/study.

setwd("~/work/qtl/EBI_eQTL")
paths = list.files(".",pattern="\\.tsv$")
X = lapply(paths,fread)
lapply(X,colnames)

#Add dataset identifiers, save by gene

save_eqtl_file = function(x_one_gene){
  out_path = paste0("by_gene/",x_one_gene$eqtl_file[1])
  fwrite(x_one_gene,out_path)
  print(out_path)
}

for(i in 1:length(X)){
  print(paste0(i,"/",length(X)))
  X[[i]][,tabix_file:=paths[i]]
  X[[i]][,eqtl_file:=paste0(gene_id,"_FROM-",tabix_file)]
  z = copy(X[[i]])
  by(z,factor(z$gene_id),save_eqtl_file)
}
