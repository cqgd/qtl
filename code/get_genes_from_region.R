#Get genes from bed to granges region, base don https://www.biostars.org/p/312666/
#See also regions_aorund_hits_bed.R
library(data.table)
library(rtracklayer)
library(ensembldb)
library(EnsDb.Hsapiens.v79)
edb = EnsDb.Hsapiens.v79
# gr <- GRanges(seqnames = 16, IRanges(30768000, 30770000), strand = "*")
# genes(edb, filter = gr)
# 
# GA <- ensembldb::select(EnsDb.Hsapiens.v79, keys= C$gene_ensembl, keytype = "GENEID", columns = c("SYMBOL"))

setwd("~/work/qtl")

bed = import("resource/hits.bed") #update this with the actual COVID-19 hits

get_region_genes = function(bed_region, edb){
  gr = GRangesFilter(bed_region)
  region_genes = genes(edb, filter = gr)$gene_id
}

region_genes = data.table()

for(i in seq_along(bed)){
  SYMBOL = get_region_genes(bed[i],edb)
  HIT = rep(bed[i]$name, length(SYMBOL))
  data.table(SYMBOL,HIT)
  region_genes = rbind(region_genes,
                        data.table(SYMBOL,HIT))
}
print(region_genes)

#checking all are still there
#genes_interest = c("ENSG00000173578", "ENSG00000163820", "ENSG00000172215", "ENSG00000173585", "ENSG00000163818", "ENSG00000163817") ##replace with get_genes_from_region.R
#genes_interest%in%region_genes

fwrite(region_genes,"resource/region_genes.txt")

genes_interest = fread("resource/region_genes.txt")$SYMBOL
