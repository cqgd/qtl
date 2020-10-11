x = fread("/home/cq/work/qtl/GSE86189_LG.po.all.txt")
x[,c("chr1","start1","end1"):=tstrsplit(frag1,"\\.")]
x[,c("chr2","start2","end2"):=tstrsplit(frag2,"\\.")]
x[chr1!=chr2]
