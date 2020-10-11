fam = fread("/home/cq/work/qtl/COMBAT_genotyped_qcd.fam")
fam[,V1:=V2] #Replace 1...N by COMBAT_ID
fwrite(fam,"/home/cq/work/qtl/COMBAT_genotyped_qcd_fixed.fam",col.names = FALSE,sep=" ") 

S = fread("/home/cq/work/qtl/data/COMBAT_topmed_filtered_chr3p21.sample")

#First IDs all match (after the 0 0 entry in the .sample)
all(S$ID_1[2:nrow(S)]==fam$V1)

#Concatenate IDs
fam[,IDs:=paste(V1,V2,sep="_")]
S[,IDs:=paste(ID_1,ID_2,sep="_")]

#Both IDs match, actually
all(S$IDs[2:nrow(S)]==fam$IDs)

#Could be something with the delimeter
