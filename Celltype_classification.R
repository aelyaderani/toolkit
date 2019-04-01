library(readxl)
library(dplyr)
library(Matrix)
library(fBasics)
library(rpgm)
#L22272_Gene <- read.csv(file = "Desktop/untitled9/LRNA_01737_1_BR_FFNuclei_C2_X5SCR_L22272_HVTTJBBXX_10pcs_70percVariability_diffExpr.csv", header=TRUE, sep = ",")
#Stan_marker <- read.csv(file = "Desktop/Stan_avg_marker_final_c copy.csv", header=TRUE, sep = ",")
#out = "/Users/aelyaderani/Desktop/clusterclassifier"
#setwd(out)
genelist <- read_excel("/Users/aelyaderani/Desktop/celltypeclassifier/fileName.xlsx",col_names = FALSE)
#genelist2 <- read_excel("/Users/aelyaderani/Desktop/celltypeclassifier/testing.xlsx")
Stan_marker <- read_excel("/Users/aelyaderani/Desktop/celltypeclassifier/ref_list.xlsx")
filefirstpath <- "Desktop/new_seurat/"
filelastpath <- "Desktop/new_seurat_out/"
Rnum <- 1
csvext <- ".csv"  
xlsx <- ".xlsx"



while (Rnum <= nrow(genelist)){
filesecondpath <- genelist[Rnum,1]

fileName <- paste(filefirstpath,filesecondpath, sep = "")
fileNamelast <- paste(fileName,xlsx, sep = "")
L22272_Gene <- read_excel(fileNamelast)

RnumUNk <- 1
while(RnumUNk <= nrow(L22272_Gene)){
  L22272_Gene[RnumUNk,6] <- "(NA)"
  L22272_Gene[RnumUNk,7] <- "(NA)"
  RnumUNk = RnumUNk + 1
}
colnames(L22272_Gene)[colnames(L22272_Gene)=="..6"] <- "Cell Type"
colnames(L22272_Gene)[colnames(L22272_Gene)=="V6"] <- "Cell Type"
colnames(L22272_Gene)[colnames(L22272_Gene)=="V7"] <- "source of data"

sizeofSTAN <- (nrow(Stan_marker))
sizeofL2 <- (nrow(L22272_Gene))
Snum <- 1
while(Snum<=sizeofSTAN){ 
    if (tolower(Stan_marker[Snum,1])==tolower(L22272_Gene[1,1])){
          L22272_Gene[1,6] <- Stan_marker[Snum,2]}
    Snum = Snum+1}

Snum <- 1
Lnum <- 2
while(Lnum<=sizeofL2){ 
  Snum <- 1
  print(Lnum)
  while(Snum<=sizeofSTAN){
  
     # if (Stan_marker[Snum,1]==L22272_Gene[Lnum,1]){
      #  if (L22272_Gene[Lnum,2]==L22272_Gene[Lnum-1,2]){
       #   if(Stan_marker[Snum,2]==L22272_Gene[Lnum-1,6]){
        #    L22272_Gene[Lnum,6] <- Stan_marker[Snum,2]}
         # else{
          #  L22272_Gene[Lnum,6] <- Stan_marker[Snum,2]}
        #}
      #}
    if (Stan_marker[Snum,1]==L22272_Gene[Lnum,1]){
      L22272_Gene[Lnum,6] <- Stan_marker[Snum,2]
      L22272_Gene[Lnum,7] <- "Dr.Bakken, PMCID: PMC6306246"}
    
    Snum = Snum+1}
  Lnum = Lnum+1}
#savepath <- genelist2[Rnum,1]
savepath2 <- paste(filelastpath,filesecondpath, sep = "")

write.csv(L22272_Gene,file = paste(savepath2,".csv",sep=""), row.names=TRUE)
Rnum = Rnum + 1
}

