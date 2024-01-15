library(devtools)
library(rJava)
library(rlang)
options(unzip = "internal")
Sys.setenv(LANGUAGE = "en") #Display English error messages
options(stringsAsFactors = FALSE) #Prohibit the conversion of chr to factor


devtools::install_github('dviraran/xCell')
library(xCell)

tpm_xCell<-xCellAnalysis(GSE_RPKM)

gsub(" ","_",row.names(tpm_xCell))->row.names(tpm_xCell)
gsub("\\+","",row.names(tpm_xCell))->row.names(tpm_xCell)
gsub("-","_",row.names(tpm_xCell))->row.names(tpm_xCell) # row name conversion
tpm_xCell[1:64,] -> tpm_xCell
as.data.frame(tpm_xCell) -> tpm_xCell
apply(tpm_xCell, 2, function(x){x/sum(x)}) -> percentage

library('e1071')
library('parallel')
# BiocManager::install('preprocessCore')
library('preprocessCore')
library(dplyr)
library(reshape2)

source("cibersort_ann.R") #Loading the cibersort program

res_cibersort <- CIBERSORT("LM22.txt", "GSE106804.txt", perm=1000, QN=F);


