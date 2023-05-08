library(devtools)
library(rJava)
library(rlang)
options(unzip = "internal")
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor


devtools::install_github('dviraran/xCell')
library(xCell)

tpm_xCell<-xCellAnalysis(GSE_RPKM)

gsub(" ","_",row.names(tpm_xCell))->row.names(tpm_xCell)
gsub("\\+","",row.names(tpm_xCell))->row.names(tpm_xCell)
gsub("-","_",row.names(tpm_xCell))->row.names(tpm_xCell)
tpm_xCell[1:64,] -> tpm_xCell
as.data.frame(tpm_xCell) -> tpm_xCell
apply(tpm_xCell, 2, function(x){x/sum(x)}) -> yue

write.table(yue, "GSE106804_xcell.txt", quote = F, row.names = T,sep='\t')


library('e1071')
library('parallel')
# BiocManager::install('preprocessCore')
library('preprocessCore')
library(dplyr)
library(reshape2)

write.table(GSE_RPKM, file = "GSE93070_RPKM.txt", quote = F, sep = "\t")

source("cibersort_ann.R")

res_cibersort <- CIBERSORT("LM22.txt", "GSE106804.txt", perm=1000, QN=F);