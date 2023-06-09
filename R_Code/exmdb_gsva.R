setwd("I:/2021_UNFINSH/exmdb2022/GSE_data/new_data/data_exp/GSE104926_ESCC")

aggregate(GSE104926_ESOPHA_Limma_res,by=list(GSE104926_ESOPHA_Limma_res$Gene),FUN=mean,na.rm=TRUE) -> res

res$log2FoldChange = log2(res$esophagitis/res$control)
res$logp <- -log10(res$p.value)


write.table(res, file = "res_GSE106804_limma_ESOPHA.txt", quote = F,row.names = F, sep = "\t")

res <- res_GSE106804_limma

library(msigdbr)
library(dplyr)
library(data.table)
library(GSVA)
library(limma)
library(stringr)
library(ggplot2)

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
h <- msigdbr(species = "Homo sapiens", # 物种拉丁名
             category = "H") #此处以hallmark为例，你也可以选择MSigDB的其他注释

h <- dplyr::select(h, gs_name, gene_symbol) %>% #或entrez_gene
  as.data.frame %>% 
  split(., .$gs_name) %>% 
  lapply(., function(x)(x$gene_symbol)) #或entrez_gene

# 在每个geneset里面去掉重复的基因
gs <- lapply(h, unique)
# 接下来去掉那些在两个或更多个pathways里出现过的genes
count <- table(unlist(gs))
keep <- names(which(table(unlist(gs)) < 2))
gs <- lapply(gs, function(x) intersect(keep, x))
# 过滤之后，很多pathway一个gene都不剩了，去掉这些
gs <- gs[lapply(gs, length) > 0]
# 预览过滤后的结果
head(gs)
res[,c(2,3)] -> res_exp

rownames(res_exp) <- res$Group.1
gsva_es <- gsva(as.matrix(res_exp), gs)

as.data.frame(gsva_es) -> gsva_es

gsva_es$score <- gsva_es$ESCC/gsva_es$control

write.table(gsva_es, file = "gsva_es_ESCA.txt", quote = F,row.names = T, sep = "\t")

