
limma_file <- function(x,y)
{
  exp<-x
  exp<-as.matrix(exp)
  a<-y
  y[,1]<-colnames(exp)
  y<-as.data.frame(y)
  library(limma)##???ذ?
  eset=exp 
  targets<-y

  targets$V1=c(paste0("F",c(targets$V1),collapse = NULL,sep=""))
  colnames(targets)=c("FileName","Target")
  lev<-unique(targets$Target)#unique
  f <- factor(targets$Target, levels=lev) 
  design <- model.matrix(~0+f) # divide into groups
  colnames(design) <- c("control","case") #levels


  cont.wt <- makeContrasts("case-control",
                           levels=design) 
  fit <- lmFit(eset, design)
  fit2 <- contrasts.fit(fit, cont.wt) 
  fit2 <- eBayes(fit2) 
  tT=topTable(fit2, adjust="BH",sort.by="logFC",n=Inf)
  tT = subset(tT, select=c("adj.P.Val","P.Value","logFC"))
  -log2(tT$P.Value) -> tT$Log2p
  colnames(tT)=c("FDR","P.Value","log2FC","Log2(p)")
  tT$name <- rownames(tT)
  tT
}  # Methods for calculating Limma


list.files("D:list.files("F:/2021unf/exmdb2022/GSE_data/Disease") -> Single_GSE
limma_res <- list()
data_id <- 1 #Calculate from the first data
for (i in Single_GSE) {
  data_id = data_id + 1
  paste("F:/22021unf/exmdb2022/GSE_data/Disease/",i,"/",i,".txt",sep = "") -> Data_path
  paste("F:/2021unf/exmdb2022/GSE_data/Disease/",i,"/","Group.txt",sep = "") -> Group_path
  GSE_Data <- read.delim(Data_path, row.names=1)
  Group <- read.delim(Group_path)
  paste("limma_",i,".txt",sep = "") -> outfile_name
  GSE_Data[is.na(GSE_Data)] <-  0
  limma_file(GSE_Data,Group) -> b
  colnames(Group)<-c("id","Group")
  data.frame(cbind(rowMeans(GSE_Data[,which(Group$Group=="0")]),rowMeans(GSE_Data[,which(Group$Group=="1")]))) -> a
  colnames(a) <- c("Control Mean","Case Mean")
  a[b$name,] -> a
  # write.table(cbind(a,b), file = outfile_name, quote = F,row.names = F, sep = "\t")
  limma_res[[i]] <- cbind(a,b,i,data_id)
}
/GE_Data[,which(Group$Group=="0")])
# rowMeans(GSE_Data[,which(Group$Group=="1")])


do.call(rbind,limma_res) -> #Limma results table


