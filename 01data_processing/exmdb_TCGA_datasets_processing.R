library(survival)

setwd("TCGA datasets")
files=list.files()
index<-read.table('gencode.annotation.file',header = T)
#Statistical cancer type
cancer.type=c()
for(i in files){
  temp=strsplit(i,split = "\\.")[[1]][1]
  cancer.type=c(cancer.type,temp)
}
cancer.type=unique(cancer.type)

error_datasets=c()

for(i in cancer.type){
  # read data --------------------------------------------------------------------
  
  data=read.table( paste(i,".mirna.tsv.gz",sep = "") ,header = T,row.names = 1,check.names = F)
  #rownames processing
  rowname=rownames(data)
  rowname=gsub("-", "_", rowname)
  rownames(data)=rowname
  
  survival.data=read.table( paste(i,".survival.tsv",sep = "") ,header = T)
  #mirna.data=read.table( paste(i,".mirna.tsv.gz",sep = "") ,header = T,row.names = 1,check.names = F)
  
  print( paste(i,"imported") )
  
  # Data cleaning --------------------------------------------------------------------
  #The para-cancer samples were deleted
  sample.names=colnames(data)
  if( sum(grep("-11", sample.names)) != 0  ){
    sample.names = sample.names[-grep("-11", sample.names)] 
  }
  
  data = data[,sample.names]
  
  #Filter out some low-quality genes (not expressed in 80% of samples)
  gene.proportion=apply(data,1,function(x){
    rate = sum(x == 0)/ncol(data)
  })
  data = data[!gene.proportion>0.8,] 
  
  
  data=as.data.frame(t(data))
  data$sample=rownames(data)
  data=merge(data,survival.data[,-3],by="sample")
  rownames(data)=data[,1]
  data=data[,-1]
  
  print( paste(i,"calculate Cox") )
  # cox ------------------------------------------------------------------
  varU=colnames(data[ ,1:(dim(data)[2]-2) ]) #
  Result<-c()
  for(j in 1:length(varU)){
    formula = as.formula( paste("Surv(OS.time,OS)~",varU[j],sep = "") ) 
    fit<-coxph(formula,data=data)
    fitSum<-summary(fit)
    result1<-c()
    result1<-rbind(result1,cbind(fitSum$coef,fitSum$conf.int)) 
    Result<-rbind(Result,result1)
  }
  
  Result=cbind(Result,id=rownames(Result))  
  write.table(Result, paste("Cox Results(miRNA)\\",i,"-Cox.txt",sep="") ,quote = F,sep = "\t", row.names=F)
  
  print( paste(i,"Done!") )
  
  if(nrow(data) <50 ){next}
# model input --------------------------------------------------------------------
  cox.result=Result[order(Result[,5] , decreasing = F),]
  
  if(nrow(data) >= 50 & nrow(data) <= 100){signatures=cox.result[1:10,]}
  if(nrow(data) > 100 & nrow(data) <= 150){signatures=cox.result[1:15,]}
  if(nrow(data) > 150& nrow(data) <= 200){signatures=cox.result[1:20,]}
  if(nrow(data) > 200& nrow(data) <= 250){signatures=cox.result[1:25,]}
  if(nrow(data) > 250& nrow(data) <= 300){signatures=cox.result[1:30,]}
  if(nrow(data) > 300){signatures=cox.result[1:30,]}
  
  data=data[,c(signatures[,"id"],"OS","OS.time")]

  alive=data[data$OS == 0,]
  dead=data[data$OS == 1,]
  seed=1234
  set.seed(seed)
  
  test = rbind(alive[sample(nrow(alive), ceiling(nrow(alive))/3 ),],
               dead[sample(nrow(dead), ceiling(nrow(dead))/3 ),]) 
  train = data[!(rownames(data) %in% rownames(test)),]
  col_OS=c(ncol(train)-1,ncol(train))
  
  # model result --------------------------------------------------------------------
  
  solo_Cindex=list()


  result <- tryCatch({
    solo_Cindex=list(
     
      rsf=RSF.model(train,test)[["Cindex"]],
      lasso=Lasso.model(train,test)[["Cindex"]],
      stepCox.both=tepwiseCox.both.model(train,test)[["Cindex"]],
      stepCox.back=stepwiseCox.backward.model(train,test)[["Cindex"]],
      coxboost=CoxBoost.model(train,test)[["Cindex"]],
      enet=Enet.model(train,test),
      ridge=ridge.model(train,test),
      plsRCox=plsRcox.model(train,test),
      superpc=superpc.model(train,test),
      survivalSVM=survivalsvm.model(train,test) ,
      stepCox.forward= stepwiseCox.forward.model(train,test)
    )
  }, error = function(e) {
   
    print("Error")
  })
  if ( sum(result == "Error")==1 ) {
    error_datasets=c(error_datasets,i)
    next
  }
  
  
  data_Cindex=t( as.data.frame(solo_Cindex))
  colnames(data_Cindex)=i
  
  data_Cindex[is.na(data_Cindex)]=0
  data_Cindex[data_Cindex<0.5]=1-data_Cindex[data_Cindex<0.5]
  data_Cindex[data_Cindex==0.5 | data_Cindex==1]=0
  data_Cindex[data_Cindex==0]=NA
  
  

  best.model=names(data_Cindex[order(data_Cindex[,1],decreasing = T),])[1]
  args <- list(train = train, test = test)
  rs=do.call(best.model,args)
  
 
  threshold=surv_cutpoint(rs,
                          time = "OS.time",
                          event = "OS",
                          variables = "RS")[1]
  Type=ifelse(rs$RS> threshold$RS$estimate,"High Risk","Low Risk")
  
  rs$Type=as.factor(Type)
  survdiff=survdiff(Surv(OS.time,OS)~Type,rs)
  
  #tROC
  tROC <- timeROC(T = rs$OS.time, 
                  delta = rs$OS, 
                  marker = rs$RS, cause = 1, 
                  weighting = "marginal", 
                  times = c(365,1095,1825), iid = T)
  
  
  save.image( paste("model Results(miRNA)/",i,".RData",sep = "") )
}








