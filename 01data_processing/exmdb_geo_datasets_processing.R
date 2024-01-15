library(survival)

setwd("GEO datasets")
index<-read.table('gencode.v22.annotation.gene.probeMap',header = T)
index$id <- sapply(strsplit(index$id, "\\."), function(x) x[1])


files=list.files()
model_error_datasets=c()
for(i in files){
  #load data
  load(paste("GEO datasets/",i,sep = ""))
  #Determine whether the data is empty
  if( any(dim(exp_matrix)==0) ){next}
  if( any(dim(this_sur)==0) ){next}
  #Determining whether miRNA expression or mRNA expression
  contains_hsa <- grepl("hsa-", rownames(exp_matrix))
  if(sum(contains_hsa) > 20 ){exp_type="miRNA"}
  if(sum(contains_hsa) == 0 ){exp_type="mRNA"}
  
  #Replace missing value NA with 0
  exp_matrix[is.na(exp_matrix)] <- 0
  #Replace the infinite value with 0
  exp_matrix[exp_matrix==-Inf]=0
  exp_matrix[exp_matrix==Inf]=0
  
  #Integrate presentation data and survival data
  data=t(exp_matrix)
  data=data.frame(data,samples=rownames(data))
  colnames(this_sur)=c("samples","OS.time","OS")
  data=merge(data,this_sur,by="samples")
  rownames(data)=data[,1]
  data=data[,-1]
  
  #Skip samples smaller than 50
  if( nrow(data) <50 ){next}
  # Univariate Cox ------------------------------------------------------------------
  # varU=colnames(data[ ,1:(dim(data)[2]-2) ]) 
  # Result<-c()
  # for(j in 1:length(varU)){
  #   formula = as.formula( paste("Surv(OS.time,OS)~",varU[j],sep = "") ) 
  #   fit<-coxph(formula,data=data)
  #   fitSum<-summary(fit)
  #   result1<-c()
  #   result1<-rbind(result1,cbind(fitSum$coef,fitSum$conf.int))  
  #   Result<-rbind(Result,result1)
  # }
  # 
  # Result=cbind(Result,id=rownames(Result))  
  # 
  # #ID conversion (mRNA is required)
  # if(exp_type == "mRNA"){Result=merge(index[,c(1,2)],Result,by="id")}
  # 
   dataset=strsplit(i,split = "_")[[1]][2]
  Result=read.table(paste("Cox Results(GEO)\\",dataset,"_",exp_type,"_Cox.txt",sep="") , header=T,sep="\t",check.names=F)
  
  # model input --------------------------------------------------------------------
  cox.result=Result[order(Result[,"Pr(>|z|)"]  , decreasing = F),]#Sort by p-value
  #Select the signature of the input model：
  if(nrow(data) >= 50 & nrow(data) <= 100){signatures=cox.result[1:10,]}
  if(nrow(data) > 100 & nrow(data) <= 150){signatures=cox.result[1:15,]}
  if(nrow(data) > 150& nrow(data) <= 200){signatures=cox.result[1:20,]}
  if(nrow(data) > 200& nrow(data) <= 250){signatures=cox.result[1:25,]}
  if(nrow(data) > 250& nrow(data) <= 300){signatures=cox.result[1:30,]}
  if(nrow(data) > 300){signatures=cox.result[1:30,]}
  
  data=data[,c(signatures[,"id"],"OS","OS.time")]
  data=data[!data$OS.time==0,]# Delete samples with a lifetime of 0
  #Select the training set and the test set
  alive=data[data$OS == 0,]
  dead=data[data$OS == 1,]
  seed=1234
  set.seed(seed)
  
  test = rbind(alive[sample(nrow(alive), ceiling(nrow(alive))/3 ),],#1/3 of the living and dead samples were selected as the test set
               dead[sample(nrow(dead), ceiling(nrow(dead))/3 ),]) 
  train = data[!(rownames(data) %in% rownames(test)),]#Select the remaining samples as the training set
  col_OS=c(ncol(train)-1,ncol(train))
  
  
  # run model --------------------------------------------------------------------
  
  solo_Cindex=list()
  
  
  model_result <- tryCatch({
    solo_Cindex=list(
     
      rsf=RSF.model(train,test)[["Cindex"]],
      lasso=Lasso.model(train,test)[["Cindex"]],
      stepCox.both=stepwiseCox.both.model(train,test)[["Cindex"]],
      stepCox.back=stepwiseCox.backward.model(train,test)[["Cindex"]],
      coxboost=CoxBoost.model(train,test)[["Cindex"]],
      enet=Enet.model(train,test),
      ridge=ridge.model(train,test),
      plsRCox=plsRcox.model(train,test),
      superpc=superpc.model(train,test),
      survivalSVM=survivalsvm.model(train,test),
      stepCox.forward=stepwiseCox.forward.model(train,test)
    )
  }, error = function(e) {
    
    print("Error")
  })
  
  if ( sum(model_result == "Error")==1 ) {
    model_error_datasets=c(model_error_datasets,i)
    next
  }
  
  
  data_Cindex=t( as.data.frame(solo_Cindex))
  colnames(data_Cindex)=i
  # Cindex adjustment
  data_Cindex[is.na(data_Cindex)]=0
  data_Cindex[data_Cindex<0.5]=1-data_Cindex[data_Cindex<0.5]
  data_Cindex[data_Cindex==0.5 | data_Cindex==1]=0
  data_Cindex[data_Cindex==0]=NA
  
  
  # survival result --------------------------------------------------------------------
  #Select the prediction results of the best model
  best.model=names(data_Cindex[order(data_Cindex[,1],decreasing = T),])[1]
  args <- list(train = train, test = test)
  rs=do.call(best.model,args)
  
  threshold=surv_cutpoint(rs,
                          time = "OS.time",
                          event = "OS",
                          variables = "RS")[1]
  Type=ifelse(rs$RS> threshold$RS$estimate,"High Risk","Low Risk")
  
  rs$Type=as.factor(Type)
  survdiff=survdiff(Surv(OS.time,OS)~Type,rs)#$pvalue
  
  #tROC
  tROC_result <- tryCatch({
    tROC <- timeROC(T = rs$OS.time, 
                    delta = rs$OS, 
                    marker = rs$RS, cause = 1, 
                    weighting = "marginal", 
                    times = c(365,1095,1825), iid = T)
  }, error = function(e) {
    print("Error")
  })
  if ( sum(tROC_result == "Error")==1 ) {
    tROC <- timeROC(T = rs$OS.time, 
                    delta = rs$OS, 
                    marker = rs$RS, cause = 1, 
                    weighting = "marginal", 
                    times = c(365,1095), iid = T)
  }
  

  
  save.image( paste("model Results(GEO)/",dataset,"_",exp_type,".RData",sep = "") )
 }


