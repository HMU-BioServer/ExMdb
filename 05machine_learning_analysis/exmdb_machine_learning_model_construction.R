library(pROC)
library(ggplot2)
library(Hmisc)
library(caret)
library(randomForestSRC)
library(e1071)
library(nnet)
library(adabag)
library(Boruta)
library(glmnet)
library(CoxBoost)
library(plsRcox)
library(superpc)
library(gbm)
library(survivalsvm)
library(timeROC)
library(survival)
library(survminer)
library(UpSetR)

#Data input requirements: the last two columns are OS , OS.time and both are numeric; there can not be too many features.

# model construction ------------------------------------------------------

## RSF ##

RSF.model=function(train,test,screeing.gene=colnames(train[,-col_OS])){
  if(length(screeing.gene) == ncol(train)-2 ){#As a solo model 
    set.seed(seed)
    RSF <- rfsrc(Surv(OS.time,OS)~.,data = train,
                 ntree = 400,  
                 splitrule = 'logrank',
                 importance = T,
                 proximity = T,
                 forest = T,
                 seed = seed)
  }
  if(length(screeing.gene) != ncol(train)-2  ){#As a composite model
    set.seed(seed)
    RSF <- rfsrc(Surv(OS.time,OS)~.,data = train[,c(screeing.gene,"OS.time","OS")],
                 ntree = 400,  
                 splitrule = 'logrank',
                 importance = T,
                 proximity = T,
                 forest = T,
                 seed = seed)
  }
  rs <- cbind(test[, c("OS.time","OS")],RS=predict(RSF,newdata = test[,screeing.gene])$predicted)
  cc <- as.numeric(summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1])
  
  #Filtering features based on importance
  screening.RSF = names( RSF$importance[RSF$importance>0] )
  
  RSF_reslut=list(Cindex=cc,
                  screening.RSF=screening.RSF)
  return(RSF_reslut)
}


## Lasso ##

Lasso.model=function(train,test,screeing.gene=colnames(train[,-col_OS])){
  if(length(screeing.gene) == ncol(train)-2 ){
    x1 <- as.matrix(train[,-col_OS])#Input data other than OS.time and OS
    x2 <- as.matrix(Surv(train$OS.time,train$OS))
  }
  if(length(screeing.gene) != ncol(train)-2  ){
    x1 <- as.matrix(train[,screeing.gene])
    x2 <- as.matrix(Surv(train$OS.time,train$OS))
  }
  set.seed(seed)
  Lasso = cv.glmnet(x1, x2,family = "cox",alpha=1,nfolds = 5)
  rs <- cbind(test[,c("OS.time","OS")],RS=as.numeric(predict(Lasso,type='link',newx=as.matrix(test[,screeing.gene]),s=Lasso$lambda.min)))    
  cc <- as.numeric(summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1])
  
  
  screening.Lasso = data.frame(features= rownames(coef(Lasso,s="lambda.min")) , coef=as.vector(coef(Lasso,s="lambda.min")) ) 
  screening.Lasso = screening.Lasso[screening.Lasso[,2] != 0,1] 
  Lasso_reslut=list(Cindex=cc,
                    screening.Lasso=screening.Lasso)
  return(Lasso_reslut)
}



## stepwiseCox.both ##
stepwiseCox.both.model=function(train,test,screeing.gene=colnames(train[,-col_OS])){ 
  if(length(screeing.gene) != ncol(train)-2  ){
    train=train[,c(screeing.gene,"OS.time","OS")]
    test=test[,c(screeing.gene,"OS.time","OS")]
  }
  set.seed(seed)
  stepwiseCox.both <- step(coxph(Surv(OS.time,OS)~.,train),direction = "both")
  rs <- cbind(test[,c("OS.time","OS")],RS=predict(stepwiseCox.both,type = 'risk',newdata = test))
  cc <- as.numeric(summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1])
  
  
  screening.stepwiseCox.both = summary(stepwiseCox.both)$coefficients
  screening.stepwiseCox.both = rownames(screening.stepwiseCox.both) 
  
  stepwiseCox.both.reslut=list(Cindex=cc,
                               screening.stepwiseCox.both=screening.stepwiseCox.both)
  return(stepwiseCox.both.reslut)
}



## stepwiseCox.backward ##
stepwiseCox.backward.model=function(train,test,screeing.gene=colnames(train[,-col_OS])){ 
  if(length(screeing.gene) != ncol(train)-2  ){#有输入筛选基因
    train=train[,c(screeing.gene,"OS.time","OS")]
    test=test[,c(screeing.gene,"OS.time","OS")]
  }
  set.seed(seed)
  stepwiseCox.backward <- step(coxph(Surv(OS.time,OS)~.,train),direction = "backward")
  rs <- cbind(test[,c("OS.time","OS")],RS=predict(stepwiseCox.backward,type = 'risk',newdata = test))
  cc <- as.numeric(summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1])
  
  
  screening.stepwiseCox.backward = summary(stepwiseCox.backward)$coefficients
  screening.stepwiseCox.backward = rownames(screening.stepwiseCox.backward) 
  
  stepwiseCox.backward.result=list(Cindex=cc,
                                   screening.stepwiseCox.backward=screening.stepwiseCox.backward)
}



## CoxBoost ##
CoxBoost.model=function(train,test,screeing.gene=colnames(train[,-col_OS])){
  if( length(screeing.gene) != ncol(train)-2  ){#有输入筛选基因
    train=train[,c(screeing.gene,"OS.time","OS")]
    test=test[,c(screeing.gene,"OS.time","OS")]
  }
  set.seed(seed)
  pen <- optimCoxBoostPenalty(train[,'OS.time'],train[,'OS'],as.matrix(train[,screeing.gene]),
                              trace=TRUE,start.penalty=500,parallel = T)
  cv.res <- cv.CoxBoost(train[,'OS.time'],train[,'OS'],as.matrix(train[,screeing.gene]),
                        maxstepno=500,K=10,type="verweij",penalty=pen$penalty)
  CoxBoost <- CoxBoost(train[,'OS.time'],train[,'OS'],as.matrix(train[,screeing.gene]),
                       stepno=cv.res$optimal.step,penalty=pen$penalty)
  rs <- cbind(test[,c("OS.time","OS")],RS=as.numeric(predict(CoxBoost,newdata=test[,screeing.gene], newOS.time=test[,"OS.time"], newOS=test[,"OS"], type="lp")))
  cc <- as.numeric(summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1])
  
  
  if( length(screeing.gene) < ncol(train)-2   ){
    CoxBoost.result=list(Cindex=cc)
    return(CoxBoost.result)
  }
  
  if( length(screeing.gene) == ncol(train)-2  ){
    coxboost.predicted = predict(CoxBoost,newdata=test, newOS.time=test[,"OS.time"], newOS=test[,"OS"], type="lp")
    screening.CoxBoost = CoxBoost$xnames[CoxBoost$coefficients[dim(CoxBoost$coefficients)[1],] !=0 ] 
    CoxBoost.result=list(Cindex=cc,
                         screening.CoxBoost=screening.CoxBoost)
    return(CoxBoost.result)
  }
  
}



## Enet(0.1 - 0.9) ##
Enet.model=function(train,test,screeing.gene=colnames(train[,-col_OS])){ 
  Enet.Cindex=list()
  if(length(screeing.gene) == ncol(train)-2 ){
    x1 <- as.matrix(train[,-col_OS])
    x2 <- as.matrix(Surv(train$OS.time,train$OS))
  }
  if(length(screeing.gene) != ncol(train)-2  ){
    x1 <- as.matrix(train[,screeing.gene])#
    x2 <- as.matrix(Surv(train$OS.time,train$OS))
  }
  
  for (alpha in seq(0.1,0.9,0.1)) {
    set.seed(seed)
    fit = cv.glmnet(x1, x2,family = "cox",alpha=alpha,nfolds = 10)
    rs <- cbind(test[,c(dim(test)[2]-1,dim(test)[2] )],
                RS=as.numeric(predict(fit,type='link',
                                      newx=as.matrix(test[,screeing.gene]),s=fit$lambda.min)))
    cc <- as.numeric(summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1])
    
    model.name=paste0('α',alpha) 
    Enet.Cindex[[model.name]] <- cc
  }
  return(Enet.Cindex)
}



## ridge ##
ridge.model=function(train,test,screeing.gene=colnames(train[,-col_OS])){
  if(length(screeing.gene) == ncol(train)-2 ){
    x1 <- as.matrix(train[,-col_OS])
    x2 <- as.matrix(Surv(train$OS.time,train$OS))
  }
  if(length(screeing.gene) != ncol(train)-2  ){
    x1 <- as.matrix(train[,screeing.gene])
    x2 <- as.matrix(Surv(train$OS.time,train$OS))
  }
  
  set.seed(seed)
  Ridge = cv.glmnet(x1, x2,family = "cox",alpha=0,nfolds = 10)
  rs <- cbind(test[,c(dim(test)[2]-1,dim(test)[2])],
              RS=as.numeric(predict(Ridge,type='link',
                                    newx=as.matrix(test[,screeing.gene]),s=Ridge$lambda.min)))
  cc <- as.numeric(summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1])
  
  
  return(cc)
}



## plsRcox ##
plsRcox.model=function(train,test,screeing.gene=colnames(train[,-col_OS])){ 
  
  if( length(screeing.gene) == ncol(train)-2 ){
    set.seed(seed)
    cv.plsRcox.res=cv.plsRcox(list(x=train[,-col_OS],time=train$OS.time,status=train$OS,nt=10))
    fit <- plsRcox(train[,-col_OS],time=train$OS.time,event=train$OS,nt=as.numeric(cv.plsRcox.res[5]))
    rs <- cbind(test[,col_OS],RS=as.numeric(predict(fit,type="lp",newdata=test[,-col_OS])))
    cc <- as.numeric(summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1])
    
  }
  if( length(screeing.gene) != ncol(train)-2  ){
    set.seed(seed)
    cv.plsRcox.res=cv.plsRcox(list(x=train[,screeing.gene],time=train$OS.time,status=train$OS),nt=10)
    fit <- plsRcox(train[,screeing.gene],time=train$OS.time,event=train$OS,nt=as.numeric(cv.plsRcox.res[5]))
    rs <- cbind(test[,c("OS.time","OS")],RS=as.numeric(predict(fit,type="lp",newdata=test[,screeing.gene])))
    cc <- as.numeric(summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1])
    
  }    
  return(cc)
}



## superpc ##
superpc.model=function(train,test,screeing.gene=colnames(train[,-col_OS])){ 
  if(length(screeing.gene) == ncol(train)-2 ){
    data <- list(x=t(train[,-col_OS]),y=train$OS.time,censoring.status=train$OS,featurenames=colnames(train)[-col_OS])
  }
  if(length(screeing.gene) != ncol(train)-2  ){
    data <- list(x=t(train[,screeing.gene]),y=train$OS.time,censoring.status=train$OS,featurenames=colnames(train)[screeing.gene])
  }
  set.seed(seed)
  fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
  cv.fit <- superpc.cv(fit,data,n.threshold = 20,#default 
                       n.fold = 10,
                       n.components=3,
                       min.features=length(screeing.gene),
                       max.features=nrow(data$x),
                       compute.fullcv= TRUE,
                       compute.preval=TRUE)
  
  test.superpc <- list(x=t(test[,screeing.gene]),y=test$OS.time,censoring.status=test$OS,featurenames=colnames(test)[screeing.gene])
  ff <- superpc.predict(fit,data,test.superpc,threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rs <- cbind(test[,c("OS.time","OS")],RS=rr)
  cc <- as.numeric(summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1])
  return(cc)
}



## survivalsvm ##
survivalsvm.model=function(train,test,screeing.gene=colnames(train[,-col_OS])){ 
  if(length(screeing.gene) == ncol(train)-2 ){
    surv.SVM = survivalsvm(Surv(OS.time,OS)~., data= train, gamma.mu = 1)
  }
  if(length(screeing.gene) != ncol(train)-2  ){
    surv.SVM = survivalsvm(Surv(OS.time,OS)~., data= train[,c(screeing.gene,"OS.time","OS")], gamma.mu = 1)
  }
  rs <- cbind(test[,c("OS.time","OS")],RS=as.numeric(predict(surv.SVM,test[,c(screeing.gene,"OS.time","OS")])$predicted))
  cc <- as.numeric(summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1])
  return(cc)
}



## stepwiseCox.forward ##
stepwiseCox.forward.model=function(train,test,screeing.gene=colnames(train[,-col_OS])){
  if(length(screeing.gene) == ncol(train)-2 ){
    set.seed(seed)
    stepwiseCox.forward <- step(coxph(Surv(OS.time,OS)~.,train),direction = "forward")
  }
  if(length(screeing.gene) != ncol(train)-2  ){
    set.seed(seed)
    stepwiseCox.forward <- step(coxph(Surv(OS.time,OS)~.,train[,c(screeing.gene,"OS.time","OS")]),direction = "forward")
  }
  rs <- cbind(test[,c("OS.time","OS")],RS=predict(stepwiseCox.forward,type = 'risk',newdata = test[,c(screeing.gene,"OS.time","OS")]))
  cc <- as.numeric(summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1])
  
  return(cc)
}



# best model select -------------------------------------------------------


rsf=function(train,test,screeing.gene=colnames(train[,-col_OS])){
  if(length(screeing.gene) == ncol(train)-2 ){
    set.seed(seed)
    RSF <- rfsrc(Surv(OS.time,OS)~.,data = train,
                 ntree = 400,  
                 splitrule = 'logrank',
                 importance = T,
                 proximity = T,
                 forest = T,
                 seed = seed)
  }
  if(length(screeing.gene) != ncol(train)-2  ){
    set.seed(seed)
    RSF <- rfsrc(Surv(OS.time,OS)~.,data = train[,c(screeing.gene,"OS.time","OS")],
                 ntree = 400,  
                 splitrule = 'logrank',
                 importance = T,
                 proximity = T,
                 forest = T,
                 seed = seed)
  }
  rs <- cbind(test[, c("OS.time","OS")],RS=predict(RSF,newdata = test[,screeing.gene])$predicted)
  cc <- as.numeric(summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1])
  
  
  screening.RSF = names( RSF$importance[RSF$importance>0] )
  
  RSF_reslut=list(Cindex=cc,
                  screening.RSF=screening.RSF)
  return(rs)
}




lasso=function(train,test,screeing.gene=colnames(train[,-col_OS])){
  if(length(screeing.gene) == ncol(train)-2 ){
    x1 <- as.matrix(train[,-col_OS])
    x2 <- as.matrix(Surv(train$OS.time,train$OS))
  }
  if(length(screeing.gene) != ncol(train)-2  ){
    x1 <- as.matrix(train[,screeing.gene])
    x2 <- as.matrix(Surv(train$OS.time,train$OS))
  }
  set.seed(seed)
  Lasso = cv.glmnet(x1, x2,family = "cox",alpha=1,nfolds = 5)
  rs <- cbind(test[,c("OS.time","OS")],RS=as.numeric(predict(Lasso,type='link',newx=as.matrix(test[,screeing.gene]),s=Lasso$lambda.min)))#newx输入数据也是除test的OS.time和OS以外的数据      
  cc <- as.numeric(summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1])
  
  
  screening.Lasso = data.frame(features= rownames(coef(Lasso,s="lambda.min")) , coef=as.vector(coef(Lasso,s="lambda.min")) ) 
  screening.Lasso = screening.Lasso[screening.Lasso[,2] != 0,1] 
  Lasso_reslut=list(Cindex=cc,
                    screening.Lasso=screening.Lasso)
  return(rs)
}




stepCox.both=function(train,test,screeing.gene=colnames(train[,-col_OS])){ 
  if(length(screeing.gene) != ncol(train)-2  ){
    train=train[,c(screeing.gene,"OS.time","OS")]
    test=test[,c(screeing.gene,"OS.time","OS")]
  }
  set.seed(seed)
  stepwiseCox.both <- step(coxph(Surv(OS.time,OS)~.,train),direction = "both")
  rs <- cbind(test[,c("OS.time","OS")],RS=predict(stepwiseCox.both,type = 'risk',newdata = test))
  cc <- as.numeric(summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1])
  
  
  screening.stepwiseCox.both = summary(stepwiseCox.both)$coefficients
  screening.stepwiseCox.both = rownames(screening.stepwiseCox.both) 
  
  stepwiseCox.both.reslut=list(Cindex=cc,
                               screening.stepwiseCox.both=screening.stepwiseCox.both)
  return(rs)
}




stepCox.back=function(train,test,screeing.gene=colnames(train[,-col_OS])){ 
  if(length(screeing.gene) != ncol(train)-2  ){
    train=train[,c(screeing.gene,"OS.time","OS")]
    test=test[,c(screeing.gene,"OS.time","OS")]
  }
  set.seed(seed)
  stepwiseCox.backward <- step(coxph(Surv(OS.time,OS)~.,train),direction = "backward")
  rs <- cbind(test[,c("OS.time","OS")],RS=predict(stepwiseCox.backward,type = 'risk',newdata = test))
  cc <- as.numeric(summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1])
  
  
  screening.stepwiseCox.backward = summary(stepwiseCox.backward)$coefficients
  screening.stepwiseCox.backward = rownames(screening.stepwiseCox.backward) 
  
  stepwiseCox.backward.result=list(Cindex=cc,
                                   screening.stepwiseCox.backward=screening.stepwiseCox.backward)
  return(rs)
}




coxboost=function(train,test,screeing.gene=colnames(train[,-col_OS])){
  if( length(screeing.gene) != ncol(train)-2  ){
    train=train[,c(screeing.gene,"OS.time","OS")]
    test=test[,c(screeing.gene,"OS.time","OS")]
  }
  set.seed(seed)
  pen <- optimCoxBoostPenalty(train[,'OS.time'],train[,'OS'],as.matrix(train[,screeing.gene]),
                              trace=TRUE,start.penalty=500,parallel = T)
  cv.res <- cv.CoxBoost(train[,'OS.time'],train[,'OS'],as.matrix(train[,screeing.gene]),
                        maxstepno=500,K=10,type="verweij",penalty=pen$penalty)
  CoxBoost <- CoxBoost(train[,'OS.time'],train[,'OS'],as.matrix(train[,screeing.gene]),
                       stepno=cv.res$optimal.step,penalty=pen$penalty)
  rs <- cbind(test[,c("OS.time","OS")],RS=as.numeric(predict(CoxBoost,newdata=test[,screeing.gene], newOS.time=test[,"OS.time"], newOS=test[,"OS"], type="lp")))
  cc <- as.numeric(summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1])
  
  
  if( length(screeing.gene) < ncol(train)-2   ){
    CoxBoost.result=list(Cindex=cc)
    return(rs)
  }
  
  if(length(screeing.gene) == ncol(train)-2  ){
    coxboost.predicted = predict(CoxBoost,newdata=test, newOS.time=test[,"OS.time"], newOS=test[,"OS"], type="lp")
    screening.CoxBoost = CoxBoost$xnames[CoxBoost$coefficients[dim(CoxBoost$coefficients)[1],] !=0 ] 
    CoxBoost.result=list(Cindex=cc,
                         screening.CoxBoost=screening.CoxBoost)
    return(rs)
  }
  
}




enet.α0.1=function(train,test,screeing.gene=colnames(train[,-col_OS])){ 
  x1 <- as.matrix(train[,-col_OS])
  x2 <- as.matrix(Surv(train$OS.time,train$OS))
  set.seed(seed)
  fit = cv.glmnet(x1, x2,family = "cox",alpha=0.1,nfolds = 10)
  rs <- cbind(test[,c(dim(test)[2]-1,dim(test)[2] )],
              RS=as.numeric(predict(fit,type='link',
                                    newx=as.matrix(test[,screeing.gene]),s=fit$lambda.min)))
  return(rs)
}
enet.α0.2=function(train,test,screeing.gene=colnames(train[,-col_OS])){ 
  x1 <- as.matrix(train[,-col_OS])
  x2 <- as.matrix(Surv(train$OS.time,train$OS))
  set.seed(seed)
  fit = cv.glmnet(x1, x2,family = "cox",alpha=0.2,nfolds = 10)
  rs <- cbind(test[,c(dim(test)[2]-1,dim(test)[2] )],
              RS=as.numeric(predict(fit,type='link',
                                    newx=as.matrix(test[,screeing.gene]),s=fit$lambda.min)))
  return(rs)
}
enet.α0.3=function(train,test,screeing.gene=colnames(train[,-col_OS])){ 
  x1 <- as.matrix(train[,-col_OS])
  x2 <- as.matrix(Surv(train$OS.time,train$OS))
  set.seed(seed)
  fit = cv.glmnet(x1, x2,family = "cox",alpha=0.1,nfolds = 10)
  rs <- cbind(test[,c(dim(test)[2]-1,dim(test)[2] )],
              RS=as.numeric(predict(fit,type='link',
                                    newx=as.matrix(test[,screeing.gene]),s=fit$lambda.min)))
  return(rs)
}
enet.α0.4=function(train,test,screeing.gene=colnames(train[,-col_OS])){ 
  x1 <- as.matrix(train[,-col_OS])
  x2 <- as.matrix(Surv(train$OS.time,train$OS))
  set.seed(seed)
  fit = cv.glmnet(x1, x2,family = "cox",alpha=0.1,nfolds = 10)
  rs <- cbind(test[,c(dim(test)[2]-1,dim(test)[2] )],
              RS=as.numeric(predict(fit,type='link',
                                    newx=as.matrix(test[,screeing.gene]),s=fit$lambda.min)))
  return(rs)
}
enet.α0.5=function(train,test,screeing.gene=colnames(train[,-col_OS])){ 
  x1 <- as.matrix(train[,-col_OS])
  x2 <- as.matrix(Surv(train$OS.time,train$OS))
  set.seed(seed)
  fit = cv.glmnet(x1, x2,family = "cox",alpha=0.1,nfolds = 10)
  rs <- cbind(test[,c(dim(test)[2]-1,dim(test)[2] )],
              RS=as.numeric(predict(fit,type='link',
                                    newx=as.matrix(test[,screeing.gene]),s=fit$lambda.min)))
  return(rs)
}
enet.α0.6=function(train,test,screeing.gene=colnames(train[,-col_OS])){ 
  x1 <- as.matrix(train[,-col_OS])
  x2 <- as.matrix(Surv(train$OS.time,train$OS))
  set.seed(seed)
  fit = cv.glmnet(x1, x2,family = "cox",alpha=0.1,nfolds = 10)
  rs <- cbind(test[,c(dim(test)[2]-1,dim(test)[2] )],
              RS=as.numeric(predict(fit,type='link',
                                    newx=as.matrix(test[,screeing.gene]),s=fit$lambda.min)))
  return(rs)
}
enet.α0.7=function(train,test,screeing.gene=colnames(train[,-col_OS])){ 
  x1 <- as.matrix(train[,-col_OS])
  x2 <- as.matrix(Surv(train$OS.time,train$OS))
  set.seed(seed)
  fit = cv.glmnet(x1, x2,family = "cox",alpha=0.1,nfolds = 10)
  rs <- cbind(test[,c(dim(test)[2]-1,dim(test)[2] )],
              RS=as.numeric(predict(fit,type='link',
                                    newx=as.matrix(test[,screeing.gene]),s=fit$lambda.min)))
  return(rs)
}
enet.α0.8=function(train,test,screeing.gene=colnames(train[,-col_OS])){ 
  x1 <- as.matrix(train[,-col_OS])
  x2 <- as.matrix(Surv(train$OS.time,train$OS))
  set.seed(seed)
  fit = cv.glmnet(x1, x2,family = "cox",alpha=0.1,nfolds = 10)
  rs <- cbind(test[,c(dim(test)[2]-1,dim(test)[2] )],
              RS=as.numeric(predict(fit,type='link',
                                    newx=as.matrix(test[,screeing.gene]),s=fit$lambda.min)))
  return(rs)
}
enet.α0.9=function(train,test,screeing.gene=colnames(train[,-col_OS])){ 
  x1 <- as.matrix(train[,-col_OS])
  x2 <- as.matrix(Surv(train$OS.time,train$OS))
  set.seed(seed)
  fit = cv.glmnet(x1, x2,family = "cox",alpha=0.1,nfolds = 10)
  rs <- cbind(test[,c(dim(test)[2]-1,dim(test)[2] )],
              RS=as.numeric(predict(fit,type='link',
                                    newx=as.matrix(test[,screeing.gene]),s=fit$lambda.min)))
  return(rs)
}




ridge=function(train,test,screeing.gene=colnames(train[,-col_OS])){
  if(length(screeing.gene) == ncol(train)-2 ){
    x1 <- as.matrix(train[,-col_OS])
    x2 <- as.matrix(Surv(train$OS.time,train$OS))
  }
  if(length(screeing.gene) != ncol(train)-2  ){
    x1 <- as.matrix(train[,screeing.gene])
    x2 <- as.matrix(Surv(train$OS.time,train$OS))
  }
  
  set.seed(seed)
  Ridge = cv.glmnet(x1, x2,family = "cox",alpha=0,nfolds = 10)
  rs <- cbind(test[,c(dim(test)[2]-1,dim(test)[2])],
              RS=as.numeric(predict(Ridge,type='link',
                                    newx=as.matrix(test[,screeing.gene]),s=Ridge$lambda.min)))
  cc <- as.numeric(summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1])
  
  
  return(rs)
}




plsRCox=function(train,test,screeing.gene=colnames(train[,-col_OS])){ 
  
  if( length(screeing.gene) == ncol(train)-2 ){
    set.seed(seed)
    cv.plsRcox.res=cv.plsRcox(list(x=train[,-col_OS],time=train$OS.time,status=train$OS,nt=10))
    fit <- plsRcox(train[,-col_OS],time=train$OS.time,event=train$OS,nt=as.numeric(cv.plsRcox.res[5]))
    rs <- cbind(test[,col_OS],RS=as.numeric(predict(fit,type="lp",newdata=test[,-col_OS])))
    cc <- as.numeric(summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1])
    
  }
  if( length(screeing.gene) != ncol(train)-2  ){
    set.seed(seed)
    cv.plsRcox.res=cv.plsRcox(list(x=train[,screeing.gene],time=train$OS.time,status=train$OS),nt=10)
    fit <- plsRcox(train[,screeing.gene],time=train$OS.time,event=train$OS,nt=as.numeric(cv.plsRcox.res[5]))
    rs <- cbind(test[,c("OS.time","OS")],RS=as.numeric(predict(fit,type="lp",newdata=test[,screeing.gene])))
    cc <- as.numeric(summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1])
    
  }    
  return(rs)
}




superpc=function(train,test,screeing.gene=colnames(train[,-col_OS])){ 
  if(length(screeing.gene) == ncol(train)-2 ){
    data <- list(x=t(train[,-col_OS]),y=train$OS.time,censoring.status=train$OS,featurenames=colnames(train)[-col_OS])
  }
  if(length(screeing.gene) != ncol(train)-2  ){
    data <- list(x=t(train[,screeing.gene]),y=train$OS.time,censoring.status=train$OS,featurenames=colnames(train)[screeing.gene])
  }
  set.seed(seed)
  fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) 
  cv.fit <- superpc.cv(fit,data,n.threshold = 20, 
                       n.fold = 10,
                       n.components=3,
                       min.features=length(screeing.gene),
                       max.features=nrow(data$x),
                       compute.fullcv= TRUE,
                       compute.preval=TRUE)
  
  test.superpc <- list(x=t(test[,screeing.gene]),y=test$OS.time,censoring.status=test$OS,featurenames=colnames(test)[screeing.gene])
  ff <- superpc.predict(fit,data,test.superpc,threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rs <- cbind(test[,c("OS.time","OS")],RS=rr)
  cc <- as.numeric(summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1])
  return(rs)
}



GBM=function(train,test,screeing.gene=colnames(train[,-col_OS])){
  if(length(screeing.gene) == ncol(train)-2 ){
    set.seed(seed)
    GBM <- gbm(formula = Surv(OS.time,OS)~.,data = train,distribution = 'coxph',
               n.trees = 1000,
               shrinkage = 0.01,
               cv.folds = 10,
               verbose = FALSE)
    
    best <- which.min(GBM$cv.error)
    set.seed(seed)
    GBM <- gbm(formula = Surv(OS.time,OS)~.,data = train,distribution = 'coxph',
               n.trees = best,
               shrinkage = 0.001,
               cv.folds = 10,
               verbose = FALSE)
  }
  
  if(length(screeing.gene) != ncol(train)-2  ){
    set.seed(seed)
    GBM <- gbm(formula = Surv(OS.time,OS)~.,data = train[, c(screeing.gene,"OS.time","OS") ],distribution = 'coxph',
               n.trees = 1000,
               shrinkage = 0.01,
               cv.folds = 10,
               verbose = FALSE)
    
    best <- which.min(GBM$cv.error)
    if(best==1){
      second_min = unique(GBM$cv.error[order(GBM$cv.error,decreasing = F)][2]) 
      best <- which(GBM$cv.error==second_min)
    }
    set.seed(seed)
    GBM <- gbm(formula = Surv(OS.time,OS)~.,data = train[,c(screeing.gene,"OS.time","OS")],distribution = 'coxph',
               n.trees = best,
               shrinkage = 0.001,
               cv.folds = 10,
               verbose = FALSE)
  }
  rs <- cbind(test[,c("OS.time","OS")],RS=as.numeric(predict(GBM,test[,c(screeing.gene,"OS.time","OS")],n.trees = best,type = 'link')))
  cc <- as.numeric(summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1])
  
  return(rs)
}




survivalSVM=function(train,test,screeing.gene=colnames(train[,-col_OS])){ 
  if(length(screeing.gene) == ncol(train)-2 ){
    surv.SVM = survivalsvm(Surv(OS.time,OS)~., data= train, gamma.mu = 1)
  }
  if(length(screeing.gene) != ncol(train)-2  ){
    surv.SVM = survivalsvm(Surv(OS.time,OS)~., data= train[,c(screeing.gene,"OS.time","OS")], gamma.mu = 1)
  }
  rs <- cbind(test[,c("OS.time","OS")],RS=as.numeric(predict(surv.SVM,test[,c(screeing.gene,"OS.time","OS")])$predicted))
  cc <- as.numeric(summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1])
  return(rs)
}




stepCox.forward=function(train,test,screeing.gene=colnames(train[,-col_OS])){
  if(length(screeing.gene) == ncol(train)-2 ){
    set.seed(seed)
    stepwiseCox.forward <- step(coxph(Surv(OS.time,OS)~.,train),direction = "forward")
  }
  if(length(screeing.gene) != ncol(train)-2  ){
    set.seed(seed)
    stepwiseCox.forward <- step(coxph(Surv(OS.time,OS)~.,train[,c(screeing.gene,"OS.time","OS")]),direction = "forward")
  }
  rs <- cbind(test[,c("OS.time","OS")],RS=predict(stepwiseCox.forward,type = 'risk',newdata = test[,c(screeing.gene,"OS.time","OS")]))
  cc <- as.numeric(summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1])
  
  return(rs)
}








