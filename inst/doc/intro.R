## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
pages=data.frame(v1=c(1,1,1,2,2,3,4),v2=c(2,3,4,3,4,4,2))

pagerk <-function(pages){
  #计算邻接矩阵
  adjacencyMatrix<-function(pages){
    n<-max(apply(pages,2,max))
    A <- matrix(0,n,n)
    for(i in 1:nrow(pages)) A[pages[i,]$dist,pages[i,]$src]<-1
    A
  }
  #变换概率矩阵
  probabilityMatrix<-function(G){
    cs <- colSums(G)
    cs[cs==0] <- 1
    n <- nrow(G)
    A <- matrix(0,nrow(G),ncol(G))
    for (i in 1:n) A[i,] <- A[i,] + G[i,]/cs
    A
  }

  #考虑d时的变换概率矩阵
  dProbabilityMatrix<-function(G,d=0.85){
    cs <- colSums(G)
    cs[cs==0] <- 1
    n <- nrow(G)
    delta <- (1-d)/n
    A <- matrix(delta,nrow(G),ncol(G))
    for (i in 1:n) A[i,] <- A[i,] + d*G[i,]/cs
    A
  }
  
  #递归计算矩阵特征值
  eigenMatrix<-function(G,iter=100){
    iter<-10
    n<-nrow(G)
    x <- rep(1,n)
    for (i in 1:iter) x <- G %*% x
    x/sum(x)
  }
  #无阻尼系数下的pagerank值
  names(pages)<-c("src","dist")
  A<-adjacencyMatrix(pages)
  G<-probabilityMatrix(A)
  q1<-eigenMatrix(G,100)
  
  #阻尼系数d=0.85下的pagerank值
  A<-adjacencyMatrix(pages)
  G<-dProbabilityMatrix(A)
  q2<-eigenMatrix(G,100)
  q<- cbind(q1,q2)
  return(q)
}
pagerk(pages)

## -----------------------------------------------------------------------------
library(ROCR)
get_confusion_stat <- function(pred,y_real,threshold=0.5){
  # auc
  tmp <- prediction(as.vector(pred),y_real)
  auc <- unlist(slot(performance(tmp,'auc'),'y.values'))
  # statistic
  pred_new <- as.integer(pred>threshold) 
  tab <- table(pred_new,y_real)
  if(nrow(tab)==1){
    print('preds all zero !')
    return(0)
  }
  TP <- tab[2,2]
  TN <- tab[1,1]
  FP <- tab[2,1]
  FN <- tab[1,2]
  accuracy <- round((TP+TN)/(TP+FN+FP+TN),4)
  recall_sensitivity <- round(TP/(TP+FN),4)
  precision <- round(TP/(TP+FP),4)
  specificity <- round(TN/(TN+FP),4)
  # 添加，预测的负例占比（业务解释：去除多少的样本，达到多少的recall）
  neg_rate <- round((TN+FN)/(TP+TN+FP+FN),4)
  re <- list('AUC' = auc,
             'Confusion_Matrix'=tab,
             'Statistics'=data.frame(value=c('accuracy'=accuracy,
                                             'recall_sensitivity'=recall_sensitivity,
                                             'precision'=precision,
                                             'specificity'=specificity,
                                             'neg_rate'=neg_rate)))
  return(re)
}
pred=c(-1,1,1,1,1,1,1,-1,-1,1)
y_real=c(1,-1,1,1,1,1,1,-1,1,1)
get_confusion_stat(pred,y_real)

