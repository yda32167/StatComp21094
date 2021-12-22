#' This is some descriptio of this function.
#' @title Use R to calculate AUC value.
#' 
#' @description AUC (Area Under Curve) is defined as the area under the ROC curve and the coordinate axis.
#' 
#' @details Considering and not considering the situation of the damping system.
#' 
#' @param pred is a vector with two classes of data, and it is the predicted value of the result.
#'
#' @param y_real is a vector with two classes of data, and it is the true value of the result.
#'
#' @import knitr
#'
#' @import ROCR
#'
#' @return a list
#' @export
#' @examples
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