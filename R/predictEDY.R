#' Prediction of Extreme Downregulation of chromosome Y (EDY) 
#' from methylation data
#'
#' Predicts individuals' EDY status
#' @title predictEDY
#' @param x A matrix with CpGs in rows, samples in columns and 
#' CpG names in the rowname or a data.frame with individuals/samples
#' in columns, being the first column the CpG name or an ExpressionSet.ontrol")
#' @import Biobase
#' @export predictEDY
#' @return A vector with EDY status (No, Yes).
#' 
predictEDY <- function(x, ...){
  if (is.matrix(x)) {
    sel <- intersect(rownames(x), colnames(train))
    x.sel <- t(x[sel, ])
  }
  else if (inherits(x, "ExpressionSet")){ 
    sel <- intersect(featureNames(x), colnames(train))
    x.sel <- t(exprs(x)[sel, ])
  }
  else if (is.data.frame(x)){
    sel <- intersect(x[,1], colnames(train))
    x.sel <- as.matrix(t(x[x[,1]%in%sel, -1]))
  }
  else {
    stop("'x' must be a matrix, data.frame or an ExpressionSet")
  }
  
  if (length(sel) < 20)
    stop("There are few or no CpGs in chromosome Y")
  
  train.subset <- as.matrix(train[,sel])
  mod <- glmnet::glmnet(x = train.subset,
                        y = train[,1],
                        family="binomial",
                        alpha=0.5, lambda=0.02)
  edy.test <- mod %>% predict(as.matrix(test[,sel]), 
                            type="class") %>% as.factor()
  
  cat("\n")
  cat("EDY predictive capacity of your model (TCGA test dataset) \n")
  cat("--------------------------------------------------------- \n \n")
  tt <- confusionMatrix(edy.test, test$EDY)
  print(tt$table)
  accuracy <- round(tt$overall*100, 1)
  cat("\n")
  cat(paste0("   Accuracy : ", accuracy[1], "% \n"))
  cat(paste0("      CI95% : (", accuracy[3], ", ", accuracy[4], ")"))
  
  edy.pred <- mod %>% predict(x.sel, type="class") %>% as.factor()
  return(edy.pred)
}
  