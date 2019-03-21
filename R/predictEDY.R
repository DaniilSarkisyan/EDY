
predictEDY <- function(x, ...){
  if (is.matrix(x)) {
    sel <- intersect(rownames(x), colnames(train))
    x.sel <- t(x[sel, ])
  }
  else if (inherits(x, "ExpressionSet")){ 
    sel <- intersect(feaureNames(x), colnames(train))
    x.sel <- t(exprs(x)[sel, ])
  }
  else if (is.data.frame(x)){
    sel <- intersect(x[,1], colnames(train))
    x.sel <- as.matrix(t(x[x[,1]%in%sel, -1]))
  }
  else {
    stop("'x' must be a matrix, data.frame or an ExpressionSet")
  }
  
  train.subset <- train[,sel]
  mod <- glmnet::glmnet(x=as.matrix(train.subset), 
                      y=train[,1],
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
  