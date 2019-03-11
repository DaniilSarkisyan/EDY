##' Extreme Downregulation of chromosome Y (EDY)
##'
##' Detection of individuals with EDY
##' @title EDY
##' @param x 
##' @param 
##' 
##' @export getEDY
##' @return ...

getEDY <- function(x, gender.var, male.key, gene.key, log = TRUE, group.var, control.key, ...){
  
  #Filter males in the expression set
  if (!missing(gender.var)) {
    ii <- which(varLabels(x)==gender.var)
    x <- x[,pData(x)[,ii]==male.key]
  } else {warning("male.key not specified. Performing analysis with the whole dataset")}
  
  #Add column with the hgnc symbol information to fData
  if (gene.key != "hgnc_symbol"){
      fData(x) <- merge(annot, fData(x), by.x = "hgnc_symbol", by.y = gene.key)}
  
  annot.expr <- fData(x)
  
  #Select from genes in gene.expr those that we know the hgnc symbol
  gene.expr <- exprs(x)[rownames(exprs(x))%in%annot.expr[, "ID"],]
  #Replace gene ID for hgnc symbol
  rownames(gene.expr) <- annot.expr[, 'hgnc_symbol']
  #Select those genes that belong to chrY
  exprY <- gene.expr[annot.expr[, 'hgnc_symbol']%in%chrY$hgnc_symbol,]
  #Select those genes that belong to the rest of the genome
  exprRef <- gene.expr[annot.expr[, 'hgnc_symbol']%in%chrRef$hgnc_symbol,]
  
  #Apply EDY formulae: 
  if (log) {
    Ry <- sweep(exprY, 2, FUN="-", apply(exprRef, 2, mean))
  }else{ 
    Ry <- sweep(log2(exprY), 2, FUN="-", 
                apply(log2(exprRef), 2, mean))
  }
  
  EDYcontinuous <- apply(Ry, 2, mean)
  
  if (!missing(group.var)&&!missing(control.key)){
    controls <- EDYcontinuous[pData(x)[,group.var]==control.key]
  } else {
    controls <- EDYcontinuous
    warning("No control group specified")
    }
  
  thresh <- median(controls) - 1.2*IQR(controls)
  EDY <- cut(EDYcontinuous, c(-Inf, thresh, Inf), 
             labels = c("Yes", "No"))
  
  #output
  names(EDY) <- names(EDYcontinuous)
  ans <- list(Ry=t(Ry), EDY=EDY, threshold=thresh, eSet=x)
  ans
}




