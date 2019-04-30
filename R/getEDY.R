#' Detection of individuals with Extreme Downregulation of chromosome Y (EDY)
#' from transcriptomic data
#'
#' To assess EDY, we calculate the mean expression of chrY genes divided by the
#' mean expression of genes in the rest of the genome (for each individual).
#' Now, the following formulae is applied only with individuals in the control
#' group, or with all the individuals in case it does not exist such group:
#' \deqn{Threshold = median - 1.2 · IQR}
#' \emph{*IQR: Interquartile range}
#' 
#' It indicates the trheshold from which down to consider EDY.
#'
#' @title getEDY
#' @param x An ExpressionSet from microarray experiment.
#' @param gender.var A string indicating the name of the column in the table of
#'   phenotype (accessed by \code{pData}) that contains the gender of the
#'   individual.
#' @param male.key A string indicating the symbol that idendtifies males (eg.:
#'   "male", "M", ...).
#' @param gene.key A string indicating the name of the column that contains the
#'   Gene Symbol in the table of features (accessed by \code{fData}).
#' @param coef Numerical. Value to consider an outlier when calling LOY. Default
#'   correspond to 1.2 which keeps 5\% of data in the normal case.
#' @param log Logical. It is set to \code{TRUE} if the gene expression values
#'   are given in a logarithmic scale in the ExpressionSet and \code{FALSE}
#'   otherwise. Default is \code{TRUE}.
#' @param group.var A string indicating the name of the column that contains the
#'   information about if the individual is case or control.
#' @param control.key A string indicating the symbol that identifies control
#'   group (eg.: "control").
#'   
#' @import Biobase
#' @export getEDY
#' @return A list containing 4 objects: 
#' \itemize{ 
#'   \item \code{EDY}: a factor with two levels, \code{NO} and \code{YES} , 
#'     that indicates whether the individual has EDY or not.
#'   \item \code{EDYcontinuous}: a vector with individual relative expression
#'     of chromosome Y with regard the autosomal genes.
#'   \item \code{threshold}: a number indicating the threshold from which down 
#'     an individual is considered to have EDY. 
#'   \item \code{eSet}: an ExpressionSet that contains only males from the
#'   initial ExpressionSet. It contains a new column named `hgnc_symbol` in
#'   `fData` in case it did not have it from before, which has the HUGO gene
#'   nomenclature comittee symbols.
#'   
#'   }

getEDY <- function(x, gender.var, male.key, gene.key, coef=1.2, 
                   log = TRUE, group.var, control.key, experiment.type, ...){
  
  object.type <- tolower(class(x)[1])
  experiment.type <- tolower(experiment.type)
  
  if (object.type == "expressionset"){
    
      #Filter males in the expression set
      if (!missing(gender.var)) {
        ii <- which(varLabels(x)==gender.var)
        x <- x[,pData(x)[,ii]==male.key]
        x <- x[,!is.na(pData(x)[,ii])]
      } else {warning("male.key not specified. Performing analysis with the whole dataset")}
      
      #Add column with the hgnc symbol information to fData
      if (gene.key != "hgnc_symbol"){
        fData(x)$id.feature <- featureNames(x)
        fData(x) <- merge(EDY::annot, fData(x), 
                          by.x = "hgnc_symbol", 
                          by.y = gene.key)
        }
      
      annot.expr <- fData(x)
      
      #Select from genes in gene.expr those that we know the hgnc symbol
      gene.expr <- exprs(x)[rownames(exprs(x))%in%annot.expr$id.feature,]
      #Replace gene ID for hgnc symbol
      rownames(gene.expr) <- annot.expr[, 'hgnc_symbol']
      #Select those genes that belong to chrY
      exprY <- gene.expr[annot.expr[, 'hgnc_symbol']%in%chrY$hgnc_symbol,]
      exprY <- exprY[complete.cases(exprY),]
      #Select those genes that belong to the rest of the genome
      exprRef <- gene.expr[annot.expr[, 'hgnc_symbol']%in%chrRef$hgnc_symbol,]
      exprRef <- exprRef[complete.cases(exprRef),]
      
  }
  else if (object.type == "rangedsummarizedexperiment"){
    
    #Filter males in the expression set
    if (!missing(gender.var)){
      ii <- which(names(colData(x)) == gender.var)
      x <- x[, !is.na(colData(x)[, ii])]
      x <- x[, colData(x)[, ii] == male.key]
    } else {warning("male.key not specified. Performing analysis with the whole dataset")}
    
    #Add column with the hgnc symbol information to fData
    if (gene.key != "hgnc_symbol"){
      rowData(x)$id.feature <- rownames(assay(x))
      #rowData(x) <- merge(EDY::annot, rowData(x), 
               #         by.x = "hgnc_symbol", 
                #        by.y = gene.key)
    }
    
    annot.expr <- data.frame(rowData(x))
    
    #Select from genes in gene.expr those that we know the hgnc symbol
    gene.expr <- assay(x)[rownames(assay(x))%in%annot.expr$id.feature,]
    #Replace gene ID for hgnc symbol
    rownames(gene.expr) <- annot.expr[, gene.key]
    #Select those genes that belong to chrY
    exprY <- gene.expr[annot.expr[, gene.key]%in%EDY::chrY$hgnc_symbol,]
    exprY <- exprY[complete.cases(exprY),]
    #Select those genes that belong to the rest of the genome
    exprRef <- gene.expr[annot.expr[, gene.key]%in%EDY::chrRef$hgnc_symbol,]
    exprRef <- exprRef[complete.cases(exprRef),]
    
  }
  
  #Apply EDY formulae: 
  if (experiment.type=="rnaseq"){
    
    if (log) {
      Ry <- sweep((exprY+1), 2, FUN="-", 
                  apply((exprRef+1), 2, mean))
      }else{ 
        Ry <- sweep(log2(exprY+1), 2, FUN="-", 
                  apply(log2(exprRef+1), 2, mean))}
    } else if (experiment.type == "microarray"){
        if (log) {
        Ry <- sweep(exprY, 2, FUN="-", 
                  apply(exprRef, 2, mean))
      }else{ 
        Ry <- sweep(log2(exprY), 2, FUN="-", 
                  apply(log2(exprRef), 2, mean))
      }
    }
  
  
  EDYcontinuous <- apply(Ry, 2, mean)
  
  if (!missing(group.var)&&!missing(control.key)){
    if (object.type == "expressionset"){
    controls <- EDYcontinuous[pData(x)[,group.var]==control.key]
    } else if (object.type == "rangedsummarizedexperiment"){
    controls <- EDYcontinuous[colData(x)[,group.var]==control.key]
    }
  } else {
    controls <- EDYcontinuous
    warning("No control group specified")
  }
  
  thresh <- median(controls, na.rm=TRUE) - coef*IQR(controls, na.rm=TRUE)
  EDY <- cut(EDYcontinuous, c(-Inf, thresh, Inf), 
             labels = c("Yes", "No"))
  EDY <- relevel(EDY, 2)
  
  #output
  if (object.type == "expressionset"){
    pData(x)$EDY <- EDY
    names(EDY) <- names(EDYcontinuous)
    ans <- list(EDY=EDY, EDYcontinuous=EDYcontinuous, 
                threshold=thresh, eSet=x)
  } else if (object.type == "rangedsummarizedexperiment"){
    colData(x)$EDY <- EDY
    names(EDY) <- names(EDYcontinuous)
    ans <- list(EDY=EDY, EDYcontinuous=EDYcontinuous, 
                threshold=thresh, RangedSummarizedExperiment=x)
  }
  class(ans) <- "EDY"
  ans
}




