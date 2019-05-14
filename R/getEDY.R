#' Detection of individuals with Extreme Downregulation of chromosome Y (EDY)
#' from transcriptomic data
#'
#' To assess EDY, the function measures the relative expression of the entire
#' chromosome Y with respect to the autosomes for each individual. For \emph{n}
#' probes (exons) in chromosome Y, with \eqn{x{i}} intensity (read counts) for
#' the \emph{i-th} probe, it computes \eqn{y = 1/n  \sumlog_2(x{i}) } as a
#' measure of the average expression of chromosome Y (summation limits between
#' \emph{i=1} and \emph{N}). Likewise, for \emph{m} probes in the autosomes, the
#' function computes the mean expression of autosomes \eqn{a = 1/m
#' \sumlog_2(x{i})} (summation limits between \emph{i=1} and \emph{M})
#' [\strong{NOTE:} for RNAseq data \eqn{log_2(x{i} + 1)} is computed to avoid
#' problems with zero counts]. The relative amount of an individual's Y
#' expression with respect to the individual's autosomes is then \eqn{Ry = y -
#' a}, and, in a population sample, the individual \emph{j} is considered with
#' EDY if \deqn{Ry{j} < median(Ry) - 1.2 · IQR(Ry)}
#'
#' where IQR is the inter-quartile range across the sample.
#'
#' @title getEDY
#' @param x An ExpressionSet or RangedSummarizedExperiment from microarray or
#'   RNAseq experiments.
#' @param gender.var A string indicating the name of the column in the table of
#'   phenotype (accessed by \code{pData} or \code{colData}) that contains the
#'   gender of the individuals.
#' @param male.key A string indicating the symbol that idendtifies males (e.g.,
#'   "male", "M", ...).
#' @param gene.key A string indicating the name of the column that contains the
#'   Gene Symbol in the table of features (accessed by \code{fData} or
#'   \code{rowData}).
#' @param coef Numerical. Value to consider an outlier when calling getEDY. Default
#'   correspond to 1.2 which keeps 5\% of data in the normal case.
#' @param log Logical. It is set to \code{TRUE} if the gene expression values
#'   are given in a logarithmic scale in the data set and \code{FALSE}
#'   otherwise. Default is \code{TRUE}.
#' @param group.var A string indicating the name of the column that contains the
#'   information about whether the individual is case or control (accessed by
#'   \code{pData} or \code{colData}).
#' @param control.key A string indicating the symbol that identifies control
#'   group (e.g., "control").
#' @param experiment.type A string indicating whether the data set is a
#'   \strong{microarray} or a \strong{RNAseq} experiment
#'   
#' @import Biobase
#' @import SummarizedExperiment
#' @export getEDY
#' @return A list containing 4 objects: 
#' \itemize{ 
#'   \item \code{EDY}: a factor with two levels, \code{NO} and \code{YES} , 
#'     that indicates whether the individual has EDY or not.
#'   \item \code{EDYcontinuous}: a vector with individual relative expression
#'     of chromosome Y with regard the autosomal genes.
#'   \item \code{threshold}: a number indicating the threshold from which down 
#'     an individual is considered to have EDY. 
#'   \item \code{eSet} or \code{RangedSummarizedExperiment}: an ExpressionSet 
#'     or RangedSummarizedExperiment that contains only males from the initial
#'     data set. It contains a new column in \code{pData} or \code{colData}
#'     containing the EDY status of each individual.
#'   }

getEDY <- function(x, gender.var, male.key, gene.key, coef=1.2, 
                   log = TRUE, group.var, control.key, experiment.type, ...){
  
  object.type <- tolower(class(x)[1])
  experiment.type <- tolower(experiment.type)
  
  if (experiment.type!="microarray" && experiment.type!="rnaseq"){
    stop("Invalid experiment.type. Allowed types are 'microarray' or 'RNAseq'")
  }
  
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
        }
      
      annot.expr <- fData(x)
      
      #Select from genes in gene.expr those that we know the hgnc symbol
      gene.expr <- exprs(x)[rownames(exprs(x))%in%annot.expr$id.feature,]
      #Replace gene ID for hgnc symbol
      rownames(gene.expr) <- annot.expr[, gene.key]
      #Select those genes that belong to chrY
      exprY <- gene.expr[annot.expr[, gene.key]%in%EDY::chrY.genes$hgnc_symbol,]
      exprY <- exprY[complete.cases(exprY),]
      #Select those genes that belong to the rest of the genome
      exprRef <- gene.expr[annot.expr[, gene.key]%in%EDY::autosomal.genes$hgnc_symbol,]
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
    }
    
    annot.expr <- data.frame(rowData(x))
    
    #Select from genes in gene.expr those that we know the hgnc symbol
    gene.expr <- assay(x)[rownames(assay(x))%in%annot.expr$id.feature,]
    #Replace gene ID for hgnc symbol
    rownames(gene.expr) <- annot.expr[, gene.key]
    #Select those genes that belong to chrY
    exprY <- gene.expr[annot.expr[, gene.key]%in%EDY::chrY.genes$hgnc_symbol,]
    exprY <- exprY[complete.cases(exprY),]
    #Select those genes that belong to the rest of the genome
    exprRef <- gene.expr[annot.expr[, gene.key]%in%EDY::autosomal.genes$hgnc_symbol,]
    exprRef <- exprRef[complete.cases(exprRef),]
    
  }
  
  #Apply EDY formulae: 
  if (experiment.type=="rnaseq"){
    
    if (log) {
      Ry <- sweep((exprY), 2, FUN="-", 
                  apply((exprRef), 2, mean))
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




