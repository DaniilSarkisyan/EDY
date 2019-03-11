library('biomaRt')
library('org.Hs.eg.db')

##' Obtaining hgnc symbol and entrezgene from genBank accession number
##'
##' Get HUGO gene nomenclature committee (hgnc) symbol from genBank accession
##' number in an ExpressionSet table of features
##' @title genBank AN to hgnc symbol
##' @param x The table of features of an ExpressionSet or any data.frame.
##' @param genBank.col The name of the column that contains the genBank
##'   accession numbers.
##'
##' @export GB_to_hgnc
##' @return A data.frame containing the previous information in \code{x} plus 5
##'   more columns: \code{start_position} (of the gene), \code{end_position} (of
##'   the gene), \code{chromosome_name}, \code{hgnc_symbol} and
##'   \code{entrezgene}

GB_to_hgnc <- function(fdata, genBank.col, ...){
  
  list_entrez_GB <- as.list(org.Hs.egACCNUM2EG)
  query <- fdata[, genBank.col]
  GB_entrezgenes <- list_entrez_GB[query]
  GB_entrezgenes <- GB_entrezgenes[!is.na(names(GB_entrezgenes))]

  GB_ACC <- names(GB_entrezgenes)

  entrezgenes <- c()
  for (i in 1:length(GB_entrezgenes)){
    entrezgenes[i] <- GB_entrezgenes[[i]]
  }
  
  matrix_entr_GB <- matrix(c(entrezgenes, GB_ACC), ncol = 2)
  colnames(matrix_entr_GB) <- c("entrezgene", "GB_ACC")
  
  #Join to fData
  fdata <-  merge(matrix_entr_GB, fdata, by.x = "GB_ACC", by.y = "GB_ACC")
  
  # Get hgnc_symbol from entrezgene
  ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")

  hgnc_symbols <- getBM(attributes = c("start_position", "end_position", "chromosome_name", 
                                       "hgnc_symbol", "entrezgene"), filters = "entrezgene",
                        values = entrezgenes, mart = ensembl)

  #Join to fData
  fdata <- merge(hgnc_symbols, fdata, by.x = "entrezgene", by.y = "entrezgene")
  #Delete duplicated GB_ACCs
  fdata <- fdata[!duplicated(fdata[,"GB_ACC"]),]
  
  #output
  fdata
}
