##' Obtaining hgnc symbol and entrezgene from genBank accession number
##'
##' Get HUGO gene nomenclature committee (hgnc) symbol from genBank accession
##' number, as well as the start and end position of the gene, the chromosome
##' and the entrezgene (NCBI gene ID)
##' @title genBank AN to hgnc symbol
##' @param x The table of features of an ExpressionSet or any data.frame or
##'   matrix.
##' @param genBank.col The name of the column that contains the genBank
##'   accession numbers.
##' 
##' @import biomaRt
##' @import org.Hs.eg.db
##' @export GB_to_hgnc
##' @return A \code{data.frame} containing the previous information in \code{x} plus 5
##'   more columns: \code{start_position} (of the gene), \code{end_position} (of
##'   the gene), \code{chromosome_name}, \code{hgnc_symbol} and
##'   \code{entrezgene}

GB_to_hgnc <- function(x, genBank.col, ...){
  
  list_entrez_GB <- as.list(org.Hs.egACCNUM2EG)
  query <- x[, genBank.col]
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
  x <-  merge(matrix_entr_GB, x, by.x = "GB_ACC", by.y = genBank.col)
  
  # Get hgnc_symbol from entrezgene
  ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")

  hgnc_symbols <- getBM(attributes = c("start_position", "end_position", "chromosome_name", 
                                       "hgnc_symbol", "entrezgene"), filters = "entrezgene",
                        values = entrezgenes, mart = ensembl)

  #Join to fData
  x <- merge(hgnc_symbols, x, by.x = "entrezgene", by.y = "entrezgene")
  #Delete duplicated GB_ACCs
  x <- x[!duplicated(x[,"GB_ACC"]),]
  
  #output
  x
}
