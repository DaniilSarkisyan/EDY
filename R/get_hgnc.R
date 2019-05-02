#' Obtaining hgnc symbol.
#'
#' Get HUGO gene nomenclature committee (hgnc) symbol as well as the start and
#' end position of the gene and the chromosome name from genBank accession
#' numbers, entrezgenes (NCBI gene IDs) or Ensembl stable IDs.
#' @title Get hgnc symbol
#' @param x An ExpressionSet.
#' @param key.type A string indicating the type of gene ID that we want to get
#'   hgnc symbol from. It must be one of \code{'ENTREZID', 'EXONID',
#'   'GENEBIOTYPE', 'GENEID', 'GENENAME', 'PROTDOMID', 'PROTEINDOMAINID',
#'   'PROTEINDOMAINSOURCE', 'PROTEINID', 'SEQNAME', 'SEQSTRAND', 'SYMBOL',
#'   'TXBIOTYPE', 'TXID', 'TXNAME', 'UNIPROTID' and 'GENBANK'}.
#' @param key.col A string indicating the name of the column that contains the
#'   \code{key.type} in the table of features of the ExpressionSet (accessed by
#'   \code{fData}).
#' 
#' @import org.Hs.eg.db
#' @import EnsDb.Hsapiens.v86
#' @export get_hgnc 
#' @return An \code{ExpressionSet} containing the previous information in
#'   \code{x} plus the column \code{hgnc_symbol} containing the gene symbols.

get_hgnc <- function(x, key.type, key.col, ...){
  
  object.type <- class(x)[1]
  if (object.type == "ExpressionSet"){
    fData(x)$id.feature <- featureNames(x)
    query <- fData(x)[, key.col]
  } else if (object.type == "RangedSummarizedExperiment"){
    colData(x)$id.feature <- rownames(assay(x))
    query <- colData(x)[, key.col]
  }
  
  key.type <- tolower(key.type)
  key.types <- c("refseq", "uniprot", "ensembl", "entrezgene" )
  if (!(key.type%in%key.types) && key.type!="genbank"){
    stop("Invalid key.type. Allowed choices are: 'refseq', 'uniprot', 'ensembl', 'entrezgene' 
         and 'genbank'")
  }
  #Get entrezgene from id
  else if (key.type%in%key.types){
    dictionary <- EDY::genome.annot$hgnc_symbol
    names(dictionary) <- EDY::genome.annot[, key.type]
    hgnc_symbol <- dictionary[query]
    hgnc_symbol <- unname(hgnc_symbol)
    #Join to fData
    if (object.type == "ExpressionSet"){
      fData(x) <- cbind(hgnc_symbol, fData(x))
    } else if (object.type == "RangedSummarizedExperiment"){
      colData(x) <- cbind(hgnc_symbol, colData(x))
      }
    }
  else if (key.type=="genbank"){ 
      list_entrez_id <- as.list(org.Hs.egACCNUM2EG) 
      
      GB.ids <- names(list_entrez_id)  
      dictionary <- unlist(list_entrez_id)
      names(dictionary) <- GB.ids
    
      entrezgenes <- unname(dictionary[query])
      
      dictionary2 <- EDY::genome.annot$hgnc_symbol
      names(dictionary2) <- EDY::genome.annot$entrezgene
      
      hgnc_symbol <- dictionary2[entrezgenes]
     
      #Join to fData
      if (object.type == "ExpressionSet"){
        fData(x) <- cbind(hgnc_symbol, fData(x))
        } else if (object.type == "RangedSummarizedExperiment"){
        colData(x) <- cbind(hgnc_symbol, colData(x))
        }
  }
  
  #output
  x
}

