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
  
  fData(x)$id.feature <- featureNames(x)
  query <- fData(x)[, key.col]
  query <- query[!is.na(query)]
  key.type <- toupper(key.type)
  key.types <- keytypes(EnsDb.Hsapiens.v86)
  if (!(key.type%in%key.types) && key.type!="GENBANK"){
    stop("Invalid key.id. Allowed choices are: 'ENTREZID', 'EXONID', 'GENEBIOTYPE', 'GENEID', 'GENENAME', 'PROTDOMID', 'PROTEINDOMAINID', 'PROTEINDOMAINSOURCE', 'PROTEINID', 'SEQNAME', 'SEQSTRAND', 'SYMBOL', 'TXBIOTYPE', 'TXID', 'TXNAME', 'UNIPROTID' and 'GENBANK'")
  }
  #Get entrezgene from id
  else if (key.type%in%key.types){
    hgnc_symbols <- select(EnsDb.Hsapiens.v86, keys = query, keytype = key.type,
           columns = c(key.type, "SYMBOL"))
    #Join to fData
    fData(x) <- merge(hgnc_symbols, fData(x), by.x = key.type, by.y = key.col)
    number <- which(names(fData(x))=="SYMBOL")
    names(fData(x))[number] <- "hgnc_symbol"
    #Delete duplicated
    fData(x) <- fData(x)[!duplicated(fData(x)[, key.type]),]
    }
  else if (key.type=="GENBANK"){ 
      list_entrez_id <- as.list(org.Hs.egACCNUM2EG) 
  
      id_entrezgenes <- list_entrez_id[query]
      id_entrezgenes <- id_entrezgenes[!is.na(names(id_entrezgenes))]
    
      id <- names(id_entrezgenes)
    
      entrezgenes <- c()
      for (i in 1:length(id_entrezgenes)){
        entrezgenes[i] <- id_entrezgenes[[i]]
      }
      
      matrix_entr_id <- matrix(c(entrezgenes, id), ncol = 2)
      colnames(matrix_entr_id) <- c("ENTREZID", key.col)
      
      hgnc_symbols <- select(EnsDb.Hsapiens.v86, keys = entrezgenes, keytype = "ENTREZID",
                             columns = c("ENTREZID", "SYMBOL"))
      
      #Join to previous matrix
      matrix_entr_id <- merge(hgnc_symbols, matrix_entr_id, by.x = "ENTREZID", by.y ="ENTREZID")
      #Join to fData
      fData(x) <- merge(matrix_entr_id, fData(x), by.x = key.col, by.y = key.col)
      number <- which(names(fData(x))=="SYMBOL")
      names(fData(x))[number] <- "hgnc_symbol"
      #Delete duplicated
      fData(x) <- fData(x)[!duplicated(fData(x)[, "ENTREZID"]),]
  }
  
  #output
  x
}

