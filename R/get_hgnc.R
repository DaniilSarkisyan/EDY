#' Obtaining hgnc symbol.
#'
#' Get HUGO gene nomenclature committee (hgnc) symbol as well as the start and
#' end position of the gene and the chromosome name from genBank accession
#' numbers, entrezgenes (NCBI gene IDs) or Ensembl stable IDs.
#' @title Get hgnc symbol
#' @param x The table of features of an ExpressionSet (accessed by
#'   \code{fData()}) or any data.frame or matrix.
#' @param gene.id The type of gene ID that we want to get hgnc symbol from. It
#'   must be one of \code{genBank}, \code{entrezgene} or \code{ensembl}.
#' @param gene.col The name of the column that contains the \code{gene.id}.
#' 
#' @importFrom biomaRt useMart getBM
#' @import org.Hs.eg.db
#' @export get_hgnc 
#' @return A \code{data.frame} containing the previous information in \code{x}
#'   plus 5 more columns: \code{start_position} (of the gene),
#'   \code{end_position} (of the gene), \code{chromosome_name},
#'   \code{hgnc_symbol} and \code{entrezgene}.

get_hgnc <- function(x, gene.id, gene.col, ...){
  
  query <- x[, gene.col]
  type <- charmatch(tolower(gene.id), c("entrezgene", "genbank", "ensembl"))
  if (is.na(type) || type==0){
    stop("Invalid gene.id. Try one from 'entrezgene', 'genbank' or 'ensembl'")
  }
  #Get entrezgene from id
  else if (type == 2 || type == 3){
    if (type == 2){
      list_entrez_id <- as.list(org.Hs.egACCNUM2EG)
    }
    else { 
      list_entrez_id <- as.list(org.Hs.egENSEMBL2EG) 
    }
    id_entrezgenes <- list_entrez_id[query]
    id_entrezgenes <- id_entrezgenes[!is.na(names(id_entrezgenes))]
  
    id <- names(id_entrezgenes)
  
    entrezgenes <- c()
    for (i in 1:length(id_entrezgenes)){
      entrezgenes[i] <- id_entrezgenes[[i]]
    }
    
    matrix_entr_id <- matrix(c(entrezgenes, id), ncol = 2)
    colnames(matrix_entr_id) <- c("entrezgene", gene.col)
    
    #Join to fData
    x <-  merge(matrix_entr_id, x, by.x = gene.col, by.y = gene.col)
  } 
  else if (type == 1) {
    entrezgenes <- query
    position <- which(names(x)==gene.col)
    names(x)[position] <- "entrezgene"
    gene.col <- "entrezgene"
  }
  # Get hgnc_symbol from entrezgene
  ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

  hgnc_symbols <- getBM(attributes = c("start_position", "end_position", "chromosome_name", 
                                       "hgnc_symbol", "entrezgene"), filters = "entrezgene",
                        values = entrezgenes, mart = ensembl)

  #Join to fData
  x <- merge(hgnc_symbols, x, by.x = "entrezgene", by.y = "entrezgene")
  #Delete duplicated ids
  x <- x[!duplicated(x[,gene.col]),]
  
  #output
  x
}

