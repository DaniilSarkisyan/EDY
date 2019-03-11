library('biomaRt')
library('org.Hs.eg.db')

##' Obtaining hgnc symbol from genBank accession number
##'
##' Obtaining HUGO gene nomenclature committee (hgnc) symbol from genBank accession number
##' @title GB_ACC to hgnc_symbol
##' @param x 
##' @param 
##' 
##' @export getEDY
##' @return ...

GB_to_hgnc <- function(fdata, GB.column, ...){
  
  list_entrez_GB <- as.list(org.Hs.egACCNUM2EG)
  query <- fdata[,GB.column]
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
                                       "hgnc_symbol", "entrezgene", ), filters = "entrezgene",
                        values = entrezgenes, mart = ensembl)

  #Join to fData
  fdata <- merge(hgnc_symbols, fdata, by.x = "entrezgene", by.y = "entrezgene")
  #Delete duplicated GB_ACCs
  fdata <- fdata[!duplicated(fdata[,"GB_ACC"]),]
  
  #output
  fdata
}
