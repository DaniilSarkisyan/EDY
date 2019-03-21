print.EDY <- function(x, ...){
  cat("EDY object\nContains information about males only from the original ExperessionSet\n$Ry: matrix. Relative expression of chrY genes.\n")
  cat(dim(x$Ry))
  cat("\n$EDY:\n$threshold: ")
  cat(x$threshold)
  cat("\n$eSet:\n")
  cat(x$eSet)
}
