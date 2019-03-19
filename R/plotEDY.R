plot.EDY <- function(x, ...){
  edy <- apply(x$Ry, 1, mean)
  plot(edy, xlab = "Individuals", ylab = "chrY / autosomomes",
       type="n")
  mycol <- ifelse(x$EDY=="No", "blue", "red")
  points(edy, pch=16, col=mycol)
  legend("bottomright", c("normal", "EDY"), pch=16,
         col=c("blue", "red"))
}
