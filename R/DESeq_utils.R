###############################################################
# R functions for plotting using the new version of deseq 
# taken from the vignette
###############################################################


plotDispEsts <- function( cds , main = "")
{

plot(
rowMeans( counts( cds, normalized=TRUE ) ),
fitInfo(cds)$perGeneDispEsts,
main = main, pch = 16, log="xy", cex.main = 4, cex = 2, cex.lab = 4, cex.axis = 4, mgp = c(6,0,0), mkh = 100)
xg <- 10^seq( -.5, 5, length.out=300 )
lines( xg, fitInfo(cds)$dispFun( xg ), col="red" , lwd = 3)
}



plotDE <- function( res, main = "" ){
plot(
res$baseMean,
res$log2FoldChange,
log="x", pch=16, , main = main,cex.main = 4, cex = 2, cex.lab = 4, cex.axis = 4, mgp = c(6,0,0), mkh = 100
,col = ifelse( res$padj < .05, "red", "black" ) )
}



