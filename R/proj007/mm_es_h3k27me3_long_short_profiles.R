#Read table
ek27_lg <- read.table(file="mm_es-cap.replicated.long.genes.mm_es_H3K27Me3.profile.tsv.gz", header=TRUE, stringsAsFactors=F)

ek27_l <- read.table(file="mm_es-cap.replicated.long.mm_es_H3K27Me3.profile.tsv.gz", header=TRUE, stringsAsFactors=F)

ek27_s <- read.table(file="mm_es-cap.replicated.short.mm_es_H3K27Me3.profile.tsv.gz", header=TRUE, stringsAsFactors=F)

pdf(file='mm_escNMI_escH3K4me3_longgene_long_short.pdf', height=6, width=6, onefile=TRUE, family='Helvetica', paper='A4', pointsize=12)

plot(ek27_lg[,1],ek27_lg[,2], col=1, lwd=3, xlim=c(0,3000), ylim=c(0,1000), xlab="Distance", ylab="H3K4Me3", main="ESC NMI Interval")
points(ek27_lg[,1]+1000,ek27_lg[,3],col=1, lwd=3)
points(ek27_lg[,1]+2000,ek27_lg[,4],col=1, lwd=3)

points(ek27_l[,1],ek27_l[,2], col=2, lwd=3)
points(ek27_l[,1]+1000,ek27_l[,3],col=2, lwd=3)
points(ek27_l[,1]+2000,ek27_l[,4],col=2, lwd=3)

points(ek27_s[,1],ek27_s[,2], col=3, lwd=3)
points(ek27_s[,1]+1000,ek27_s[,3],col=3, lwd=3)
points(ek27_s[,1]+2000,ek27_s[,4],col=3, lwd=3)

abline(v=c(1000,2000), lty=2,col=4)
leg <- c("long_gene","long","short")
legend("topright", legend=leg, col=c(1,2,3), lty=c(1,1,1), lwd=c(3,3,3), bty="n")
dev.off()

#pdf(file='mm_testesNMI_livertestesH3K4me3.pdf', height=6, width=6, onefile=TRUE, family='Helvetica', paper='A4', pointsize=12)

plot(tl[,1],tl[,2], col=1, lwd=3, xlim=c(0,3000), ylim=c(0,35000),xlab="Distance", ylab="H3K4Me3", main="Testes NMI Interval")
points(tl[,1]+1000,tl[,3],col=1, lwd=3)
points(tl[,1]+2000,tl[,4],col=1, lwd=3)
points(tt[,1],tt[,2], col=2, lwd=3)
points(tt[,1]+1000,tt[,3],col=2, lwd=3)
points(tt[,1]+2000,tt[,4],col=2, lwd=3)
abline(v=c(1000,2000), lty=2,col=4)
leg <- c("liver","testes")
legend("topright", legend=leg, col=c(1,2), lty=c(1,1), lwd=c(3,3), bty="n")
dev.off()

