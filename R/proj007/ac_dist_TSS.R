#Plot interval distance from TSS
#Read table
x <- read.table(file="ac_liver-cap.replicated.gene.tss", header=TRUE, stringsAsFactors=F)

x2 <- as.numeric(x[,3])
x3 <- x2[!is.na(x2)]
d = density(x3)

#Read table
t <- read.table(file="ac_testes-cap.replicated.gene.tss", header=TRUE, stringsAsFactors=F)

t2 <- as.numeric(t[,3])
t3 <- t2[!is.na(t2)]
d2 = density(t3)

pdf(file='ac_liver_testes_tss_distance.pdf', height=6, width=6, onefile=TRUE, family='Helvetica', paper='A4', pointsize=12)

plot(d2,xlim=c(0,500000), lwd=3)

lines(d,col=3, lwd=3)

leg <- c("liver","testes")
legend("topright", legend=leg, col=c(1,3), lty=c(1,1), lwd=c(3,3), bty="n")

dev.off()