# Read in data file 
x <- read.table(file="ac_liver-cap.transcript.tss-profile.tsv", header=TRUE)
leg <- c("TSS","TTS","TSS cap","TSS no cap")
x[,1] <- x[,1]-3000

y <- read.table(file="ac_liver-cap.transcript.tss-profile.capseq.tsv", header=TRUE)
y[,1] <- y[,1]-3000

z <- read.table(file="ac_liver-cap.transcript.tss-profile.nocapseq.tsv", header=TRUE)
z[,1] <- z[,1]-3000

pdf(file='ac_liver_tss-profile.pdf', height=6, width=6, onefile=TRUE, family='Helvetica', paper='A4', pointsize=12)

plot(x[,1],x[,2],main="",xlab="Distance (bp)", ylab="Read depth", ylim=c(0,8000000), lwd="3", type="n", )

lines(x[,1], x[,3], lwd="3", col=2)
lines(x[,1], x[,2], lwd="3", col=1)
lines(y[,1], y[,2], lwd="3", col=3)
lines(z[,1], z[,2], lwd="3", col=4)

legend("topright", legend=leg, col=c(1,2,3,4), lty=c(1,1,1,1), lwd=c(3,3,3,3), bty="n")
#dev.off

####################

# Read in data file 
x <- read.table(file="ac_testes-cap.transcript.tss-profile.tsv", header=TRUE)
leg <- c("TSS","TTS","TSS cap","TSS no cap")
x[,1] <- x[,1]-3000

y <- read.table(file="ac_testes-cap.transcript.tss-profile.capseq.tsv", header=TRUE)
y[,1] <- y[,1]-3000

z <- read.table(file="ac_testes-cap.transcript.tss-profile.nocapseq.tsv", header=TRUE)
z[,1] <- z[,1]-3000

pdf(file='ac_testes_tss-profile.pdf', height=6, width=6, onefile=TRUE, family='Helvetica', paper='A4', pointsize=12)

plot(x[,1],x[,2],main="",xlab="Distance (bp)", ylab="Read depth", ylim=c(0,8000000), lwd="3", type="n", )

lines(x[,1], x[,3], lwd="3", col=2)
lines(x[,1], x[,2], lwd="3", col=1)
lines(y[,1], y[,2], lwd="3", col=3)
lines(z[,1], z[,2], lwd="3", col=4)

legend("topright", legend=leg, col=c(1,2,3,4), lty=c(1,1,1,1), lwd=c(3,3,3,3), bty="n")
#dev.off

