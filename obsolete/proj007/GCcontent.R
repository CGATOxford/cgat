###########################################
# Read in data file and convert 'None' to NA

setwd("C:/Users/exet2865/Documents/CGAT007 Project/Capseq7/CpG and GC content/CpG O E")

#TESTES
#ac
x <- read.table(file="ac_testes-cap.replicated.cpg.export", na.string="None")
colnames(x) <- c("id","Capseq","Control","flank3","flank5")
ac <- x[,2:dim(x)[2]]

#dr
y <- read.table(file="dr_testes-cap.replicated.cpg.export", na.string="None")
colnames(y) <- c("id","Capseq","Control","flank3","flank5")
dr <- y[,2:dim(y)[2]]

#gg
z <- read.table(file="gg_testes-cap.replicated.cpg.export", na.string="None")
colnames(z) <- c("id","Capseq","Control","flank3","flank5")
gg <- z[,2:dim(z)[2]]

#hs
a <- read.table(file="hs_testes-cap.replicated.cpg.export", na.string="None")
colnames(a) <- c("id","Capseq","Control","flank3","flank5")
hs <- a[,2:dim(a)[2]]

#mm
b <- read.table(file="mm_testes-cap.replicated.cpg.export", na.string="None")
colnames(b) <- c("id","Capseq","Control","flank3","flank5")
mm <- b[,2:dim(b)[2]]

#oa
c <- read.table(file="oa_testes-cap.replicated.cpg.export", na.string="None")
colnames(c) <- c("id","Capseq","Control","flank3","flank5")
oa <- c[,2:dim(c)[2]]

#xt
d <- read.table(file="xt_testes-cap.replicated.cpg.export", na.string="None")
colnames(d) <- c("id","Capseq","Control","flank3","flank5")
xt <- d[,2:dim(d)[2]]

#d <- tapply(gc, 2, density, na.rm=TRUE)

d1 <- density(ac[,1], na.rm=TRUE)
d2 <- density(dr[,1], na.rm=TRUE)
d3 <- density(gg[,1], na.rm=TRUE)
d4 <- density(hs[,1], na.rm=TRUE)
d5 <- density(mm[,1], na.rm=TRUE)
d6 <- density(oa[,1], na.rm=TRUE)
d7 <- density(xt[,1], na.rm=TRUE)

leg <- c("ac","dr","gg","hs","mm","oa","xt")

pdf(file='HL_cpgobsexp_testes.pdf', height=6, width=6, onefile=TRUE, family='Helvetica', paper='A4', pointsize=12, colormodel='cmyk')

plot(d1, xlim=c(0,2), ylim=c(0,6), xlab="CpG obsexp", main="Testes CAPseq", col="1", lwd=3)

lines(d2, col=2, lwd=3)
lines(d3, col=3, lwd=3)
lines(d4, col=4, lwd=3)
lines(d5, col=5, lwd=3)
lines(d6, col=6, lwd=3)
lines(d7, col=7, lwd=3)

legend("topright", legend=leg, col=c(1,2,3,4,5,6,7), lty=rep(1,5), bty="n", lwd=rep(3,5))

dev.off()

###############################################################

#TESTES CONTROL
#ac
x <- read.table(file="ac_testes-cap.replicated.cpg.export", na.string="None")
colnames(x) <- c("id","Capseq","Control","flank3","flank5")
ac <- x[,2:dim(x)[2]]

#dr
y <- read.table(file="dr_testes-cap.replicated.cpg.export", na.string="None")
colnames(y) <- c("id","Capseq","Control","flank3","flank5")
dr <- y[,2:dim(y)[2]]

#gg
z <- read.table(file="gg_testes-cap.replicated.cpg.export", na.string="None")
colnames(z) <- c("id","Capseq","Control","flank3","flank5")
gg <- z[,2:dim(z)[2]]

#hs
a <- read.table(file="hs_testes-cap.replicated.cpg.export", na.string="None")
colnames(a) <- c("id","Capseq","Control","flank3","flank5")
hs <- a[,2:dim(a)[2]]

#mm
b <- read.table(file="mm_testes-cap.replicated.cpg.export", na.string="None")
colnames(b) <- c("id","Capseq","Control","flank3","flank5")
mm <- b[,2:dim(b)[2]]

#oa
c <- read.table(file="oa_testes-cap.replicated.cpg.export", na.string="None")
colnames(c) <- c("id","Capseq","Control","flank3","flank5")
oa <- c[,2:dim(c)[2]]

#xt
d <- read.table(file="xt_testes-cap.replicated.cpg.export", na.string="None")
colnames(d) <- c("id","Capseq","Control","flank3","flank5")
xt <- d[,2:dim(d)[2]]

#d <- tapply(gc, 2, density, na.rm=TRUE)

d1 <- density(ac[,2], na.rm=TRUE)
d2 <- density(dr[,2], na.rm=TRUE)
d3 <- density(gg[,2], na.rm=TRUE)
d4 <- density(hs[,2], na.rm=TRUE)
d5 <- density(mm[,2], na.rm=TRUE)
d6 <- density(oa[,2], na.rm=TRUE)
d7 <- density(xt[,2], na.rm=TRUE)

leg <- c("ac","dr","gg","hs","mm","oa","xt")

pdf(file='HL_cpgobsexp_testescontrol.pdf', height=6, width=6, onefile=TRUE, family='Helvetica', paper='A4', pointsize=12, colormodel='cmyk')

plot(d1, xlim=c(0,2), ylim=c(0,6), xlab="CpG obsexp", main="Testes Control CAPseq", col="1", lwd=3)

lines(d2, col=2, lwd=3)
lines(d3, col=3, lwd=3)
lines(d4, col=4, lwd=3)
lines(d5, col=5, lwd=3)
lines(d6, col=6, lwd=3)
lines(d7, col=7, lwd=3)

legend("topright", legend=leg, col=c(1,2,3,4,5,6,7), lty=rep(1,5), bty="n", lwd=rep(3,5))

dev.off()


###############################################################

#LIVER
#ac
x <- read.table(file="ac_liver-cap.replicated.cpg.export", na.string="None")
colnames(x) <- c("id","Capseq","Control","flank3","flank5")
ac <- x[,2:dim(x)[2]]

#dr
y <- read.table(file="dr_liver-cap.replicated.cpg.export", na.string="None")
colnames(y) <- c("id","Capseq","Control","flank3","flank5")
dr <- y[,2:dim(y)[2]]

#gg
z <- read.table(file="gg_liver-cap.replicated.cpg.export", na.string="None")
colnames(z) <- c("id","Capseq","Control","flank3","flank5")
gg <- z[,2:dim(z)[2]]

#hs
a <- read.table(file="hs_liver-cap.replicated.cpg.export", na.string="None")
colnames(a) <- c("id","Capseq","Control","flank3","flank5")
hs <- a[,2:dim(a)[2]]

#mm
b <- read.table(file="mm_liver-cap.replicated.cpg.export", na.string="None")
colnames(b) <- c("id","Capseq","Control","flank3","flank5")
mm <- b[,2:dim(b)[2]]

#oa
c <- read.table(file="oa_male_liver-cap.replicated.cpg.export", na.string="None")
colnames(c) <- c("id","Capseq","Control","flank3","flank5")
oa <- c[,2:dim(c)[2]]

#xt
d <- read.table(file="xt_liver-cap.replicated.cpg.export", na.string="None")
colnames(d) <- c("id","Capseq","Control","flank3","flank5")
xt <- d[,2:dim(d)[2]]

#d <- tapply(gc, 2, density, na.rm=TRUE)

d1 <- density(ac[,1], na.rm=TRUE)
d2 <- density(dr[,1], na.rm=TRUE)
d3 <- density(gg[,1], na.rm=TRUE)
d4 <- density(hs[,1], na.rm=TRUE)
d5 <- density(mm[,1], na.rm=TRUE)
d6 <- density(oa[,1], na.rm=TRUE)
d7 <- density(xt[,1], na.rm=TRUE)

leg <- c("ac","dr","gg","hs","mm","oa","xt")

pdf(file='HL_cpgobsexp_liver.pdf', height=6, width=6, onefile=TRUE, family='Helvetica', paper='A4', pointsize=12, colormodel='cmyk')

plot(d1, xlim=c(0,2), ylim=c(0,6), xlab="CpG obsexp", main="Liver CAPseq", col="1", lwd=3)

lines(d2, col=2, lwd=3)
lines(d3, col=3, lwd=3)
lines(d4, col=4, lwd=3)
lines(d5, col=5, lwd=3)
lines(d6, col=6, lwd=3)
lines(d7, col=7, lwd=3)

legend("topright", legend=leg, col=c(1,2,3,4,5,6,7), lty=rep(1,5), bty="n", lwd=rep(3,5))

dev.off()


###############################################################


#LIVER CONTROL
#ac
x <- read.table(file="ac_liver-cap.replicated.cpg.export", na.string="None")
colnames(x) <- c("id","Capseq","Control","flank3","flank5")
ac <- x[,2:dim(x)[2]]

#dr
y <- read.table(file="dr_liver-cap.replicated.cpg.export", na.string="None")
colnames(y) <- c("id","Capseq","Control","flank3","flank5")
dr <- y[,2:dim(y)[2]]

#gg
z <- read.table(file="gg_liver-cap.replicated.cpg.export", na.string="None")
colnames(z) <- c("id","Capseq","Control","flank3","flank5")
gg <- z[,2:dim(z)[2]]

#hs
a <- read.table(file="hs_liver-cap.replicated.cpg.export", na.string="None")
colnames(a) <- c("id","Capseq","Control","flank3","flank5")
hs <- a[,2:dim(a)[2]]

#mm
b <- read.table(file="mm_liver-cap.replicated.cpg.export", na.string="None")
colnames(b) <- c("id","Capseq","Control","flank3","flank5")
mm <- b[,2:dim(b)[2]]

#oa
c <- read.table(file="oa_male_liver-cap.replicated.cpg.export", na.string="None")
colnames(c) <- c("id","Capseq","Control","flank3","flank5")
oa <- c[,2:dim(c)[2]]

#xt
d <- read.table(file="xt_liver-cap.replicated.cpg.export", na.string="None")
colnames(d) <- c("id","Capseq","Control","flank3","flank5")
xt <- d[,2:dim(d)[2]]

#d <- tapply(gc, 2, density, na.rm=TRUE)

d1 <- density(ac[,2], na.rm=TRUE)
d2 <- density(dr[,2], na.rm=TRUE)
d3 <- density(gg[,2], na.rm=TRUE)
d4 <- density(hs[,2], na.rm=TRUE)
d5 <- density(mm[,2], na.rm=TRUE)
d6 <- density(oa[,2], na.rm=TRUE)
d7 <- density(xt[,2], na.rm=TRUE)

leg <- c("ac","dr","gg","hs","mm","oa","xt")

pdf(file='HL_cpgobsexp_livercontrol.pdf', height=6, width=6, onefile=TRUE, family='Helvetica', paper='A4', pointsize=12, colormodel='cmyk')

plot(d1, xlim=c(0,2), ylim=c(0,6), xlab="CpG obsexp", main="Liver Control CAPseq", col="1", lwd=3)

lines(d2, col=2, lwd=3)
lines(d3, col=3, lwd=3)
lines(d4, col=4, lwd=3)
lines(d5, col=5, lwd=3)
lines(d6, col=6, lwd=3)
lines(d7, col=7, lwd=3)

legend("topright", legend=leg, col=c(1,2,3,4,5,6,7), lty=rep(1,5), bty="n", lwd=rep(3,5))

dev.off()


###############################################################


#CGI
#ac
ac <- read.table(file="anoCar2.cgi.cpg.export", na.string="None")

#dr
dr <- read.table(file="danRer7.cgi.cpg.export", na.string="None")

#gg
gg <- read.table(file="galGal3.cgi.cpg.export", na.string="None")

#hs
hs <- read.table(file="hg19.cgi.cpg.export", na.string="None")

#mm
mm <- read.table(file="mm9.cgi.cpg.export", na.string="None")

#oa
oa <- read.table(file="ornAna1.cgi.cpg.export", na.string="None")

#xt
xt <- read.table(file="xenTro3.cgi.cpg.export", na.string="None")

#d <- tapply(gc, 2, density, na.rm=TRUE)

d1 <- density(ac[,2], na.rm=TRUE)
d2 <- density(dr[,2], na.rm=TRUE)
d3 <- density(gg[,2], na.rm=TRUE)
d4 <- density(hs[,2], na.rm=TRUE)
d5 <- density(mm[,2], na.rm=TRUE)
d6 <- density(oa[,2], na.rm=TRUE)
d7 <- density(xt[,2], na.rm=TRUE)

leg <- c("ac","dr","gg","hs","mm","oa","xt")

pdf(file='HL_cpgobsexp_CGI.pdf', height=6, width=6, onefile=TRUE, family='Helvetica', paper='A4', pointsize=12, colormodel='cmyk')

plot(d1, xlim=c(0,2), ylim=c(0,6), xlab="CpG obsexp", main="CGI", col="1", lwd=3)

lines(d2, col=2, lwd=3)
lines(d3, col=3, lwd=3)
lines(d4, col=4, lwd=3)
lines(d5, col=5, lwd=3)
lines(d6, col=6, lwd=3)
lines(d7, col=7, lwd=3)

legend("topright", legend=leg, col=c(1,2,3,4,5,6,7), lty=rep(1,5), bty="n", lwd=rep(3,5))

dev.off()


###############################################################


#TSS gene
#ac
ac <- read.table(file="anoCar2.tss.gene.cpg.export", na.string="None")

#dr
dr <- read.table(file="danRer7.tss.gene.cpg.export", na.string="None")

#gg
gg <- read.table(file="galGal3.tss.gene.cpg.export", na.string="None")

#hs
hs <- read.table(file="hg19.tss.gene.cpg.export", na.string="None")

#mm
mm <- read.table(file="mm9.tss.gene.cpg.export", na.string="None")

#oa
oa <- read.table(file="ornAna1.tss.gene.cpg.export", na.string="None")

#xt
xt <- read.table(file="xenTro3.tss.gene.cpg.export", na.string="None")

#d <- tapply(gc, 2, density, na.rm=TRUE)

d1 <- density(ac[,2], na.rm=TRUE)
d2 <- density(dr[,2], na.rm=TRUE)
d3 <- density(gg[,2], na.rm=TRUE)
d4 <- density(hs[,2], na.rm=TRUE)
d5 <- density(mm[,2], na.rm=TRUE)
d6 <- density(oa[,2], na.rm=TRUE)
d7 <- density(xt[,2], na.rm=TRUE)
med1 <- median(ac[,2], na.rm=TRUE)
med2 <- median(dr[,2], na.rm=TRUE)
med3 <- median(gg[,2], na.rm=TRUE)
med4 <- median(hs[,2], na.rm=TRUE)
med5 <- median(mm[,2], na.rm=TRUE)
med6 <- median(oa[,2], na.rm=TRUE)
med7 <- median(xt[,2], na.rm=TRUE)

leg <- c("ac","dr","gg","hs","mm","oa","xt")
#pdf(file='HL_cpgobsexp_TSSgene.pdf', height=6, width=6, onefile=TRUE, family='Helvetica', paper='A4', pointsize=12, colormodel='cmyk')
plot(d1, xlim=c(0,2), ylim=c(0,6), xlab="CpG obsexp", main="TSS gene", col="1", lwd=3)

lines(d2, col=2, lwd=3)
lines(d3, col=3, lwd=3)
lines(d4, col=4, lwd=3)
lines(d5, col=5, lwd=3)
lines(d6, col=6, lwd=3)
lines(d7, col=7, lwd=3)
abline(v=med1, col=1)
abline(v=med2, col=2)
abline(v=med3, col=3)
abline(v=med4, col=4)
abline(v=med5, col=5)
abline(v=med6, col=6)
abline(v=med7, col=7)

legend("topright", legend=leg, col=c(1,2,3,4,5,6,7), lty=rep(1,5), bty="n", lwd=rep(3,5))
dev.off()