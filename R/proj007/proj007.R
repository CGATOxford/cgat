speciesPlot <- function(dir=".", pattern="*testes-cap.replicated.cpg_density.export", main="CAPseq testes", xlab="CpGDensity", filename="test.pdf", plotcol=2, xlimit=c(0,5), ylimit=c(0.5))
{
    filelist <- dir(path=dir, pattern=pattern)
    leg = NULL
    pdf(file=filename, height=8, width=8, onefile=TRUE, family='Helvetica', paper='A4', pointsize=12, colormodel="cmyk")
    plot(0, xlim=xlimit, ylim=ylimit, xlab=xlab, main=main, type="n")
    for (i in 1:length(filelist)) {
        x <- read.table(file=paste(dir,filelist[i],sep="/"), na.string="None")
        data <- x[,plotcol]
        data <- data[data>0]
        d <- density(data, na.rm=TRUE)
        lines(d, col=i, lwd=3)
        med <- median(data,na.rm=T)
        abline(v=med, col=i)
        leg[i] <- substr(filelist[i],1,2)
    }
    legend("topright", legend=leg, col=seq(1:length(filelist)), lty=rep(1,length(filelist)), bty="n", lwd=rep(3,length(filelist)))
    dev.off()
}

combinedTSSPlot <- function(capseqfile="", nocapseqfile="", outfile="combined_tss_plot.pdf", ylimit=c(0,5), scale=1)
{
    capseq <- read.table(file=capseqfile, header=TRUE)
    # Adjust x axis
    capseq[,1] <- capseq[,1]-3000
    nocapseq <- read.table(file=nocapseqfile, header=TRUE)
    nocapseq[,1] <- nocapseq[,1]-3000
    ymax <- max(capseq[,2],capseq[,3],nocapseq[,2])*scale*1.1
    ylimit <- c(0,ymax)
    leg <- c("NMI TSS","NMI TTS","non-NMI TSS")
    pdf(file=outfile, height=8, width=8, onefile=TRUE, family='Helvetica', paper='A4', pointsize=12)
    plot(capseq[,1],capseq[,2]*scale,main="",xlab="Distance (bp)", ylab="Normalised Read depth", ylim=ylimit, lwd="3", col=1, type="l" )
    lines(capseq[,1], capseq[,3]*scale, lwd="3", col=2)
    lines(nocapseq[,1], nocapseq[,2]*scale, lwd="3", col=3)
    legend("topright", legend=leg, col=c(1,2,3), lty=c(1,1,1), lwd=c(3,3,3), bty="n")
    dev.off()}

sharesVsUniqueLengthPlot <- function(liver_shared="", liver_unique="", testes_shared="", testes_unique="", outfile="shared_vs_unique_interval_length.pdf", ylimit=c(0,1))
{
    liver_shared <- read.table(file=liver_shared, header=FALSE)
    testes_shared <- read.table(file=testes_shared, header=FALSE)
    liver_unique <- read.table(file=liver_unique, header=FALSE)
    testes_unique <- read.table(file=testes_unique, header=FALSE)
    d1 <- density(liver_shared[,1], na.rm=TRUE)
    d2 <- density(testes_shared[,1], na.rm=TRUE)
    d3 <- density(liver_unique[,1], na.rm=TRUE)
    d4 <- density(testes_unique[,1], na.rm=TRUE)
    ymax <- max(d1[[2]],d2[[2]],d3[[2]],d4[[2]])*1.1
    ylimit <- c(0,ymax)
    leg <- c("Liver shared","Testes shared", "Liver unique", "Testes unique")
    pdf(file=outfile, height=8, width=8, onefile=TRUE, family='Helvetica', paper='A4', pointsize=12, colormodel="cmyk")
    plot(d1, xlim=c(0,8000), ylim=ylimit, xlab="Length (bp)", main="", col="1", lwd=3)
    lines(d2, col=2, lwd=3)
    lines(d3, col=3, lwd=3)
    lines(d4, col=4, lwd=3)
    legend("topright", legend=leg, col=c(1,2,3,4), lty=rep(1,4), bty="n", lwd=rep(3,4))
    dev.off()
}

sharesVsUniqueCpgPlot <- function(liver_shared="", liver_unique="", testes_shared="", testes_unique="", outfile="shared_vs_unique_interval_cpg_obsexp.pdf", ylimit=c(0,1), xlimit=c(0,1.5), xlabel="CpG Observed/Expected")
{
    liver_shared <- read.table(file=liver_shared, header=FALSE)
    testes_shared <- read.table(file=testes_shared, header=FALSE)
    liver_unique <- read.table(file=liver_unique, header=FALSE)
    testes_unique <- read.table(file=testes_unique, header=FALSE)
    d1 <- density(liver_shared[,1], na.rm=TRUE)
    d2 <- density(testes_shared[,1], na.rm=TRUE)
    d3 <- density(liver_unique[,1], na.rm=TRUE)
    d4 <- density(testes_unique[,1], na.rm=TRUE)
    ymax <- max(d1[[2]],d2[[2]],d3[[2]],d4[[2]])*1.1
    ylimit <- c(0,ymax)
    leg <- c("Liver shared","Testes shared", "Liver unique", "Testes unique")
    pdf(file=outfile, height=8, width=8, onefile=TRUE, family='Helvetica', paper='A4', pointsize=12, colormodel="cmyk")
    plot(d1, ylim=ylimit, xlim=xlimit, xlab=xlabel, main="", col="1", lwd=3)
    lines(d2, col=2, lwd=3)
    lines(d3, col=3, lwd=3)
    lines(d4, col=4, lwd=3)
    legend("topright", legend=leg, col=c(1,2,3,4), lty=rep(1,4), bty="n", lwd=rep(3,4))
    dev.off()
}

liverTestesChromatinPlot <- function(infiles=c("",""), outfile="liver_testes_unique_chromatin.pdf")
{
    data <- list()
    length(data) <- length(infiles)
    x <- seq(1,3000)
    all <- NULL
    leg <- NULL
    for (i in 1:length(infiles)) {
        chromatin <- read.table(file=infiles[i], header=TRUE, stringsAsFactors=F)
        data[[i]] <- c(chromatin[,2],chromatin[,3],chromatin[,4])
        all <- c(all,chromatin[,2],chromatin[,3],chromatin[,4])
        pos1 <- gregexpr(".replicated.",infiles[i],fixed=T)[[1]][1]+12
        pos2 <- gregexpr(".profile.",infiles[i],fixed=T)[[1]][1]-1
        leg[i] <- substr(infiles[i], pos1, pos2)
    }
    ymax <- max(all)*1.1
    ylimit <- c(0,ymax)
    pdf(file=outfile, height=8, width=8, onefile=TRUE, family='Helvetica', paper='A4', pointsize=12)
    plot(x,data[[1]], col=2, lwd=3, ylim=ylimit, xlab="Unique Intervals", ylab="H3K4Me3", main="", type="l")
    for (i in 2:length(data)) {
        lines(x,data[[i]], col=i+1, lwd=3)
    }
    abline(v=c(1000,2000), lty=2,col=4)
    legend("topright", legend=leg, col=c(2,3), lty=c(1,1), lwd=c(3,3), bty="n")
    dev.off()
}

overlappedGenesProfilePlot <- function(overlapped="", control="", outfile="overlapped_genes_chromatin.pdf", ylabel="")
{
    overlapped <- read.table(file=overlapped, header=TRUE, stringsAsFactors=F)
    control <- read.table(file=control, header=TRUE, stringsAsFactors=F)
    ymax <- max(control[,2],control[,3],control[,4],overlapped[,2],overlapped[,3],overlapped[,4])
    ylimit <- c(0,ymax)
    pdf(file=outfile, height=8, width=8, onefile=TRUE, family='Helvetica', paper='A4', pointsize=12)
    plot(control[,1],control[,2], col=2, lwd=3, xlim=c(0,3000), ylim=ylimit, xlab="Gene Profile", ylab=ylabel, main="", type="l")
    lines(control[,1]+1000,control[,3],col=2, lwd=3)
    lines(control[,1]+2000,control[,4],col=2, lwd=3)
    lines(overlapped[,1],overlapped[,2], col=3, lwd=3)
    lines(overlapped[,1]+1000,overlapped[,3],col=3, lwd=3)
    lines(overlapped[,1]+2000,overlapped[,4],col=3, lwd=3)
    abline(v=c(1000,2000), lty=2,col=4)
    leg <- c("<10% Overlapped genes",">90% Overlapped genes")
    legend("topright", legend=leg, col=c(2,3), lty=c(1,1), lwd=c(3,3), bty="n")
    dev.off()
}

overlappedGenesSmoothedProfilePlot <- function(overlapped="", control="", outfile="overlapped_genes_smoothed_chromatin.pdf", ylabel="", smooth=0.7)
{
    overlapped <- read.table(file=overlapped, header=TRUE, stringsAsFactors=F)
    control <- read.table(file=control, header=TRUE, stringsAsFactors=F)
    ymax <- max(control[,2],control[,3],control[,4],overlapped[,2],overlapped[,3],overlapped[,4])
    ylimit <- c(0,ymax)
    #x_joined <- c(control[,1],control[,1]+length(control[,1]),control[,4]+*length(control[,1])*2))
    control_joined <- c(control[,2],control[,3],control[,4])
    control_smoothed <- smooth.spline(control_joined,spar=smooth)
    overlapped_joined <- c(overlapped[,2],overlapped[,3],overlapped[,4])
    overlapped_smoothed <- smooth.spline(overlapped_joined,spar=smooth)
    pdf(file=outfile, height=8, width=8, onefile=TRUE, family='Helvetica', paper='A4', pointsize=12)
    plot(control_smoothed, xlim=c(0,length(control_joined)), ylim=ylimit, xlab="Gene Profile", ylab=ylabel, main="", type="n")
    lines(control_smoothed,col=2, lwd=3)
    lines(overlapped_smoothed,col=3, lwd=3)
    abline(v=c(1000,2000), lty=2,col=4)
    leg <- c("<10% Overlapped genes",">90% Overlapped genes")
    legend("topright", legend=leg, col=c(2,3), lty=c(1,1), lwd=c(3,3), bty="n")
    dev.off()
}

twoWayVennPlot <- function(listsize1=1, list1name="A", listsize2=1, list2name="B", overlap=1, outfile="overlapped_genes_chromatin.pdf", ylabel="")
{
    library(VennDiagram) 
    list1 <- seq(1,listsize1)
    list2 <- seq(overlap,overlap+listsize2)
    x <- list(list1name=list1,list2name=list2)
    pdf(file=outfile, height=8, width=8, onefile=TRUE, family='Helvetica', paper='A4', pointsize=12)
    venn <- venn.diagram( x, filename=NULL, col="#58595B", fill=c("#EC1C24","#69BC45"), alpha=0.75, label.col=c("darkred", "white", "darkgreen"), cex=2.0, fontfamily="Helvetica", fontface="bold")
    grid.draw(venn)
    dev.off()
}
    
nmi_conservation <- function(infile="genelists_merged.stats", outfile="nmi_conservation") {
    # load file
    x <- read.table(file=infile,header=TRUE)
    # make contingency table
    species <- x[,1]
    non_conserved_genes <- x[,3]-x[,4]
    conserved_genes <- x[,4]
    nmi_genes <- x[,2]
    conserved_nmi_genes <- x[,5]
    nonconserved_nmi_genes <- nmi_genes-conserved_nmi_genes
    conserved_non_nmi_genes <- conserved_genes-conserved_nmi_genes
    nonconserved_non_nmi_genes <- non_conserved_genes-nonconserved_nmi_genes
    fisher_out <- NULL
    for (i in 1:length(species)) {
        mat <- matrix(c(conserved_nmi_genes[i],nonconserved_nmi_genes[i],conserved_non_nmi_genes[i],nonconserved_non_nmi_genes[i]),nrow=2,ncol=2)
        rownames(mat) <- c("conserved","nonconserved")
        colnames(mat) <- c("NMI","nonNMI")
        write.table(mat, file=paste(species[i],".",outfile,".contingency-table.tsv",sep=""), sep="\t")
        # Fisher exact test
        res <- fisher.test(mat)
        fisher_result <- data.frame(species[i],res$estimate,res$p.value,res$conf.int[1],res$conf.int[2])
        fisher_out <- rbind(fisher_out,fisher_result)
    }
    colnames(fisher_out) <- c("species","odds.ratio","p.value","conf.int.low","conf.int.high")
    write.table(fisher_out, file=paste(outfile,".fisher.test.tsv",sep=""), sep="\t", row.names=FALSE)
}

