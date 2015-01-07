wiggle2Hilbert <- function(wiggleRle, image.dir, datName) {
  
  # Requires an rle object imported using rtracklayer::import.bw - v.slow and 
  # lots of memory required for bp resolution bigwig files.
  # remove non-canonical contigs
  # generate hilbert curves and plot as png - save in image.dir
  
  suppressPackageStartupMessages(library(HilbertVis))
    
  contigs2keep <- names(wiggleRle)

  for(contig in contigs2keep){
    # convert to matrix that represents hilbert curve
    hilb <- hilbertImage(wiggleRle[[contig]])
    
    # save and plot in a png
    fileSave <- paste0(image.dir,"/", contig,"-",datName , "-hilbertImage.png")
    hilbImage <- showHilbertImage(hilb,
    	      	 palettePos=colorRampPalette(c("white", "red"))(300),
		 paletteNeg=colorRampPalette(c("white", "blue"))(300),
		 maxPaletteValue=max(abs(hilb)))
    png(file=fileSave)
    print(hilbImage)
    dev.off()
  }
}

# stolen from http://www.ebi.ac.uk/huber-srv/hilbert/gallery.html
# reads in gff to use with hilbert Vis - thanks to Simon Anders


readGFF <- function( file, chr ) {

   # read in a GFF file:
   a <- read.table( file, as.is=TRUE )

   # select the desired chromsome by comparing with column 1:
   b <- subset( a, V1==chr )

   # make a wiggle vector with the scores
   w <- makeWiggleVector( b$V4, b$V5, b$V6, max(b$V5) + 10 )

   # make another wiggle vector, this time ignoring the score and simply
   # putting a 1 everywhere. This leaves 0 at all base pairs not
   # interrogated by probes:
   wm <- makeWiggleVector( b$V4, b$V5, rep( 1, nrow(b) ), max(b$V5)+10 )

   # Replace all the non-interrogated base pairs with NA
   w[ wm == 0 ] <- NA
   
   # return this
   return(w)
}
