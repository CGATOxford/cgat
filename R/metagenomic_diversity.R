####################################################
####################################################
####################################################
# Given a data frame of taxomonc counts perform
# rarefaction analysis 
####################################################
####################################################
####################################################

library(plyr)
library(vegan)
library(gtools)
library(scales)

rarefaction <- function(dat, from = 1, to = 5000000, step = 1000000, groups = c()){
  
            # data frame with samples as rows and taxa
	    # as columns
	    sums <- apply(dat, 1, sum)
            ngroups <- length(unique(groups))
          
            # the result container
            result <- data.frame("group" = c(), "mean" = c(), "se" = c(), "sample" = c())

            # calculate diversity at each step up to "to"
	    for (i in seq(from, to, step)){
                 r <- data.frame(rarefy(dat, sample = i))
                 colnames(r) <- "richness"
                 r$group <- groups                 

                 # calculate mean and se per group
                 res <- ddply(r, .(group), summarize, mean = mean(richness), se = sd(richness)/sqrt(length(richness)))
                 res$sample <- factor(i)
                 result <- rbind(result, res)	
            }
            return (result)
}

processRarefactionResults <- function(dat, rf){

	# process rarefaction for samples that don't hit the max counts
	sums <- data.frame(colSums(dat))

	for (i in 1:nrow(rf)){
	    if (sums[rf[i,]$group,] < as.numeric(as.character(rf[i,]$sample))){
	        rf[i,]$mean <- NA
	    }
	}
	write.table(rf, file="test.tsv",sep="\t", row.names=F, quote=F)
	return (rf)
}

plotRarefaction <- function(rf, colours = c("brown", "darkGreen", "slateGrey", "darkBlue")){
           
          # plot rarefaction curve
          library(ggplot2)
	  rf$group <- as.character(rf$group)
	  rf$group <- factor(rf$group, levels=mixedsort(unique(rf$group)))
          plot1 <- ggplot(rf, aes(x = as.numeric(as.character(sample)), y = mean, colour=group, group = group))
          plot2 <- plot1 + geom_line() + geom_errorbar(aes(ymax = mean + se, ymin = mean - se), width = 0.25)
          plot2 + scale_colour_manual(values = colours) + scale_x_continuous(labels = comma) + theme(axis.text.x=element_text(angle=90))
}


plotSpecaccum <- function(dat){
    
          # plot accumulation of species as more samples
          # are added
          s <- specaccum(dat, method = "random")
          s <- s2 <- data.frame(s$sites, s$richness, s$sd)
          plot = ggplot(s2, aes(x = s.sites, y = s.richness)) + geom_line() + geom_errorbar(aes(ymax = s.richness + s.sd, ymin = s.richness - s.sd), width = 0.25)
          plot + scale_x_discrete(labels = comma) + theme(axis.text.x=element_text(angle=90))
}


buildDiversity <- function(dat, outfile, index = "shannon"){
	       d <- data.frame(diversity(dat, index = index))
	       d$sample <- rownames(d)
	       colnames(d)[1] <- index
	       d <- d[,c("sample", index)]
	       write.table(d, file=outfile, sep="\t", quote=F, row.names=F)
}


plotDiversity <- function(dat, index = "shannon", colours = c("brown", "darkGreen", "slateGrey", "darkBlue"), groups = c()){

          # plot shannon diversity index
          d <- data.frame(diversity(dat, index = index))
          colnames(d) <- "diversity"

          d$groups <- groups
          d <- ddply(d, .(groups), summarize, mean = mean(diversity), se = sd(diversity)/sqrt(length(diversity)))
	  d$groups <- as.character(d$groups)
	  d$groups <- factor(d$groups, levels=mixedsort(unique(d$groups)))
	  plot <- ggplot(d, aes(x=groups, y=mean, fill=groups)) + geom_bar(stat = "identity", position = "dodge") 
          plot + geom_errorbar(aes(ymax = mean + se, ymin = mean - se), width = 0.25) + scale_fill_manual(values = colours)
}

div.test <- function(dat, method = "kruskal", groups = c()){

         # test significance of diversity estimates
          d <- data.frame(diversity(dat))
          colnames(d) <- "diversity"
	  d$groups <- groups
          k <- kruskal.test(d$diversity, d$groups)
          k
}

richness.test <- function(dat, sample, groups = c()){

                 # test significance of richness at
                 # highest sampling
                 # rownames are taxa
                 r <- rarefy(dat, sample)
                 k <- kruskal.test(r~factor(groups))
                 k
		 }





