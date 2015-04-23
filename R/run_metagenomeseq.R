###################################################
###################################################
# metagenomeSeq tests for differences in abundance
# of features across two or more conditions. It 
# requires R version 3 so this script performs
# normalisation and analysis without the neeed for
# rpy2
###################################################
###################################################

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("metagenomeSeq"))

# make options list
option_list <- list(
               make_option(c("-c", "--count-matrix"),
                           help="input count matrix"),
               make_option(c("-p", "--prefix"),
                           help="prefix to use for outfiles"),
               make_option(c("--k"), type="integer",
                           help="number of samples over -a"),
               make_option(c("--a"), type="double",
                           help="number at which k samples should be 
                           over to be tested")
               )

# suppress warning messages
options(warn=-1)


###########################
# get command line options
###########################
opt <- parse_args(OptionParser(option_list=option_list))

#######################
# read count matrix in
#######################
print ("reading count matrix")
countMatrix <- read.csv(opt$`count-matrix`, 
	                header = T,
                        stringsAsFactors = F,
                        sep = "\t",
                        row.names = 1)


# get reads per million
rpm <- sweep(countMatrix, 2, colSums(countMatrix)/1000000, "/")

######################
# filter data
######################
print ("filtering data by RPM")
tokeep <- which(rowSums(rpm >= opt$a) >= opt$k)
countMatrix <- countMatrix[tokeep, ]

# get taxa
taxa = rownames(countMatrix)

# create list of data
input.data <- list(taxa = data.frame(taxa), 
	           counts = data.frame(countMatrix))

############################################
############################################
# get samples and associated conditions
# at present we are only considering 
# comparisons between "internal" conditions
# TODO: extend to all pairwise comparisons
############################################
############################################

print ("assigning samples to conditions")
samples = colnames(countMatrix)
samples = data.frame(samples)
samples$samples <- as.character(samples$samples)
rownames(samples) <- samples$samples
conditions <- unlist(strsplit(samples$samples, ".R[0-9]"))[seq(1, nrow(samples)*2, 2)]
conditions <- unlist(strsplit(conditions, ".", fixed = T))[seq(2, nrow(samples)*2, 2)]
samples$sampleType <- conditions
samples2 <- data.frame(samples$sampleType)
rownames(samples2) <- rownames(samples)
colnames(samples2) <- "sampleType"

############################################
# create phenoData annotated data frame
############################################
phenoData <- as(samples2, "AnnotatedDataFrame")

############################################
# create new MRexperiment object
############################################

print ("building new MR experiment")
dat <- newMRexperiment(input.data$counts,
                       phenoData=phenoData)

######################################################
# filter data based on k samples being over value a
######################################################

#print ("filtering data")
#tokeep <- which(rowSums(MRcounts(dat) >= opt$a) >= opt$k)
#dat <- dat[tokeep, ]

################################
# create normalisation factors
################################

print ("creating normalisation factors")
p = cumNormStat(dat)
dat.norm = cumNorm(dat, p = p)

################################
# output the normalised matrix
# outputs log2 counts
################################

print ("writing normalised count matrix")
mat.norm <- data.frame(MRcounts(dat.norm, norm = T, log = T))

print("removing NAs from norm matrix")
notnas <- which(!(is.na(colSums(mat.norm))))

mat.norm <- mat.norm[, notnas]
mat.norm$taxa <- rownames(mat.norm)

write.table(mat.norm, file = paste(opt$prefix, "norm.matrix", sep = "."), row.names = F, sep = "\t")

################################
# build model
################################

print ("building model and fitting")
dat.norm <- dat.norm[,notnas]

conds <- pData(dat.norm)$sampleType
mod <- model.matrix(~0+conds)
colnames(mod) <- levels(conds)

settings = zigControl(maxit = 1, verbose = F)
res <- fitZig(obj = dat.norm, mod = mod, control = settings)

zigfit <- res$fit
finalmod <- res$fit$design

print ("making contrasts")
# make contrasts
done <- c()
conts <- c()
for (i in colnames(mod)){
    for (j in colnames(mod)){
        cont <- paste(i, j, sep = "-")
        cont2 <- paste(j, i, sep = "-")
    	if (!(cont %in% done) & i != j){
            done <- append(cont, done)
            done <- append(cont2, done)
	    conts <- append(cont, conts)
         }
    }
}

contrast.matrix = makeContrasts(contrasts = conts, levels = finalmod)
fit2 <- contrasts.fit(zigfit, contrast.matrix)
fit2 <- eBayes(fit2)

######################################
# output the results of each contrast
######################################

print ("writing results")
result = data.frame()
for (cont in conts){
    c1 <- unlist(strsplit(cont, "-"))[1]
    c2 <- unlist(strsplit(cont, "-"))[2]
    res <- topTable(fit2, coef = cont, number = nrow(zigfit))
    res$group1 <- c1
    res$group2 <- c2
    res$taxa <- rownames(res)
    result <- rbind(result, res)
    }

# output the results
write.table(result, file = paste(opt$prefix, "diff.tsv", sep = "."), row.names = F, sep = "\t")

######################################
# produce various plots
######################################

###########
# mds plot
###########










