###############################################
# R plugin for Plink - detection of epistasis #
# adjusted for covariates                     #
# NOTE: must use plink development version    #
###############################################

Rplink <- function(PHENO, GENO, CLUSTER, COVAR){
  epi.log <- function(snp){
    # covar positions are based on the input data table,
    # not the position they are passed in!
    # Assume SNP is the last column
    n.covar <- dim(COVAR)[2]
    test.snp <- COVAR[, n.covar]
    if(dim(COVAR)[2] > 2){
    	covars <- COVAR[, 1:(n.covar - 1)]
    }
    else {
        covars <- as.matrix(COVAR[, 2], ncol=1)
    }
    m <- glm(PHENO == 2 ~ covars + test.snp + snp + test.snp*snp, family=binomial(link="logit"))
    len <- dim(summary(m)$coefficients)[1]
    sum.snp <- summary(m)$coefficients[len,]
    ors <- exp(sum.snp[1])
    r <- c(dim(covars)[1], ors, sum.snp[2], sum.snp[3], sum.snp[4])
    names(r) <- c("NMISS", "OR", "SE", "STAT", "P")
    
    return(c(length(r), r))
  }
  apply(GENO, 2, epi.log)
}