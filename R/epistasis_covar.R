###############################################
# R plugin for Plink - detection of epistasis #
# adjusted for covariates                     #
# NOTE: must use plink development version    #
###############################################

Rplink <- function(PHENO, GENO, CLUSTER, COVAR){
  epi.log <- function(snp){
    # covar positions are based on the input data table,
    # not the position they are passed in!
    # Assume SNP is the first column
    test.snp = COVAR[, 1]
    covars = COVAR[,(2:dim(COVAR)[2])]
    m <- glm(PHENO == 2 ~ snp + covars + test.snp*snp, family=binomial(link="logit"))
    len <- dim(summary(m)$coefficients)[1]
    sum.snp <- summary(m)$coefficients[len,]
    ors <- exp(sum.snp[1])
    r <- c(ors, sum.snp[2], sum.snp[3], sum.snp[4])
    names(r) <- c("OR", "SE", "STAT", "P")
    
    return(c(length(r), r))
  }
  apply(GENO, 2, epi.log)
}