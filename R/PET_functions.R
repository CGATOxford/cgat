############################################
# Pleiotropy Estimation and Testing (PET)  #
# Functions by Zhang et al Genet Epidemiol #
############################################

estimateTraitResidual <- function(data, trait, 
                                  covars, link, distribution) {
  # calculate the covariate adjusted residuals
  # from a linear model for a trait
  
  # check if the trait is 0/1 or 1/2 for binary traits
  if((max((data[[trait]])[!is.na(data[[trait]])]) == 2) &&
       (length(unique((data[[trait]])[!is.na(data[[trait]])])) == 2)){
    trait <- paste0(trait, " == 2")
  }
  
  t.covars <- paste(covars, collapse=" + ", sep=" + ")
  t.form <- paste(trait, t.covars, sep=" ~ ")
  t.mod <- glm(t.form, family=distribution(link), data=data)
  
  t.residuals <- t.mod$residuals
  
  return(t.residuals)  
}

constructDataMatrix <- function(trait1, trait2, geno) {
  # create the matrix for the LMM with X, T and IT components
  
  n = length(geno)
  design.vec <- c(rep(1, n), rep(2, n))
  Y <- c(trait1, trait2)
  X <- c(geno, geno)
  
  I.vec <- design.vec * X
  mu <- mean(trait1 + trait2)
  
  input.matrix <- data.frame(cbind(1, Y, X, design.vec, I.vec))
  
  colnames(input.matrix) <- c("mu", "Y", "X", "Tvec", "Ivec")
  
  return(input.matrix)
}


calculateResidualCovariance <- function(data, method="ML", model="LM"){
  # calculate residual covariance from LMM residuals
  if(model == "LMM"){
    require(nlme)
    lmm.fit <- lme(Y ~ X + Tvec + Ivec, random=~1|mu, data=data,
                   na.action=na.omit,
                   method=method)
  }
  else if(model == "LM"){
    lmm.fit <- glm(Y ~ X + Tvec + Ivec + mu, data=data, 
                   na.action=na.omit)
  }
  
  # the returned residuals is a 2-segment vector
  covar.mat <- cov(matrix(resid(lmm.fit), ncol=2))
  return(covar.mat[1, 2])
}


calcPleiotropyCorrelation <- function(trait1, trait2, residual.covar){
  # calculate the pleiotropy correlation coefficient (PCC) from
  # the difference between the observed trait covariance and
  # residual model covariance, standardized to the
  # product of the trait standard deviations
  
  trait.covar <- cov(trait1, trait2)
  trait1.sd <- sd(trait1)
  trait2.sd <- sd(trait2)
  
  delta <- trait.covar - residual.covar
  rho <- (abs(delta)/(trait1.sd * trait2.sd))

  return(rho)  
}

PleiotropyEstimationTest <- function(data.set, trait1, trait2,
                                     genotypes){
  
  # calculate the pleiotropy correlation coefficient
  t1 <- data.set[[trait1]]
  t2 <- data.set[[trait2]]
  geno <- data.set[[genotypes]]
  
  input.matrix <- constructDataMatrix(t1, t2, geno)
  resid.covar <- calculateResidualCovariance(input.matrix, "ML", "LM")
  
  pcc <- calcPleiotropyCorrelation(t1, t2, resid.covar)
  
  return(pcc)
}

BootstrapPET <- function(data, indices, trait1, trait2, genotypes){
  data.mat <- data[indices, ]
  pcc = PleiotropyEstimationTest(data.mat, trait1, trait2, genotypes)
  pcc
}


PETB <- function(data, trait1, trait2, genotypes, resamples,
                 plot=FALSE){
  # calculate the assymetric 2-sided P-value from bootstrap resampling
  require(boot)
  
  # boot is bootstrapping around the observed value
  # however, it should be bootrapping around 0 to test
  # PCC > 0
  
  boot.out <- boot(data=data, statistic=BootstrapPET, R=resamples,
                   trait1=trait1, trait2=trait2, genotypes=genotypes)
  boot.vals <- boot.out$t - mean(boot.out$t)
  
  if(plot){
    par(mfrow=c(1, 2))
    hist(boot.vals, breaks=100, xlab="PCC Null")
    abline(v=boot.out$t0, lty="dashed", main="Histogram of null PCC")
    
    qqnorm(boot.vals, xlab="Quantiles of Standard Normal PCC Null")
    qqline(boot.vals, distribution = qnorm)
  }
  
  # calculate the left and right tail p-vals
  # 2-sided p is twice the smallest one-tailed p-value
  # assumes symmetry in distribution of bootstrapped
  # pcc values
  
  pval <- sum(abs(boot.out$t0) > boot.vals)/(boot.out$R)
  
  return(c(boot.out$t0, pval))
}

# a function that loops over a SNP list and calculates
# the PCC and p-value for each
loopPET <- function(data.df, trait1, trait2, trait1.link, trait2.link,
                    trait1.mod, trait2.mod, covars, resamples, snp.list){
  pet_results <- list()
  for(s in 1:length(snp.list)){
    subsample <- data.df[(!is.na(data.df[[trait1]])) && (!is.na(data.df[[trait2]]))
                         && (!is.na(data.df[[snp.list[s]]])), ]    

    trait1.vec <- estimateTraitResidual(data=subsample, trait=trait1,
                                        covars=covars, link=trait1.link,
                                        distribution=trait1.mod)
    
    trait2.vec <- estimateTraitResidual(data=subsample, trait=trait2,
                                        covars=covars, link=trait2.link,
                                        distribution=trait2.mod)
    
    geno <- subsample[[snp.list[s]]]
    
    data.mat <- data.frame(cbind(trait1.vec, trait2.vec, geno))
    colnames(data.mat) <- c("Trait1", "Trait2", "Geno")
    
    out <- PETB(data = data.mat, trait1 = "Trait1", trait2 = "Trait2",
                genotypes = "Geno", resamples=resamples, plot=F)
    
    pet_results[[snp.list[s]]] <- out
  }
  return(pet_results)
}



# this function is from stackexchange:
# http://stats.stackexchange.com/questions/15011/generate-a-random-variable-with-a-defined-correlation-to-an-existing-variable
# it generates a correlated variable given a value for rho
getBiCop <- function(n, rho, mar.fun=rnorm, x = NULL, ...) {
  if (!is.null(x)) {X1 <- x} else {X1 <- mar.fun(n, ...)}
  if (!is.null(x) & length(x) != n) warning("Variable x does not have the same length as n!")
  
  C <- matrix(rho, nrow = 2, ncol = 2)
  diag(C) <- 1
  
  C <- chol(C)
  
  X2 <- mar.fun(n)
  X <- cbind(X1,X2)
  
  # induce correlation (does not change X1)
  df <- X %*% C
  
  ## if desired: check results
  #all.equal(X1,X[,1])
  #cor(X)
  
  return(df)
}