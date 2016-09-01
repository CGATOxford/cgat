geneListSnpColocQtl <- function(gene_list, results_table, MAF_table, eqtl_table,
                                trait_type, prev=NULL){  
  require(coloc)
  results_list = list()
  for(i in 1:length(gene_list)){
    # need to remove duplicated SNPs
    gene <- na.omit(eqtl_table[eqtl_table$Gene == gene_list[i],])
    gene <- gene[!duplicated(gene$SNP),]
    
    share.snps <- intersect(results_table$SNP, gene$SNP)
    
    # make sure they actually share some SNPs
    if(dim(gene)[1] > 0 & length(share.snps) > 6){
      res_match <- results_table[results_table$SNP %in%share.snps, ]
      gene_match <- gene[gene$SNP %in% res_match$SNP, ]
      all_snps <- intersect(res_match$SNP, gene_match$SNP)
      
      mafs <- MAF_table[MAF_table$SNP %in% all_snps, ]$MAF
      gene_match <- gene_match[order(gene_match$SNP, decreasing=T), ]
      res_match <- res_match[order(res_match$SNP, decreasing=T), ]

      # this seems to have a hissy fit over the SNPs that
      # match, even if > 1 SNP are in the two data sets.
      # Does this require a minimum number of SNPs to match??
      gene.res <- tryCatch({coloc.abf(dataset1=list(pvalues=res_match$P,
                                                    N=max(res_match$NMISS),
                                                    type=trait_type, s=prev,
                                                    snp=all_snps),
                                      dataset2=list(pvalues=gene_match$P,
                                                    N=max(gene_match$NMISS),
                                                    type="quant", snp=all_snps),
                                      MAF=mafs)$summary
                            },
                           warning = function(warn){
                             print(paste("MY_WARNNG: ", warn))
                           },
                           error = function(err){
                             print(paste("MY_ERROR: ", err))
                             gene.res <- c(length(all_snps), 0, 0, 0, 0, 0)
                             return(gene.res)
                           },
                           finally = {
                             print("Insufficient matching SNPs")
                           })
      results_list[[gene_list[i]]] <- gene.res
    }
    else {
        gene.res <- c(0, 0, 0, 0, 0, 0)
	results_list[[gene_list[i]]] <- gene.res
    } 
  }
  result.df <- data.frame(do.call(rbind, results_list))
  return(result.df)
}


TwoTraitSnpColocQtl <- function(trait1_table, trait2_table, MAF_table,
                                trait1_type, trait2_type, prev1=NULL,
                                prev2=NULL){
  require(coloc)
  results_list = list()
  # need to remove duplicated SNPs
  trait1_table <- trait1_table[!duplicated(trait1_table$SNP),]
  trait2_table <- trait2_table[!duplicated(trait2_table$SNP),]
    
  share.snps <- intersect(trait1_table$SNP, trait2_table$SNP)
  all_snps <- intersect(share.snps, MAF_table$SNP)  
    
  # make sure they actually share some SNPs
  if(length(share.snps) > 0){
    trait1_match <- trait1_table[trait1_table$SNP %in%all_snps, ]
    trait2_match <- trait2_table[trait2_table$SNP %in% trait1_match$SNP, ]
    mafs <- MAF_table[MAF_table$SNP %in% trait1_match$SNP, ]$MAF

    res <- coloc.abf(dataset1=list(pvalues=trait1_match$P,
                                        N=max(trait1_match$NMISS),
                                        type=trait1_type, s=prev1, snp=all_snps),
                          dataset2=list(pvalues=trait2_match$P,
                                        N=max(trait2_match$NMISS),
                                        type=trait2_type, s=prev2, snp=all_snps),
                          MAF=mafs)
    results_list[[1]] <-c(res$summary)
  }
  result.df <- do.call(rbind, results_list)
  return(result.df)
}

