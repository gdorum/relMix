#' Create data matrix with possible genotype combinations for specified individuals
#'
#' A data matrix of genotypes for known individuals and all possible genotypes for unknown individuals is created.
#' @param locus A list of class \code{\link[Familias]{FamiliasLocus}} containing information about the locus
#' @param knownGenos List of known genotypes. Each element is a vector with genotype for one individual. The elements must be named
#' @param idsU Vector of indices for unknown individuals
#' @return A data matrix of genotypes where each row corresponds to an individual.
#' @author Guro Dorum
#' @seealso \code{\link[Familias]{FamiliasLocus}}.
#' @examples
#' #Define alleles and frequencies
#' alleles <- 1:2
#' afreq <- c(0.5,0.5)
#' #Create locus object
#' locus <- Familias::FamiliasLocus(frequencies=afreq,name="M1",allelenames= alleles)
#' #Known genotypes of alleged father and mother, child's genotype is uknown
#' gAF <- c(1,1)
#' gMO <- c(1,1)
#' datamatrix <- createDatamatrix(locus,knownGenos=list(AF=gAF,MO=gMO),idsU=c("CH"))
#' @export
createDatamatrix <- function(locus,knownGenos,idsU=NULL){

  nU <- length(idsU)
  origAlleles <- names(locus$alleles)
  alleles <- 1:length(origAlleles)

  #If there are unknown contributors in the mixture
  if(nU > 0){

    #Find all possible genotypes for one person
    allgenos <- allGenos(alleles)
    #Find combinations for all uknown contributors
    grid.subset <- expand.grid(rep(list(1:nrow(allgenos)),nU))

    #Then, look up the genotype combinations in allgenos
    #to find the actual genotype
    datamatrixU <- NULL
    for(i in 1:nrow(grid.subset)){
      datamatrixU <- cbind(datamatrixU,allgenos[as.numeric(grid.subset[i,]),,drop=F])
    }
    #Convert back to original allele names
    datamatrixU <- t(apply(datamatrixU,1,function(y) origAlleles[y]))
    known <- t(sapply(knownGenos,function(x) rep(x,dim(datamatrixU)[2]/2)))
    datamatrix <- rbind(known,datamatrixU)

    } else{ #Only known contributors
      datamatrix <- matrix(unlist(knownGenos),nrow=length(knownGenos),byrow=TRUE)
  }


  #Check if there is a silent allele; if so, must account for a possible
  #silent allele in all known homozygotes
  if(any(origAlleles=="s")){
    nK <- length(knownGenos)
    knownGenosS <- lapply(knownGenos,function(x){
      u <- unique(unlist(x))
      if(length(u)<2) return(c(u,'s'))
      else return(u)
    })
    #Genotype combinations
    genos <- expand.grid(rep(list(1:2),nK))
    #All genotypes for known contributors
    g <- list(matrix(unlist(knownGenos),nrow=nK,byrow=TRUE),matrix(unlist(knownGenosS),nrow=nK,byrow=TRUE))
    #Look up the genotype combinations in g to find the actual genotype
    gtall <- numeric()
    for(i in 1:nrow(genos)){
      gt <- numeric()
      for(j in 1:nK){
        gt <- rbind(gt,g[[as.numeric(genos[i,j])]][j,,drop=F])
      }
      gtall <- cbind(gtall,gt)
    }
    gtall <- t(apply(gtall, 1, rep, dim(datamatrix)[2]/2))
    if(nU > 0){
      A <- matrix(datamatrixU,ncol=2,byrow=TRUE)
      datamatrix <- rbind(gtall,c(t(apply(A,2,rep,each=nrow(genos)))))
    } else datamatrix <- gtall
  }


  datamatrix <- as.data.frame(datamatrix)
  rownames(datamatrix) <- c(names(knownGenos),idsU)
  datamatrix
}

