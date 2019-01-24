#' Create a mixture of genotypes with simulated drop-in and dropout
#'
#' Takes as input genotypes and creates a mixture. Alleles drop in and out of the mixture with the specified probabilities
#' @param G List of genotypes. Each element is a vector with genotype for one individual
#' @param alleles Numeric or character Vector of allele names for the marker
#' @param afreq Numeric vector of allele frequencies for the marker
#' @param D List of dropout values (between 0 and 1) per contributor.
#' Each element is a vector containing heterozygous and homozygous dropout probability for the given contributor
#' @param di Drop-in value (between 0 and 1)
#' @return A vector of mixture alleles.
#' @author Guro Dorum
#' @examples
#' #Define alleles and frequencies
#' alleles <- 1:2
#' afreq <- c(0.5,0.5)
#' #Genotypes
#' gM <- c(1,1)
#' gC <- c(1,2)
#' #Dropout and drop-in values
#' d <- 0.1
#' di <- 0.05
#' #No drop-in for first contributor
#' D <- list(c(0,0),c(d,d^2))
#' R <- generateMix(G=list(gM,gC),alleles,afreq,D=D,di=di)
#' @export
generateMix <- function(G,alleles,afreq,D,di){

  #Probability for each allele of not being observed
  p <- .prAllMix(G,alleles,afreq,D,di)
  #Decide whether allele should be in mixture or not
  drop <- sapply(p,function(x) sample(0:1,size=1,prob=c(x,1-x)))
  #Return mixture
  alleles[which(drop==1)]
}