#Likelihoods for mixtures with dropout and dropin

#Wrapper for prAllMix() which calculates probability of seeing alleles
#in mixture based on genotypes of known contributors

#Model specified in Appendix of Haned et al. (2012)
#Given a mixture with number of contributors n known, and genotypes of contributors
#R: mixture alleles
#G: list of genotype(s) of contributor(s)
#n: number of contributors
#D: list of dropout-probabilities. Two probabilities per contributor, heterozygous and homozygous
#di: dropin probability in the interval [0,1] 
#alleles: all alleles for given locus
#afreq: frequencies
mixLikDrop <- function(R,G,D,di=0,alleleNames,afreq){
  
  if(length(D)!=length(G)) stop("You must specify dropout probabilities for all contributors")
  #if(length(G)!=n) stop("You must specify one genotype for each contributor")
  
#   #If allele names does not start from 1, rename
#   alleles <- 1:length(afreq)
#   R <- alleles[alleleNames%in%R] 
#   G <- lapply(G,function(x) { a <- alleles[alleleNames%in%x]; if(length(a)<2) a <- c(a,a);a}) 
 
  #Probability of not seeing allele
  p <- .prAllMix(G,alleleNames,afreq,D,di)
  #If allele is in mixture, 1-prob(not seeing allele)
  probs <- ifelse(alleleNames %in% R,1-p,p)
  
  #Return likelihood
  prod(probs)

}

#calculates probability of not observing alleles
#based on genotypes of known contributors
#and given dropout and drop-in
.prAllMix <- function(G,alleles,afreq,D,di){
  
  n <- length(G)
  na <- length(alleles)
  #Go through all alleses and find the probability of not observing them in the mixture
  p <- numeric(na)
  for(i in 1:na){
    #na.rm=TRUE makes make it possible for individuals to have 1 or 2 NA (silent) alleles
    nn <- sapply(G,function(x) sum(x==alleles[i],na.rm=TRUE))
    #Probability of dropout and no dropin, i.e. the allele is not in the mixture
    p[i] <- (1-di*afreq[alleles[i]]) * prod(unlist(sapply(1:n,function(j) D[[j]][nn[j]])))
  }
  p
}
