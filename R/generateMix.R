#Function that takes as input genotypes and creates a mixture
#with drop-in and dropout as specified by probabilities
generateMix <- function(G,alleles,afreq,D,di){
  
  #Probability for each allele of not being observed
  p <- .prAllMix(G,alleles,afreq,D,di)
  #Decide whether allele should be in mixture or not
  drop <- sapply(p,function(x) sample(0:1,size=1,prob=c(x,1-x)))
  #Return mixture
  alleles[which(drop==1)]
}