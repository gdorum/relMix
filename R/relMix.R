#Computes the likelihood for a mixture with relatives
#Drop-in and dropout can be included
#Can be several replicates - R is a list with one replicate per element (alternatively, if R is a vector will be converted to list)
#kinship calculations are done with relMix
#ids gives the index of contributors to mixture, corresponding with indiced in pedigree
relMix <- function(pedigrees, locus, R, datamatrix, ids, D=rep(list(c(0,0)),length(ids)),di=0, kinship=0){
  
  if(class(pedigrees)=='list'){
      for (i in 1:length(pedigrees)) {
        if (class(pedigrees[[i]]) != "FamiliasPedigree") 
          stop("Argument pedigrees should be a list of Familias pedigrees")
        
        if(any(sapply(pedigrees,function(x) any(!ids%in%x$id))))
            stop("Argument ids does not correspond with indices in pedigree")
    } 
  }else {
    if(class(pedigrees) != "FamiliasPedigree") 
      stop("Argument pedigrees should be a list of Familias pedigrees")
    
    if(any(!ids%in%pedigrees$id))
      stop("Argument ids does not correspond with indices in pedigree")
  }
  if(class(R) != "list") R <- list(R)
  if (class(locus) != "FamiliasLocus") 
    stop("Argument locus should be a Familias locus")
  if (!is.data.frame(datamatrix)) 
    stop("Argument datamatrix should be a dataframe")
  
  #Which pedigree (if several) to use as reference
  ref <- ifelse(length(pedigrees)==2, 2,1)
  #Compute probability of kinship
  #Number of possible genotype combinations
  n <- dim(datamatrix)[2]/2
  #Make one locus combination
  myloci <- vector("list", n)
  for (i in 1:n) {
    myloci[[i]] = locus
    myloci[[i]]$locusname = paste("Term", i, sep = "")
  }
  
  alt <- FamiliasPosterior(pedigrees, myloci, datamatrix, ref = ref, kinship = kinship)
  result <- apply(alt$likelihoodsPerSystem, 2, sum)
  terms <- alt$likelihoodsPerSystem
  
  #terms that gives likelihood 0 under both hypotheses we can remove
  idx <- which(apply(terms,1,function(x) any(x>0)))
  terms <- terms[idx,,drop=F]
  
  afreq <- locus$alleles
  alleles <- names(afreq)
 
  ix <- seq(1,ncol(datamatrix),2)
  ix <- ix[idx]
  #If ix is empty, all terms have likelihood 0 under both hypotheses and no need to compute mixture probabilities
  l <- numeric()
  if(length(ix)>0){
    M <- lapply(ix,function(i){datamatrix[rownames(datamatrix)%in%ids,i:(i+1)]})
    #Compute probability of mixture given contributors:
    #If there are silent alleles, these won't be sent to mixLikDrop
    #as they are indifferent to dropout and drop-in
    alleleNames <- alleles[alleles!='s']
    freqs <- afreq[alleles!='s']
    for(i in 1:length(M)){
      
      G <- sapply(ids,function(j) list(M[[i]][j,]))
      #Set NA for silent alleles
      G <- lapply(G,function(x) {x[x=='s'] <- NA;x})
      #Only one replicate
      #l[i] <- mixLikDrop(R=R,G=G,D=D,di=di,alleleNames=alleleNames,afreq=freqs)
      #Several replicates - the likelihood is just the product of the likelihood for each replicate
      l[i] <- prod(sapply(X=R,FUN=mixLikDrop,G=G,D=D,di=di,alleleNames=alleleNames,afreq=freqs))
    }
  }
 
  
  #Combine kinship and mixture calculations
  #Likelihood for each term
  termlik <- terms*l
  #terms that gives likelihood 0 under both hypotheses we can remove
  idx <- which(apply(termlik,1,function(x) any(x>0)))
  termlik <- termlik[idx,,drop=F]
  #Likelihood for each hypothesis
  res <- lapply(1:ncol(termlik),function(i) sum(termlik[,i]))
  L <- c(res, list(termlik))
  names(L) <- c(colnames(terms),"terms")
  L
}