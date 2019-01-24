#' Relationship inference based on mixtures
#'
#' Calculates likelihoods for relationship inference involving mixtures and missing reference profiles,
#' including drop-in and dropout, mutations, silent alleles and theta correction.
#'
#' @param pedigrees A list of pedigrees defined using FamiliasPedigree in Familias
#' @param locus A Familias locus. Note that a silent allele must be indicated by 's' (and not 'silent' as in Familias)
#' @param R A vector of mixture alleles, or a list of such if there are multiple replicates
#' @param datamatrix Each line corresponds to one constellation of genotypes for the involved individuals.
#'  Indices of individuals must be given as rownames and must correspond to indices in the pedigree
#' @param ids  Index vector indicating which individuals are contributors to the mixture. The indices must correspond to
#'  indices in the pedigree
#' @param D List of dropout values (between 0 and 1) per contributor.
#'  Each element is a vector containing heterozygous and homozygous dropout probability for the given contributor
#' @param di Drop-in value (between 0 and 1)
#' @param kinship Defines the theta-parameter
#' @return The likelihoods for the pedigrees and detailed output for each term considered in the calculation.
#' @details The function requires the package \pkg{Familias} and calls on the function \code{FamiliasPedigree}.
#' @references Dorum et al. (2017) Pedigree-based relationship inference from complex DNA mixtures, Int J Legal Med., doi:10.1007/s00414-016-1526-x; \cr
#' Kaur et al. (2016) Relationship inference based on DNA mixtures, Int J Legal Med.;130(2):323-9; \cr
#' Egeland, Kling, Mostad (2015) \pkg{Familias}
#' @author Navreet Kaur, Thore Egeland, Guro Dorum
#' @seealso \code{\link{relMixGUI}} for a relMix GUI, and \code{\link{FamiliasLocus}} on how to create a Familias locus.
#' @examples #Example 1: paternity trio with mixture of mother and child
#' #Define alleles and frequencies
#' alleles <- 1:2
#' afreq <- c(0.4,0.6)
#' #Define pedigrees
#' persons <- c("CH","MO","AF")
#' ped1 <- Familias::FamiliasPedigree(id=persons, dadid=c("AF",NA, NA), momid=c("MO", NA,NA),
#'                         sex=c("male", "female", "male"))
#' ped2 <- Familias::FamiliasPedigree(id=c(persons, "TF"), dadid=c("TF", NA, NA,NA),
#'                         momid=c("MO", NA, NA,NA), sex=c("male", "female", "male", "male"))
#' pedigrees <- list(isFather = ped1, unrelated=ped2)
#' #Create locus object
#'locus <- Familias::FamiliasLocus(frequencies=afreq,name="M1",
#'                       allelenames= alleles)
#' #Known genotypes of alleged father and mother
#' gAF <- c(1,1)
#' gMO <- c(1,1)
#' #Mixture alleles
#' R <- c(1,2)
#' datamatrix <- createDatamatrix(locus,knownGenos=list(AF=gAF,MO=gMO),idsU=c("CH"))
#' #Define dropout and drop-in values
#' d <- 0.1
#' di <- 0.05
#' res <- relMix(pedigrees, locus, R, datamatrix, ids=c("MO","CH"),
#'              D=list(c(0,0),c(d,d^2)),di=di, kinship=0)
#' #LR=0.054
#' res$isFather/res$unrelated
#'
#' #Example 2: Exhaustive example with silent allele, mutations, dropout and drop-in
#' #H1: Contributors are mother and child
#' #H2: Contributors are mother and unrelated
#' #Possible dropout in both contributors
#' gMO <- c(1,1) #Mother's genotype
#' R <- 1 #Mixture alleles
#' #Mother/child pedigree
#' persons <- c("CH","MO")
#' ped1 <- Familias::FamiliasPedigree(id=persons, dadid=c(NA,NA), momid=c("MO", NA),
#'                         sex=c("male", "female"))
#' ped2 <- Familias::FamiliasPedigree(id=c(persons), dadid=c(NA, NA),
#'                         momid=c( NA, NA),
#'                         sex=c("male", "female"))
#' pedigrees <- list(H1 = ped1, H2=ped2)
#' #Alleles and frequencies:
#' #When silent alleles are involved, a custom mutation matrix is required.
#' #No mutations are possible to or from silent alleles.
#' #We create the mutation model with FamiliasLocus and modify it before it is
#' #passed on to relMix
#' alleles <- c(1,2,'silent')
#' afreq <- c(0.4,0.5,0.1)
#' #Create initial locus object with mutation matrix
#' locus <- Familias::FamiliasLocus(frequencies=afreq,name='M1',
#'                       allelenames= alleles, MutationModel='Equal',
#'                       femaleMutationRate=0.1,maleMutationRate=0.1)
#' #Modify mutation matrix from Familias:
#' #Silent allele must be given as 's' (not 'silent' as in Familias)
#' newAlleles <- c(alleles[-length(alleles)],'s')
#' mm <- locus$femaleMutationMatrix
#' colnames(mm) <- rownames(mm) <- newAlleles
#' #Create new locus object with modified mutation matrix
#' locus <- Familias::FamiliasLocus(frequencies=afreq,name='M1',
#'                       allelenames= newAlleles, MutationModel='Custom', MutationMatrix=mm)
#' knownGenos <- list(gMO)
#' names(knownGenos) <- c("MO")
#' datamatrix <- createDatamatrix(locus,knownGenos,ids="CH")
#' d <- 0.1 #Dropout probability for both contributors
#' di <- 0.05
#' res2 <- relMix(pedigrees, locus, R, datamatrix, ids=c("MO","CH"),
#'               D=list(c(d,d^2),c(d,d^2)),di=di, kinship=0)
#' #LR=1.68
#' res2$H1/res2$H2
#'
#' @export
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

  alt <- Familias::FamiliasPosterior(pedigrees, myloci, datamatrix, ref = ref, kinship = kinship)
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