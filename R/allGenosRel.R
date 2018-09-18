#Finds all possible genotypes for specified individuals in a pedigree
#x: pedigree, a linkdat object
#partialmarker: a marker object compatible with x
#ids: index of pedigree members to find genotypes for
#Returns:
#datamatrix: a matrix with one row per individual. Number of columns
#             is 2 times the number of possible genotype combinations
allGenosRel <- function(x,partialmarker,ids){
  
  #Find all possible genotypes for one person
  allgenos <- allGenotypes(attr(partialmarker,"nalleles"))
  
  #grid.subset finds all possible genotype combinations for 
  #the individuals specified by ids, given as combination numbers
  grid.subset <- geno.grid.subset(x,partialmarker,ids)
  
  #Some of the genotypes may not fit with the specified relationship
  #Compute the likelihood of all genotypes to find which one does
  probs <- apply(grid.subset, 1, function(allg_rows) {
    partialmarker[ids, ] = allgenos[allg_rows, ]
    likelihood(x, locus1 = partialmarker)
  })
  
  #pick out the genotype combinations that has likelihood diff. from 0
  grid.subset <- grid.subset[probs!=0,,drop=FALSE]
  
  #Then, look up the genotype combinations in allgenos
  #to find the actual genotypes
  datamatrix <- numeric()
  for(i in 1:nrow(grid.subset)){
    datamatrix <- cbind(datamatrix,allgenos[grid.subset[i,],,drop=F])
  }
  
  rownames(datamatrix) <- ids
  as.data.frame(datamatrix)

}
