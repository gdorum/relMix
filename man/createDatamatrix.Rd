\name{createDatamatrix}
\alias{createDatamatrix}
\title{
Create data matrix with possible genotype combinations for specified individuals
}
\description{
A data matrix of genotypes for known individuals and all possible genotypes for unknown individuals is created. Can be used as input for relMix.
}
\usage{
createDatamatrix(locus, knownGenos, idsU)
}
\arguments{
  \item{locus}{
	A Familias locus containing information about the alleles
}
  \item{knownGenos}{
  List of known genotypes. Each element is a vector with genotype for one individual. The elements must be named
}
 \item{idsU}{
  Vector of indices for unknown individuals
}
}

\value{
Returns a data matrix of genotypes where each row corresponds to an individual.
}

\author{
Guro Dorum
}

\seealso{
See \code{\link{FamiliasLocus}} and \code{\link{relMix}}.
}
\examples{
require(Familias)
#Define alleles and frequencies
alleles <- 1:2
afreq <- c(0.5,0.5)
#Create locus object
locus <- FamiliasLocus(frequencies=afreq,name="M1",allelenames= alleles)
#Known genotypes of alleged father and mother, child's genotype is uknown
gAF <- c(1,1)
gMO <- c(1,1)
datamatrix <- createDatamatrix(locus,knownGenos=list(AF=gAF,MO=gMO),idsU=c("CH"))
}
