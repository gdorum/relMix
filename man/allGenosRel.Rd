\name{allGenosRel}
\alias{allGenosRel}

\title{
Finds all possible genotypes for specified individuals in a pedigree
}
\description{
Finds all possible genotypes for specified individuals in a pedigree. Pedigree and marker must \code{\link{linkdat}} objects. 
}
\usage{
allGenosRel(x, partialmarker, ids) 
}

\arguments{
  \item{x}{
  Pedigree, a linkdat object
}
\item{ids}{
  Vector of indices of pedigree members to find genotypes for
}
\item{partialmarker}{
  A marker object compatible with x
}
}
\value{
Return a data matrix with one row per individual. The number of columns
is 2 times the number of possible genotype combinations
}

\author{
Guro Dorum
}

\examples{
#Find all possible genotypes for a child conditioned on the mother's genotype
library(paramlink)
alleles <- 1:2
afreq <- c(0.5,0.5)
#Mother's genotype
gM <- c(1,1)
ped <- nuclearPed(1)
m <-  marker(ped,2,gM,alleles=alleles,afreq=afreq)
allGenosRel(ped, m, ids=3) 
}
