\name{allGenos}
\alias{allGenos}
\title{
Find all possible genotypes
}
\description{
Finds all possible genotypes based on input alleles.
}
\usage{
allGenos(alleles)
}
\arguments{
  \item{alleles}{
  Vector of allele names
}
}
\value{
Returns matrix of all possible genotypes, one row per genotype.
}
\author{
Guro Dorum
}
\examples{
alleles <- 1:3
allGenos(alleles)
}