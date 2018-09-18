\name{mixLikDrop}
\alias{mixLikDrop}
\title{
Likelihoods for DNA mixtures with dropout and dropin
}
\description{
Function used by relMix. Computes the likelihood of a mixture conditioned on a given number of known and unknown contributors, and drop-in and dropout probabilities.
}
\usage{
mixLikDrop(R,G,D,di=0,alleleNames,afreq)
}
\arguments{
  \item{R}{
Vector of mixture alleles
}
\item{G}{
  List of genotypes. Each element is a vector with genotype for one individual 
}
\item{D}{
  List of dropout values (between 0 and 1) per contributor. Each element is a vector containing heterozygous and homozygous dropout probability for the given contributor
  }
  \item{di}{
  Drop-in value
  }
 \item{alleleNames}{
  Vector of allele names for the marker
}
\item{afreq}{
  Vector of allele frequencies for the marker
}
}

\value{
Returns a likelihood
}
\references{
The model is specified in the appendix of Haned et al. (2012) Exploratory data analysis for the interpretation of low template DNA mixtures, FSI Genetics Dec;6(6):762-74
}
\author{
Guro Dorum
}

\seealso{
\code{\link{relMix}}
}
\examples{
#Define alleles and frequencies
alleles <- 1:2
afreq <- c(0.5,0.5)
#Genotypes
gM <- c(1,1)
gC <- c(1,2)
#Mixture alleles
R <- c(1,2)
#Dropout and drop-in values
d <- 0.1
di <- 0.05
#No drop-in for first contributor
D <- list(c(0,0),c(d,d^2))
mixLikDrop(R=R,G=list(gM,gC),D=D,di=di,alleleNames=alleles,afreq=afreq)
}
