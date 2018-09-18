\name{generateMix}
\alias{generateMix}
\title{
Create a mixture of genotypes with simulated drop-in and dropout
}
\description{
Takes as input genotypes and creates a mixture. Alleles drop in and out of the mixture with the specified probabilities
}
\usage{
generateMix(G,alleles,afreq,D,di)
}
\arguments{
  \item{G}{
  List of genotypes. Each element is a vector with genotype for one individual 
}
 \item{alleles}{
  Vector of allele names for the marker
}
\item{afreq}{
  Vector of allele frequencies for the marker
}
 \item{D}{
  List of dropout values (between 0 and 1) per contributor. Each element is a vector containing heterozygous and homozygous dropout probability for the given contributor
}
 \item{di}{
  Drop-in value
  }
}

\value{
Returns a vector of mixture alleles.
}

\author{
Guro Dorum
}

\examples{
#Define alleles and frequencies
alleles <- 1:2
afreq <- c(0.5,0.5)
#Genotypes
gM <- c(1,1)
gC <- c(1,2)
#Dropout and drop-in values
d <- 0.1
di <- 0.05
#No drop-in for first contributor
D <- list(c(0,0),c(d,d^2))
R <- generateMix(G=list(gM,gC),alleles,afreq,D=D,di=di)
}
