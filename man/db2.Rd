\name{db2}
\alias{db2}
\docType{data}
\title{
Allele database for 22 markers
}
\description{
Frequencies for 22 loci from the prototype 24-plex STR panel from Thermo Fisher
}
\usage{data(db)}
\format{
  A data frame with 206 observations on the following 3 variables.
  \describe{
    \item{\code{Marker}}{a factor with levels corresponding to name of markers}
    \item{\code{Allele}}{a numeric vector denoting allele}
    \item{\code{Frequency}}{a numeric vector in (0,1)}
  }
}
\details{
The format is convenient for R
}
\source{
Hill et al. (2013) U.S. population data for 29 autosomal STR loci. Forensic Sci. Int. Genet. 7, e82-e83.
Hill et al. (2006) Allele Frequencies for 26 MiniSTR Loci with U.S. Caucasian, African American, and Hispanic Populations. http://www.cstl.nist.gov/biotech/strbase/NISTpop.htm
}

\examples{
data(db2)
}
\keyword{datasets}
