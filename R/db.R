#' Allele database
#'
#' Norwegian database with 17 EXS17 markers and 6 additional markers.
#' @usage data(db)
#' @format A data frame with 324 observations on the following 3 variables:
#'  \describe{
#' \item{\code{Marker}}{a factor with levels corresponding to name of markers}
#' \item{\code{Allel}}{a numeric vector denoting allele}
#' \item{\code{Frequency}}{a numeric vector in (0,1)}}
#' @source Dupuy et al. (2013), unpublished.
#' @examples
#' data(db)
#' #Checks that frequencies add to 1
#' lapply(split(db$Frequency,db$Marker),sum)
#' #Finds number of alleles for all markers
#' unlist(lapply(split(db$Frequency,db$Marker),length))
#' #A closer look at the marker SE33
#' SE33=db[db$Marker=="SE33",]
#' barplot(SE33$Frequency)
#' @keywords datasets
"db"