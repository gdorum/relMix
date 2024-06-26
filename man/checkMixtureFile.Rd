% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/checkMixtureFile.R
\name{checkMixtureFile}
\alias{checkMixtureFile}
\title{Load and check a mixture file}
\usage{
checkMixtureFile(filename)
}
\arguments{
\item{filename}{Path of the mixture file}
}
\value{
A list containing
\itemize{
\item \code{df} The loaded data frame, NULL if errors are present.
\item \code{warning} A list of strings describing the errors that ocurred but could be fixed or that do not prevent the execution of the program.
\item \code{error} A list of strings describing the errors that occurred that made it impossible to return a valid data frame.
If this list is not empty, then the data frame item will be null.}
}
\description{
Given a mixture file name, returns the loaded data frame along with any detected errors or warnings.
}
\details{
If warnings are found, the function attempts to fix them and explains what it has done in the warning messages.
If an error is found, checking stops and a NULL data frame is returned. The error is described in the error messages.
}
\examples{
mixfile <- system.file("extdata","mixture.txt",package="relMix")
result <- checkMixtureFile(mixfile);
print(result$df);
print(result$warning);
print(result$error);
}
\author{
Elias Hernandis
}
