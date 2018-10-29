\name{checkMixtureFile}
\alias{checkMixtureFile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Load and check a mixture file
}
\description{
Given a mixture file name, returns the loaded data frame along with any detected errors or warnings.
}
\usage{
checkMixtureFile(filename)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
    path to the mixture file
}
}
\details{
If warnings are found, the function attempts to fix them and explains what it has done in the warning messages. If an error is found, checking stops and a NULL dataframe is returned. The error is described in the error messages.
}
\value{
A list containing the data frame and two lists of messages:
\item{df }{The loaded data frame, NULL if errors are present}
\item{warning }{A list of strings describing the errors that ocurred but could be fixed or that do not prevent the execution of the program.}
\item{error }{A list of strings describing the errors that ocurred that made it imposible to return a valid data frame. If this list is not empty, then the dataframe item will be null.}
}

\author{
Elias Hernandis
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
result <- checkMixtureFile('data/valid/mix.txt');
print(result$df);
print(result$warning);
print(result$error);
}