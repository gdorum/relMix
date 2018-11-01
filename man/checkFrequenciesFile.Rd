\name{checkFrequenciesFile}
\alias{checkFrequenciesFile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Load a frequency database file and check it against a mixture file
}
\description{
Given a frequency database file name, returns the loaded data frame along with any detected errors or warnings.
}
\usage{
checkFrequenciesFile(filename)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
    path to the frequency database file
}
  \item{mix}{
    a dataframe with the mixture data
  }
}
\details{
The mixture data is used to perform more advanced checks, such as to make sure all alleles present in the mixture file have an entry in the frequency database.

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
\code{\link{checkMixtureFile}} for information on how to load a mixture file
}
\examples{
mix <- checkMixtureFile('data/valid/mix.txt');

# note: the mixture dataframe is passed as an argument
# if the previous check failed, the program should not continue
# with the frequencies file check
freqs <- checkFrequenciesFile('data/valid/freq.txt', mix$df);
}