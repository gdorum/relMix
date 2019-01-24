#' Load and check a mixture file
#'
#' Given a mixture file name, returns the loaded data frame along with any detected errors or warnings.
#' @param filename Path of the mixture file
#' @return A list containing
#' \itemize{
#' \item{df} {The loaded data frame, NULL if errors are present.}
#' \item{warning} {A list of strings describing the errors that ocurred but could be fixed or that do not prevent the execution of the program.}
#' \item{error} {A list of strings describing the errors that ocurred that made it imposible to return a valid data frame.
#' If this list is not empty, then the dataframe item will be null.}}
#' @details If warnings are found, the function attempts to fix them and explains what it has done in the warning messages.
#' If an error is found, checking stops and a NULL dataframe is returned. The error is described in the error messages.
#' @examples
#' \dontrun{
#' mixfile <- system.file("extdata","mixture.txt",package="relMix")
#' result <- checkMixtureFile(mixfile);
#' print(result$df);
#' print(result$warning);
#' print(result$error);
#' }
#' @author Elias Hernandis
#' @export
checkMixtureFile <- function(filename) {
    r <- commonChecks(filename, "reference file");
    df <- r$df;
    warning <- r$warning;
    error <- r$error;

    # check header column count
    if (!is.null(df) && ncol(df) < 3) {
        error <- "Incorrect number of columns"
    }

    referenceHeader <- c("SampleName", "Marker")

    if (length(error) == 0 && any(referenceHeader != names(df)[1:2])) {
        warning <- append(warning, "Column titles must be written without spaces and with the first letter uppercased (fixed).");

        fixedHeader <- sapply(names(df)[1:2], titleize);
        fixedHeader[1] <- "SampleName"; # special case for the first column header, which is "SampleName" instead of "Sample Name"

        if (all(fixedHeader != names(df)[1:2])) {
            error <- append(error, paste("There are problems with the header row of the mixture table. Please make sure the first two colums are named as follows: ", paste(referenceHeader,collapse=", ")));
        } else {
            # if the errors are minor, fix the header for the user
            names(df) <- referenceHeader;
        }
    }

    # If there are two or more sample names, i.e. several replicates,
    #check that they have equal number of rows. Otherwise return an error
    sampleNames <- unique(df$SampleName);
    if (length(error) == 0 && length(sampleNames) > 1) {
      if(abs(max(table(df$SampleName)) - min(table(df$SampleName)))>0){
          error <- paste("Samples", paste(sampleNames,collapse=" and "), "have unequal number of rows in mixture file.");
      }
    }

    # Check all allele reference data is numeric
    numCheck <- apply(df[,3:ncol(df),drop=FALSE],2,function(x) all(sapply(x, is.numeric)))
    if (length(error) == 0 && !all(numCheck)) {
           error <- c("There are values that are not numeric in the following columns of the mixture file",
                      paste(names(df)[3:ncol(df)][!numCheck],collapse=", "));
    }
    else if (!all(numCheck)) {
      error <- append(error,c("There are values that are not numeric in the following columns of the mixture file",
                      paste(names(df)[3:ncol(df)][!numCheck],collapse=", ")))
    }

    if (length(error) > 0) {
        return(list(df=NULL, warning=NULL, error=error));
    }

    return(list(df=df, warning=warning, error=NULL));
}
