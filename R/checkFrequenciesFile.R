#' Load and check a frequency file
#'
#' Loads a frequency database file and compares it against mixture data to check for common errors.
#'
#' @param filename Path of the frecuency database file
#' @param mix Data frame with mixture data. See relMix vignette for description of the format
#' @return A list containing
#' \itemize{
#' \item {\code{df}} {Data frame with frequencies}
#' \item {\code{warning}} {List of strings describing the errors that ocurred but could be fixed or that do not prevent
#' the execution of the program.}
#' \item {\code{error}} {List of strings describing the errors that ocurred that made it imposible to return a valid data frame.
#' If this list is not empty, then the dataframe item will be NULL}}
#' @details
#' The mixture data is used to perform more advanced checks, such as to make sure all alleles present
#' in the mixture file have an entry in the frequency database.
#' If warnings are found, the function attempts to fix them and explains what it has done in the warning messages.
#' If an error is found, checking stops and a NULL dataframe is returned. The error is described in the error messages.
#' @seealso \code{\link{checkMixtureFile}} for information on how to load a mixture file.
#' @examples
#' \dontrun{
#' mixfile <- system.file("extdata","mixture.txt",package="relMix")
#' mix <- checkMixtureFile(mixfile)
#' # note: the mixture dataframe is passed as an argument
#' # if the previous check failed, the program should not continue
#' # with the frequencies file check
#' freqfile <- system.file('extdata','frequencies22Markers.txt',package='relMix')
#' freqs <- checkFrequenciesFile(freqfile, mix$df)
#' }
#' @author Elias Hernandis
#' @importFrom utils combn
#' @export
checkFrequenciesFile <- function(filename, mix) {
    r <- commonChecks(filename, "frequencies file");
    df <- r$df;
    warning <- r$warning;
    error <- r$error;

    if (is.null(df)) {
        return(list(df=NULL, error=error, warning=NULL));
    }

    # Make sure it says "Allele" on cell (1,1)
    if (length(error) == 0 && names(df)[1] != "Allele") {
        if (titleize(names(df)[1]) == "Allele") {
            warning <- append(warning, "The first column must be named \"Allele\" (fixed)");
            names(df)[1] <- "Allele";
        }
    }

    # Check that all alleles are numeric
    if (length(error) == 0 && !all(sapply(df$Allele, is.numeric))) {
        error <- append(error, "There are values that are not numeric in the Allele column of the frequency file");
    }

    # Check that all frequencies are either numeric or NA
    if (length(error) == 0) {
        for (i in 2:ncol(df)) {
            if (!all(sapply(df[,i], is.numeric))) {
                error <- append(error, paste("There are non-numeric frequencies in column", i, "of the frequency file."));
            }
        }
    }

    #If more than one marker
    if(ncol(df)>2){
    # Check for duplicate markers
    # TODO: this error will never trigger b/c R renames duplicate columns when
    # loading w/ headers
    if (length(error) == 0 && unique(names(df)[-1]) != names(df)[-1]) {
        error <- append(error, paste("There are duplicate markers in your frequency table."));
    }

      comb <- combn(union(names(df)[-1], mix$Marker), 2);
      for (i in 1:ncol(comb)) {
          m1 <- comb[1,i];
          m2 <- comb[2, i];
          if (levenshteinDistance(m1, m2) == 1) {
              warning <- append(warning, paste("Found two markers with very close names: did you mean", m1, "or", m2, "?"));
          }
      }
    }

    # Check that all marker names present in the mixture file are present in the frequency file
    if (length(error) == 0 && !all(mix$Marker %in% names(df))) {
        missingFreqs <- paste(setdiff(mix$Marker, names(df)), collapse=", ");
        error <- append(error, paste("The frequency database does not contain all markers present in the mixture file. The following are missing: ", missingFreqs));
    }

    if (length(error) > 0) {
        return(list(df=NULL, warning=NULL, error=error));
    }

    return(list(df=df, warning=warning, error=NULL));
}

