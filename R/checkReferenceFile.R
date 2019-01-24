#' Check a reference profiles file
#'
#' Given a reference profile file name the function attempts to load it and compare it to the mixture file to detect possible errors.
#' @param filename Path of the reference profiles file
#' @param mix Data frame with mixture data
#' @return A list containing
#' \itemize{
#'\item{\code{df}} {The loaded data frame, NULL if errors are present}
#' \item{\code{warning}} {A list of strings describing the errors that ocurred but could be fixed or that do not prevent the execution of the program.}
#' \item{\code{error}} {A list of strings describing the errors that ocurred that made it impossible to return a valid data frame.
#' If this list is not empty, then the data frame item will be null.}}
#' @details See the relMix vignette for a description of the format of the reference file. The data frame with mixture data is used to compare
#  marker names and detect possible misspellings.
#' If warnings are found, the function attempts to fix them and explains what it has done in the warning messages.
#' If an error is found, checking stops and a NULL dataframe is returned. The error is described in the error messages.
#' @seealso \code{\link{checkMixtureFile}} for information on how to load a mixture file.
#' @examples
#' \dontrun{
#' #Load a mixture file
#' mixfile <- system.file("extdata","mixture.txt",package="relMix")
#' mix <- checkMixtureFile(mixfile);
#' #Note: the mixture dataframe is passed as an argument. If the previous check failed,
#' #the program should not continue with the reference file check
#' reffile <- system.file("extdata","references.txt",package="relMix")
#' checkReferenceFile(reffile, mix$df);
#' }
#' @author Elias Hernandis
#' @export
checkReferenceFile <- function(filename, mix) {
    r <- commonChecks(filename, "reference file");
    df <- r$df;
    warning <- r$warning;
    error <- r$error;

    if (is.null(df)) {
        return(list(df=NULL, error=error, warning=NULL));
    }

    # check header column count
    if (!is.null(df) && ncol(df) != 4) {
        error <- "Incorrect number of columns"
    }

    referenceHeader <- c("SampleName", "Marker", "Allele1", "Allele2")

    if (length(error) == 0 && any(referenceHeader != names(df))) {
        warning <- append(warning, "Column titles must be written without spaces and with the first letter uppercased (fixed).");

        fixedHeader <- sapply(names(df), titleize);
        fixedHeader[1] <- "SampleName"; # special case for the first column header, which is "SampleName" instead of "Sample Name"

        if (all(fixedHeader != names(df))) {
            error <- append(error, paste("There are problems with the header row of the reference table. Please make sure it is the following: ", referenceHeader));
        } else {
            # if the errors are minor, fix the header for the user
            names(df) <- referenceHeader;
        }
    }

    # Kinship checking is deferred until the end

    # Check that marker names are the ones present in the mixture file
    if (!setequal(df$Marker, mix$Marker)) {
        markerError <- "The marker names in your reference file are different from those in your mixture file.";
        comb <- combn(union(df$Marker, mix$Marker), 2);
        for (i in 1:ncol(comb)) {
            m1 <- comb[1,i];
            m2 <- comb[2, i];
            if (levenshteinDistance(m1, m2) == 1) {
                markerError <- paste(markerError, "Found two markers with very close names: did you mean", m1, "or", m2, "?")
              }
        }
    error <- append(error, markerError);
    }

    #check that all profiles have equal number of rows. Otherwise return an error
    sampleNames <- unique(df$SampleName);
    if (length(error) == 0 && length(sampleNames) > 1) {
      if(abs(max(table(df$SampleName)) - min(table(df$SampleName)))>0){
        error <- paste("Profiles", paste(sampleNames,collapse=" and "), "have unequal number of rows.");
      }
    }


    # Check all allele reference data is numeric
    if (length(error) == 0 && !all(sapply(df[,3], is.numeric))) {
        error <- "There are values that are not numeric in the Allele1 column";
    }

    if (!all(sapply(df[,4], is.numeric))) {
        error <- append(error, "There are values that are not numeric in the Allele2 column");
    }

    if (length(error) > 0) {
        return(list(df=df, warning=NULL, error=error));
    }

    return(list(df=df, warning=warning, error=NULL));
}
