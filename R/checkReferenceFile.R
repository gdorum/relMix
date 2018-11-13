# File checkReferenceFile.R
#
# Elias Hernandis <eliashernandis@gmail.com>
#
# Given a reference filepath this function attempts to load it.
# The second parameter is a dataframe with mixture data. It is used to compare
# Marker names and detect possible misspellings.
#
# It returns a list containing a dataframe, a list of warnings and a list of
# errors. If there are fatal errors, the dataframe will be FALSE. If there are
# fixable errors (see the README for details), the dataframe will contain the
# data from the mixture file but with those errors already fixed. The original
# data file will not be updated. In the event of ambiguous or possibly wrong
# marker names or sample names, the function will report them as warnings but
# not fix them automatically.


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
