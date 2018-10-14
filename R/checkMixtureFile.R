# File checkMixtureFile.R
#
# Elias Hernandis <eliashernandis@gmail.com>
#
# Given a mixture filepath this function attempts to load it.
#
# It returns a list containing a dataframe, a list of warnings and a list of
# errors. If there are fatal errors, the dataframe will be FALSE. If there are
# fixable errors (see the README for details), the dataframe will contain the
# data from the mixture file but with those errors already fixed. The original
# data file will not be updated. In the event of ambiguous or possibly wrong
# marker names or sample names, the function will report them as warnings but
# not fix them automatically.

# source('util.R')

checkMixtureFile <- function(filename) {
    r <- commonChecks(filename, "reference file");
    df <- r$df;
    warning <- r$warning;
    error <- r$error;

    # check header column count
    if (!is.null(df) && ncol(df) != 5) {
        error <- "Incorrect number of columns"
    }

    referenceHeader <- c("SampleName", "Marker", "Allele1", "Allele2", "Allele3")

    if (length(error) == 0 && any(referenceHeader != names(df))) {
        warning <- append(warning, "Column titles must be written without spaces and with the first letter uppercased (fixed).");

        fixedHeader <- sapply(names(df), titleize);
        fixedHeader[1] <- "SampleName"; # special case for the first column header, which is "SampleName" instead of "Sample Name"

        if (all(fixedHeader != names(df))) {
            error <- append(error, paste("There are problems with the header row of the mixture table. Please make sure it is the following: ", referenceHeader));
        } else {
            # if the errors are minor, fix the header for the user
            names(df) <- referenceHeader;
        }
    }


    # Detect possible mistakes with similar SampleNames
    sampleNames <- unique(df$SampleName);
    if (length(sampleNames) > 1) {
        comb <- combn(sampleNames, 2)
        for (i in 1:ncol(comb)) {
            if (levenshteinDistance(comb[1,i], comb[2,i]) == 1) {
                warning <- append(warning, paste("Two very similar sample names were found in the mixture file. Did you mean", comb[1,i], "or", comb[2,i],"?"));
            }
        }
    }

    # Check all allele reference data is numeric
    if (length(error) == 0 && !all(sapply(df[,3], is.numeric))) {
        error <- "There are values that are not numeric in the Allele1 column of the mixture file";
    }

    if (!all(sapply(df[,4], is.numeric))) {
        error <- append(error, "There are values that are not numeric in hte Allele2 column of the mixture file");
    }

    if (!all(sapply(df[,5], is.numeric))) {
        error <- append(error, "There are values that are not numeric in hte Allele3 column of the mixture file");
    }

    if (length(error) > 0) {
        return(list(df=NULL, warning=NULL, error=error));
    }

    return(list(df=df, warning=warning, error=NULL));
}
