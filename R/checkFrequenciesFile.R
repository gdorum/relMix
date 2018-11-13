# File checkFrequenciesFile.R
#
# Elias Hernandis <eliashernandis@gmail.com>
#
# Given a frequency database filepath this function attempts to load it.
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


library(sets)

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

    # Check that are alleles are numeric
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

