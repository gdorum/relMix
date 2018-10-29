# File checkReferenceFile.R
#
# Elias Hernandis <eliashernandis@gmail.com>
#
# Given a pedigree filepath this function attempts to load it.
# The second parameter is a dataframe with reference data. It is used to compare
# Marker names and detect possible misspellings.
#
# It returns a list containing a dataframe, a list of warnings and a list of
# errors. If there are fatal errors, the dataframe will be FALSE. If there are
# fixable errors (see the README for details), the dataframe will contain the
# data from the mixture file but with those errors already fixed. The original
# data file will not be updated. In the event of ambiguous or possibly wrong
# marker names or sample names, the function will report them as warnings but
# not fix them automatically.


checkPedigreeFile <- function(filename, df) {
    warning <- c();
    error <- c();

    allowedKinships = c("Mother", "Father", "Child");
    if (filename != "") {
        # load custom pedigrees
        library(Familias);
        source(filename);
        allowedKinships <- pedigrees$id;
    }

    # check the first "header" column for uppercase names
    if (length(error) == 0) {
        fixedHeaderColumn <- sapply(as.character(df[,1]), titleize)
        if (!all(fixedHeaderColumn %in% allowedKinships)) {
            if (filename != "") {
                error <- append(error, "There are sample names in the reference file that are not defined in the custom pedigree file.");
            } else {
                error <- append(error, "Sample names in reference profile must correspond to one of the standard (Father, Mother, Child). Alternatively, you may use a custom pedigree file.");
            }
        }
        if (length(error) == 0 && any(fixedHeaderColumn != df[,1])) {
            warning <- append(warning, "Kinship names must be written with only the first letter uppercase, such as in 'Father' (fixed).");
            df[,1] <- fixedHeaderColumn;

        }
    }

    # Check for very similar kinships
    comb <- combn(df$Marker, 2);
    for (i in 1:ncol(comb)) {
        m1 <- comb[1,i];
        m2 <- comb[2,i];
        if (levenshteinDistance(m1, m2) == 1) {
            warning <- append(warning, paste("Found two kinships with very close names in the reference file: did you mean", m1, "or", m2, "?"));
        }
    }

    return(list(df=df, warning=warning, error=error));
}
