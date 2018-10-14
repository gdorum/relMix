# File util.R
#
# Elias Hernandis <eliashernandis@gmail.com>

# Given a string it returns it all lowercase with only the first letter of each
# word uppercase. It leaves non alphabetical characters untouched.
titleize <- function(x) {
    s <- strsplit(x, " ")[[1]];
    paste(toupper(substring(s, 1,1)), tolower(substring(s, 2)), sep="", collapse=" ")
}

# Given a TSV file this function performs some sanity checks on it, reporting
# fatal errors if they arise. See the README for details.
commonChecks <- function(filename, fileType) {
    error <- c();
    warning <- c();
    df <- tryCatch(
                   read.table(filename, header=TRUE, sep="\t", stringsAsFactors=FALSE),
                   error = function(e) { NULL });

    if (is.null(df)) {
        error <- append(error, paste("The", fileType, "you provided is not a well-formed tab-separated value (TSV) file. This can be avoided by generating the mixture file with a spreadsheet program instead of composing the file by hand."));
    }

    return(list(warning = warning, error = error, df = df));
}

# Given a pair of strings this function computes the Levenshtein distance
# between them, that is, an indication of the number of mistakes it takes to
# mistype one given the other is the correct one.
levenshteinDistance <- function(s1, s2) {
    m <- nchar(s1);
    n <- nchar(s2);
    t1 <- strsplit(s1, '')[[1]]
    t2 <- strsplit(s2, '')[[1]]
    d <- matrix(0, m+1, n+1);
    d[,1] <- matrix(0:m);
    d[1,] <- matrix(0:n);

    for (i in 1:m) {
        for (j in 1:n) {
            substCost <- 1;
            if (t1[i] == t2[j]) {
                substCost <- 0;
            }

            d[i+1,j+1] <- min(
                              d[i, j+1] + 1,
                              d[i+1, j] + 1,
                              d[i, j] + substCost
                              );
        }
    }
    return(d[m+1,n+1]);
}
