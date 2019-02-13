#' Check a pedigree file
#'
#' Given a pedigree file path the function attempts to load it and compare it
#' to the reference profiles to detect possible errors.
#'
#'@param filename Path of the pedigree file
#'@param df Data frame with reference profiles
#' @return A list containing
#' \itemize{
#' \item{df} {Pedigree, NULL if errors are present.}
#' \item{\code{warning}} {A list of strings describing the errors that ocurred but could be fixed or that do not prevent the execution of the program.}
#' \item{\code{error}} {A list of strings describing the errors that ocurred that made it imposible to return a valid data frame.
#' If this list is not empty, then the dataframe item will be null.}}
#' @details The pedigree file must be a .R file defining a pedigree (see the relMix vignette for an example). The data frame with reference data is used to compare names of individuals and detect possible misspellings.
#' If warnings are found, the function attempts to fix them and explains what it has done in the warning messages.
#' If an error is found, checking stops. The error is described in the error messages.
#' @examples
#' \dontrun{
#' #First load mixture file
#' mixfile <- system.file("extdata","mixture_silent_ex.txt",package="relMix")
#' mix <- checkMixtureFile(mixfile);
#' #Load reference file
#' reffile <- system.file("extdata","references_silent.txt",package="relMix")
#' ref <- checkReferenceFile(reffile, mix$df)
#' #Check pedigree file
#' pedfile <- system.file("extdata","custom_pedigree_maternity_duo.R",package="relMix")
#' checkPedigreeFile(pedfile,ref$df);
#' }
#' @author Elias Hernandis
#' @export

checkPedigreeFile <- function(filename, df) {
    warning <- c();
    error <- c();

    allowedKinships <- c("Father", "Mother", "Child")
    if (filename != "") {
        # load custom pedigrees
        #source(filename,local=TRUE);
      localEnv <- new.env()
      source(filename,local=localEnv)

    #     if(!"ped1"%in%ls(envir=localEnv)) {error <- append(error,"The pedigree file must define a pedigree of type 'FamiliasPedigree' in a variable named 'ped1'.")
    #     } else  allowedKinships <- get("ped1",envir=localEnv)$id;
    # }
      if(!"ped1"%in%ls(envir=localEnv)) {error <- append(error,"The pedigree file must define a pedigree of type 'FamiliasPedigree' in a variable named 'ped1'.")
      } else{
        ped1 <- get("ped1",envir=localEnv)
        allowedKinships <- ped1$id;
        #Make names of individuals in pedigree with first letter uppercase
        allowedKinships <- sapply(allowedKinships, titleize)
        ped1$id <- allowedKinships
        assign("ped1",ped1,envir=localEnv)
      }
    }



    # check the first "header" column for uppercase names
    if (length(error) == 0) {
      fixedHeaderColumn <- sapply(as.character(df[,1]), titleize)
        if (!all(fixedHeaderColumn %in%  allowedKinships)) {
            if (filename != "") {
                error <- append(error, "There are sample names in the reference file that are not defined in the custom pedigree file.");
            } else {
                 error <- append(error, "Sample names in reference profile must correspond to one of the standard (Father, Mother, Child). Alternatively, you may use a custom pedigree file.");
            }
        }

        # if (length(error) == 0 && any(fixedHeaderColumn != df[,1])) {
        #     warning <- append(warning, "Kinship names must be written with only the first letter uppercase, such as in 'Father' (fixed).");
        #     df[,1] <- fixedHeaderColumn;
        #
        # }
    }

    # # Check for very similar kinships
    # comb <- combn(df$Marker, 2);
    # for (i in 1:ncol(comb)) {
    #     m1 <- comb[1,i];
    #     m2 <- comb[2,i];
    #     if (levenshteinDistance(m1, m2) == 1) {
    #         warning <- append(warning, paste("Found two kinships with very close names in the reference file: did you mean", m1, "or", m2, "?"));
    #     }
    # }

    #return(list(warning=warning, error=error));

    if (length(error) > 0) {
      return(list(df=NULL, warning=NULL, error=error));
    }
    return(list(df=get("ped1",localEnv),warning=warning, error=error));
}
