# relMix 1.3

#### Major changes in v. 1.3

* Input files in ```relMixGUI``` are now being checked for most common errors, such as:
    + Incorrect table structure
    + Incorrect headers
    + Typos in marker and sample names
    + Non-numeric allele frequencies
    + Duplicate markers
  
* Four new functions, incorporated in ```relMixGUI```, provide the input check:
    + checkFrequenciesFile()
    + checkMixtureFile()
    + checkReferenceFile()
    + checkPedigreeFile()
   
##### Minor changes in v. 1.3.1

Package has been updated to use `gWidgets2` and `gWidgets2tcltk` instead of `gWidgets` and `gWidgetstcltk` because the latter are being removed from CRAN.

##### Minor changes in v. 1.3.2

* Replace deprecated argument `initialFilename` with updated name `initial.filename` in call to `gWidgets2::gfile` function.
* Fix a bug where database alleles were being all set to NA.
* Avoid saving loaded profiles when an error occurs in `f_importprof`.
* Fix an issue where if an allele was repeated in the mixture file, but not present in the database file, it would get added twice causing an error when Familias is called.

##### Minor changes in v. 1.3.3

* Make `tkrplot` a *suggests* dependency so that users which cannot install it can still use the non-plotting functionality of relMix. Now, if `tkrplot` is not available, the pedigrees will not be plotted in the results screen and instead a message explaining the problem will be shown.
