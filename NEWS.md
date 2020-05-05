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