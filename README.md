
<!-- README.md is generated from README.Rmd. Please edit that file -->

# relMix

relMix makes relationship inference involving DNA mixtures with unknown
profiles and interprets DNA mixtures with related contributors. The main
function is the graphical user interface `relMixGUI`.

## Installation

Please follow the instructions here https://www.familias.name/openfamilias.html
to install relMix.


To provide pedigree plots, relMix uses the package `tkrplot`. However,
this package has some compatibility issues with MacOS and hence is not
included as a *hard* dependency. Users who wish to see pedigree plots in
the results screen have to install the package `tkrplot` manually with

``` r
install.packages("tkrplot")
```

The `tkrplot` package will be loaded by relMix automatially, so users do
not need to run `library("tkrplot")` in advance.
