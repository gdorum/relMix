## ---- eval=FALSE--------------------------------------------------------------
#  relMixGUI()

## ---- out.width =  400, fig.retina = NULL,echo=FALSE,fig.cap="Figure 1: The RelMix main window."----
knitr::include_graphics("relMix_screenshot.png")

## ---- echo=FALSE, results='asis'----------------------------------------------
M <- read.table("mixture.txt",sep="\t",header=TRUE)
requireNamespace("pander",quietly=TRUE)
pander::pandoc.table(rbind(head(M, 5),data.frame(SampleName='...',Marker='...',Allele1='...',Allele2='...',Allele3='...')),missing="",caption="Table 1: The mixture profile in example 1.")

## ---- echo=FALSE, results='asis'----------------------------------------------
G <- read.table("references.txt",sep="\t",header=TRUE)
G2 <- rbind(G[1:2,],data.frame(SampleName='...',Marker='...',Allele1='...',Allele2='...'),G[22:24,],data.frame(SampleName='...',Marker='...',Allele1='...',Allele2='...'),G[44,])
rownames(G2) <- NULL
knitr::kable(G2,caption="Table 2: Reference profiles for mother and alleged father in example 1.")

## ---- out.width = '30%', fig.align='center', echo=FALSE, fig.cap='Figure 2: Database window.'----
knitr::include_graphics("database_crop.png")

## ---- echo=FALSE, results='asis'----------------------------------------------
freqs <- read.table("frequencies22Markers.txt",sep="\t",header=T,stringsAsFactors=F)
freqs2 <- freqs[11:18,1:5]
rownames(freqs2) <- NULL
pander::pandoc.table(cbind(rbind(freqs2,rep('...',5)),'...'=rep('...',9)),missing="",digits=4,caption="Table 3: Excerpt of the allele frequency database in example 1. Dots are used as decimal separator.")

## ---- out.width =  350, fig.retina = NULL,echo=FALSE,fig.cap="Figure 3: Database options. Theta value, silent allele frequency and minimum allele frequency."----
knitr::include_graphics("database_options.png")

## ---- out.width =  300, fig.retina = NULL,echo=FALSE,fig.cap="Figure 4: Define the pedigrees in example 1."----
knitr::include_graphics("pedigrees_ex1.png")

## ---- out.width =  150, fig.retina = NULL,echo=FALSE,fig.cap="Figure 5: Specify mother and child as the contributors in both pedigrees in example 1."----
knitr::include_graphics("contributors_ex1.png")

## ---- out.width =  150, fig.retina = NULL,echo=FALSE,fig.cap="Figure 6: In example 1 we assume dropout probabilities 0 for the mother and 0.1 for the child and drop-in 0.05."----
knitr::include_graphics("dropout_ex1.png")

## ---- out.width = 500, fig.retina = NULL,echo=FALSE,fig.cap="Figure 7: Results for example 1."----
knitr::include_graphics("results_ex1.png")

## ---- echo=FALSE, results='asis'----------------------------------------------
M <- read.table("mixture_silent_ex.txt",sep="\t",header=TRUE)
pander::pandoc.table(head(M, 5),missing="",caption="Table 4: Mixture file for example 2.")

## ---- echo=FALSE, results='asis'----------------------------------------------
G <- read.table("references_silent.txt",sep="\t",header=TRUE)
rownames(G) <- NULL
knitr::kable(G,caption="Table 5: Reference profile for the known contributor (C1) in example 2.")

## ---- echo=FALSE, results='asis'----------------------------------------------
freqs <- read.table("freqsSilent.txt",sep="\t",header=T,stringsAsFactors=F)
rownames(freqs) <- NULL
pander::pandoc.table(freqs,missing="",digits=4,caption="Table 6: Allele frequencies.")

## ---- out.width =  350, fig.retina = NULL,echo=FALSE,fig.cap="Figure 8: Database options for example 2."----
knitr::include_graphics("database_options_ex2.png")

## ---- out.width =  350, fig.retina = NULL,echo=FALSE,fig.cap="Figure 9: Scale frequencies or add rest allele."----
knitr::include_graphics("scaling_ex2.png")

## ---- out.width =  200, fig.retina = NULL,echo=FALSE,fig.cap="Figure 10: Mutations in example 2."----
knitr::include_graphics("mutations_ex2.png")

## ---- eval=FALSE--------------------------------------------------------------
#  persons <- c("C2","C1")
#  ped1 <- FamiliasPedigree(id=persons, dadid=c(NA,NA), momid=c("C1", NA),
#                           sex=c("male", "female"))

## ---- eval=FALSE--------------------------------------------------------------
#  persons <- c("C2","C1")
#  ped1 <- FamiliasPedigree(id=c(persons), dadid=c(NA, NA),
#                           momid=c( NA, NA), sex=c("male", "female"))

## ---- out.width =  300, fig.retina = NULL,echo=FALSE,fig.cap="Figure 11: We import custom pedigrees from R scripts in example 2."----
knitr::include_graphics("custom_ped_ex2.png")

## ---- out.width =  150, fig.retina = NULL,echo=FALSE,fig.cap="Figure 12: Tick both individuals as contributors to the mixture in both pedigrees."----
knitr::include_graphics("contributors_ex2.png")

## ---- out.width =  150, fig.retina = NULL,echo=FALSE,fig.cap="Figure 13: Possible dropout for both contributors in example 2."----
knitr::include_graphics("dropout_ex2.png")

## ---- out.width =  500, fig.retina = NULL,echo=FALSE,fig.cap="Figure 14: Computed LR for example 2."----
knitr::include_graphics("results_ex2.png")

