#Start the Familias package
#It must be already downloaded with e.g.
#install.packages(Familias)
library(Familias)
#Mother/child
persons <- c("C2","C1")
ped1 <- FamiliasPedigree(id=persons, dadid=c(NA,NA), momid=c("C1", NA),
                         sex=c("male", "female"))
