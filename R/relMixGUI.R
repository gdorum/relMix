#' GUI for relMix
#'
#' User-friendly graphical user interface for relMix.
#' @details Includes error checking for the input files.
#' @author Guro Dorum, Elias Hernandis
#' @seealso \code{\link{relMix}} for the main function implemented in \code{relMixGUI}.
#' @examples
#' #Examples can be found in the vignette and example data files can be found
#' #in the folder "inst\extdata" in the installation folder for relMix
#' @importFrom graphics plot
#' @importFrom utils read.table write.table packageVersion
#' @export
relMixGUI <- function(){

  options("guiToolkit"="tcltk")

  # Files must have '.' as decimal separator
  tableReader <- function(filename) {
    tab <- read.table(filename,header=TRUE,sep="\t",stringsAsFactors=FALSE,na.strings=c(NA,""))
    tryCatch( {  if(ncol(tab)==1) tab <- read.table(filename,header=TRUE,sep=",",stringsAsFactors=FALSE,na.strings=c(NA,"")) } ,error=function(e) e)
    tryCatch( {  if(ncol(tab)==1) tab <- read.table(filename,header=TRUE,sep=";",stringsAsFactors=FALSE,na.strings=c(NA,"")) } ,error=function(e) e)
    if(ncol(tab)==1) tab <- read.table(filename,header=TRUE,sep=";",stringsAsFactors=FALSE,na.strings=c(NA,""))
    return(tab) #need dataframe to keep allele-names correct!!
  }

  tableWriter <- function(filename,obj){
    write.table(obj,file=filename,sep="\t",row.names=FALSE,quote=FALSE)
  }

  # Receives a the output of an error checking function and displays
  # the appropriate error messages. Returns the data frame on
  # successfull load or NULL on error.
  process_errors <- function(result) {
      message = "";
      if (length(result$error) > 0) {
          message <- paste(message, "Error(s):\n", result$error, sep="\n");
          message <- paste(message, "\nThese errors are fatal. You will need to fix the file yourself.", sep="\n");
          f_errorWindow(message);
          return(NULL); # Terminate early if there are errors, returning a NULL dataframe
      }

      if (length(result$warning) > 0) {
          message <- paste(message, "Warning(s):\n", result$warning, sep="\n");
          message <- paste(message, "\nThese warnings are not fatal. You may continue using the program but please be aware that results may be incorrect.", sep="\n");
          f_warningWindow(message)
      }

      return(result$df);
  }

  f_importprof <- function(h,...) {
    type=h$action #get type of profile

    # enforce loading files in order to facilitate error checking
    if (type=="reference") {
      if(!exists('mixture',envir=mmTK)) {f_errorWindow("Import mixture profile before reference profiles"); return();}
    }
    if (type == "frequencies") {
      if(!exists('mixture',envir=mmTK)) {f_errorWindow("Import mixture profile before frequencies"); return();}
      if(!exists('reference',envir=mmTK)) {f_errorWindow("Import reference profiles before frequencies"); return();}
    }

    proffile = gWidgets2::gfile(text=paste("Open ",type," file",sep=""),type="open",
                     filter=list("text"=list(patterns=list("*.txt","*.csv","*.tab")),"all"=list(patterns=list("*"))))

    # check errors
    tryCatch({
    if (type == "mixture") {
      Data <- process_errors(checkMixtureFile(proffile));
    } else if (type == "reference") {
      # get mixture file from environment
      mix <- get("mixture",envir=mmTK);
      Data <- process_errors(checkReferenceFile(proffile, mix));
    } else if (type == "frequencies") {
      mix <- get("mixture",envir=mmTK);
      Data <- process_errors(checkFrequenciesFile(proffile, mix));
    }}, error = function(e) {
      # Handle fatal errors due to bad file formats
      # For example, if when loading a mixture file a frequency file is provided, the error will be caught here
      f_errorWindow(paste("There was an error loading the file. It does not look like a ", type, " file. Please make sure the file is correct and try again."))
    })
    assign(h$action,Data,envir=mmTK) #save object
  }

  f_export <- function(obj) {
    savefile = gWidgets2::gfile(text=paste("Save file as",sep=""),type="save", initialfilename = "LR.txt")
    #filter=list("text"=list(patterns=list("*.txt","*.csv","*.tab")),"all"=list(patterns=list("*")))
    tableWriter(savefile,obj) #load profile
  }

  # Make pedigree
  f_pedigree <- function(h,...){
    if(gWidgets2::svalue(h$obj)=="Paternity"){
      #Define the persons involved in the case
      persons <- c("Mother", "Father", "Child")
      sex <- c("female", "male", "male")
      #Pedigree files exported from Familias use notation "ped1"
      ped1 <- Familias::FamiliasPedigree(id=persons, dadid=c(NA,NA,"Father"), momid=c(NA,NA,"Mother"), sex=c("female", "male", "male"))
    } else if(gWidgets2::svalue(h$obj)=="Unrelated"){
      #Define the persons involved in the case
      persons <- c("Mother", "Father", "Child")
      sex <- c("female", "male", "male")
      ped1 <- Familias::FamiliasPedigree(id=persons, dadid=c(NA,NA,NA), momid=c(NA,NA,"Mother"), sex=c("female", "male", "male"))
    } else {
      pedfile = gWidgets2::gfile(text=paste("Open pedigree R file",sep=""),type="open",
                      filter=list("text"=list(patterns=list("*.R")),"all"=list(patterns=list("*"))))
      ped1 <- process_errors(checkPedigreeFile(pedfile,get("reference",mmTK)))
      #source(pedfile)
      #persons <- ped1$id #Pedigree files exported from Familias use notation "ped1"
      persons <- ped1$id
    }
    assign(paste("ped",h$action,sep=""),ped1,envir=mmTK)
    assign(paste("persons_ped",h$action,sep=""),persons,envir=mmTK)
  }

  # f_pedigree <- function(h,...){
  #
  #   if(gWidgets2::svalue(h$obj)=="Paternity"){
  #     #Define the persons involved in the case
  #     #If pedigree 2, use same individuals as in first pedigree
  #     if(gWidgets2::svalue(h$action)=='2'){
  #       firstPed <- get('ped1',envir=mmTK)
  #       persons1 <- get('persons_ped1',envir=mmTK)
  #       persons <- persons1
  #     } else{
  #       persons <- c("Mother", "Father", "Child")
  #       sex <- c("female", "male", "male")
  #     }
  #     #persons <- c("Mother", "Father", "Child")
  #     #sex <- c("female", "male", "male")
  #     ped1 <- Familias::FamiliasPedigree(id=persons, dadid=c(NA,NA,"Father"), momid=c(NA,NA,"Mother"), sex=c("female", "male", "male"))
  #   } else if(gWidgets2::svalue(h$obj)=="Unrelated"){
  #     #Define the persons involved in the case
  #     if(gWidgets2::svalue(h$action)=='2'){
  #       firstPed <- get('ped1',envir=mmTK)
  #       persons1 <- get('persons_ped1',envir=mmTK)
  #       persons <- persons1
  #     } else{
  #       persons <- c("Mother", "Father", "Child")
  #       sex <- c("female", "male", "male")
  #     }
  #     #persons <- c("Mother", "Father", "Child")
  #     #sex <- c("female", "male", "male")
  #     ped1 <- Familias::FamiliasPedigree(id=persons, dadid=c(NA,NA,NA), momid=c(NA,NA,"Mother"), sex=c("female", "male", "male"))
  #   } else {
  #     pedfile = gWidgets2::gfile(text=paste("Open pedigree R file",sep=""),type="open",
  #                     filter=list("text"=list(patterns=list("*.R")),"all"=list(patterns=list("*"))))
  #     source(pedfile)
  #     if(!all(exists("persons"),exists("ped1"))) f_errorWindow("File should define both pedigree and persons")#stop("File should define both pedigree and persons")
  #   }
  #   assign(paste("ped",h$action,sep=""),ped1,envir=mmTK)
  #   assign(paste("persons_ped",h$action,sep=""),persons,envir=mmTK)
  # }

  #Get values from object and assign to new object in mmTK environment
  f_values <- function(h,...) {
    r <- gWidgets2::svalue(h$obj)
    assign(h$action,r,envir=mmTK)
  }
  #Pop-up error window
  f_errorWindow <- function(message){
    # errorWindow <- gWidgets2::gwindow("Error", )
    # gWidgets2::glabel(message,container=errorWindow, expand=TRUE)
    # gWidgets2::gbutton("ok", container=errorWindow,handler=function(h,...){
    #   gWidgets2::dispose(h$obj)
    # })
    gWidgets2::gmessage(message, title="Error", icon="error")
  }

  f_warningWindow <- function(message) {
    gWidgets2::gmessage(message, title="Warning", icon="warning")
  }

  #Makes pop-up window for database
  f_database <- function(h,...){
    freqWindow <- gWidgets2::gwindow("Frequency database", height = 50, width = 250, visible=FALSE)
    freqGroup2 <- gWidgets2::ggroup(horizontal = FALSE, container=freqWindow,spacing=7)
    gWidgets2::gbutton(text="Import allele frequencies",container=freqGroup2,handler=f_importprof,action="frequencies")
    optButton <- gWidgets2::gbutton("Options", handler=f_options, container=freqGroup2)
    saveButton <- gWidgets2::gbutton(text = 'ok',container=freqGroup2, handler = function(h,...){

      #Get frequencies
      if(!exists('frequencies',envir=mmTK)) {f_errorWindow("Alle frequencies not imported")}
      # TODO remove these
      if(!exists('mixture',envir=mmTK)) {f_errorWindow("Import mixture profile before frequencies")}
      if(!exists('reference',envir=mmTK)) {f_errorWindow("Import reference profiles before frequencies")}
      A <- get("frequencies",envir=mmTK)
      M <- get("mixture",envir=mmTK)
      G <- get("reference",envir=mmTK)

      freqs <- A[,-1,drop=FALSE]
      colnames(freqs) <- colnames(A)[2:ncol(A)]
      rownames(freqs) <- A[,1]
      optPar <- get('optPar',envir=mmTK)

      ix1 <- M[,2]%in%colnames(freqs)
      ix2 <- G[,2]%in%colnames(freqs)
      m <- unique(M[!ix1,2],G[!ix1,2])
      if(length(m)>0) {
        errorWindow <- gWidgets2::gwindow("Error")
        gWidgets2::glabel(paste("Markers not found in database:",paste(m,collapse=", ")),container=errorWindow)
        stop(paste("Markers not found in database:",paste(m,collapse=", ")))
        #f_errorWindow(paste("Markers not found in database:",paste(m,collapse=", ")))
      }


      freqsS <- f_silent(freqs,optPar$silent) #Add silent allele
      freqsS <- get('freqsS',envir=mmTK)
      db2 <- f_unobserved(freqsS,M,G,optPar$MAF) #Add unobserved alleles

      gWidgets2::dispose(freqWindow)
    }) #end handler savebutton

    gWidgets2::visible(freqWindow) <- TRUE
  }


  #Makes pop-up window to fill in mutation details
  f_mutations <- function(h,...){
    mutWindow <- gWidgets2::gwindow("Mutations",visible=FALSE)
    mutGroup2 <- gWidgets2::ggroup(horizontal = FALSE, container=mutWindow,spacing=7)
    mutFrame2 <- gWidgets2::gframe("Mutation model",container=mutGroup2)
    mutFrame3 <- gWidgets2::gframe("Range",container=mutGroup2)
    mutPar <- get('mutPar',mmTK)
    #Range window
    objRange <- gWidgets2::gedit(mutPar$mutRange, width=4, container=mutFrame3)
    gWidgets2::enabled(objRange) <- FALSE
    #Common male and female mutation model
    Mut <- get('Mut',mmTK)
    mutModels <- c("Equal","Proportional","Stepwise")
    comboMut <- gcombobox(mutModels, selected=which(mutModels==Mut),
                          container=mutFrame2, handler=function(h,...) {
                            #assign("Mut",gWidgets2::svalue(h$obj),envir=mmTK)
                            #Range will not be used unless model is Stepwise
                            #assign("range1",mutPar$mutRange,envir=mmTK)
                            f_changeMutation(gWidgets2::svalue(h$obj),objRange)})
    #List to store parameters in
    objMut <- vector('list',3)
    names(objMut) <- c("mutRange","fMutRate","mMutRate")
    objMut[[1]] <- objRange
    mutFrame4 <- gWidgets2::gframe("Mutation rates",container=mutGroup2,horizontal=TRUE)
    #Separate male and female mutation rates
    gWidgets2::glabel("Female",container=mutFrame4)
    objMut[[2]] <- gWidgets2::gedit(mutPar$fMutRate, width=5, container=mutFrame4)
    gWidgets2::glabel("Male",container=mutFrame4)
    objMut[[3]] <- gWidgets2::gedit(mutPar$mMutRate, width=5,container=mutFrame4)
    saveButton <- gWidgets2::gbutton(text = "Save",container=mutGroup2, handler = function(h,...) {
      assign("Mut",gWidgets2::svalue(comboMut),envir=mmTK)
      saveParameters(objMut,mutWindow,"mutPar")}
    )
    gWidgets2::visible(mutWindow) <- TRUE
  }

  #Set mutation model
  f_changeMutation <- function(mutMod,objRange){
    if(mutMod=="Stepwise"){ #If stepwise, also give option to give a mutation range
      #Common male and female mutation range
      gWidgets2::enabled(objRange) <- TRUE
    }
    if(mutMod=="Equal"){
      gWidgets2::enabled(objRange) <- FALSE
    }
    if(mutMod=="Proportional"){
      gWidgets2::enabled(objRange) <- FALSE
    }
    #if(mutMod=="Custom"){ #If custom, read in mutation model from file
    #  gWidgets2::enabled(objRange) <- FALSE
    #  enable(objRateF)<- FALSE
    #  enable(objRateM) <- FALSE
    #  mutfile = gWidgets2::gfile("Open mutation R file",type="open",
    #                  filter=list("text"=list(patterns=list("*.R",".txt")),"all"=list(patterns=list("*"))))
    #  MM <- read.table(mutfile)
    #  assign(paste(h$action,"Mat",sep=""),MM,envir=mmTK)
    #}
  }

  #Makes option pop-up window to fill in theta and silent alleles
  f_options <- function(h,...){
    optPar <- get("optPar",mmTK)
    optWindow <- gWidgets2::gwindow("Options",visible=FALSE)
    optGroup <- gWidgets2::ggroup(horizontal = FALSE, container=optWindow,spacing=7)
    #Theta
    gWidgets2::glabel("Theta",container=optGroup)
    objOpt <- vector('list',3)
    names(objOpt) <- c("theta","silent","MAF")
    objOpt[[1]] <- gWidgets2::gedit(optPar$theta, width=5, container=optGroup)
    #Silent allele
    gWidgets2::glabel("Silent allele frequency",container=optGroup)
    objOpt[[2]] <- gWidgets2::gedit(optPar$silent, width=5, container=optGroup)
    #Minimum allele frequency
    gWidgets2::glabel("Min. allele frequency",container=optGroup)
    objOpt[[3]] <- gWidgets2::gedit(optPar$MAF, width=5, container=optGroup)
    saveButton <- gWidgets2::gbutton(text = "Save",container=optGroup, handler = function(h,...){
      saveParameters(objOpt,optWindow,"optPar")
    })
    gWidgets2::visible(optWindow) <- TRUE
  }

  #Makes pop-up window to fill in dropout values per contributor and drop-in
  f_dropout <- function(h,...){
    dropWindow <- gWidgets2::gwindow("Dropout and drop-in",width=500,visible=FALSE)
    dropGroup <- gWidgets2::ggroup(horizontal = FALSE, container = dropWindow)
    #Get pedigree data
    if(!exists("idxC1",envir=mmTK)){
      f_errorWindow("Specify contributors first")
    }
    dropFrame2 <- gWidgets2::gframe("Dropout for each contributor",container=dropGroup)
    #Get contributors
    dropliste <- get("dropliste",envir=mmTK)
    #Get default values/values that have already been specified
    idxC1 <- get("idxC1",envir=mmTK)
    idxC2 <- get("idxC2",envir=mmTK)
    idxC <- union(idxC1,idxC2)
    #If dropliste has not been updated yet, it contains only one default dropout value
    #and a drop-in value. Set the default dropout value for all contributors
    #if(length(dropliste)==2) dropliste <- c(rep(dropliste[1],length(idxC)),dropliste[2])
    dropNames <- names(dropliste)
    objDrop <- list()
    for(i in 1:length(idxC)){
      g <- gWidgets2::glabel(idxC[i], container=dropFrame2,horizontal=FALSE)
      #Check if a dropout value is already specified, if not set to 0
      val <- ifelse(idxC[i]%in%dropNames,dropliste[[which(dropNames==idxC[i])]],0)
      objDrop[[i]] <- gWidgets2::gedit(format(val,digits=4),width=4,container=dropFrame2)
    }
    dropinFrame <- gWidgets2::gframe("Drop-in",container=dropGroup)
    #objDrop[[length(dropliste)]] <- gspinbutton(0,1,by=0.01, value=dropliste[[length(dropliste)]],digits=3,container=dropinFrame)
    objDrop[[length(idxC)+1]] <- gWidgets2::gedit(format(dropliste[[length(dropliste)]],digits=4),width=4,container=dropinFrame)
    names(objDrop) <- c(idxC,"dropin")
    saveButton <- gWidgets2::gbutton(text = "Save",container=dropGroup, handler = function(h,...) {
      saveParameters(objDrop,dropWindow,"dropliste")
    })
    gWidgets2::visible(dropWindow) <- TRUE
  }
  #Save button
  saveParameters <- function(objList,window,varName){
    values <- list()
    for(i in 1:length(objList)){
      #If comma used as decimal, replace with dot
      v <- as.numeric(sub(",",".",gWidgets2::svalue(objList[[i]])))
      if(length(v)==0) {f_errorWindow("Specify all values"); stop()}
      else if(is.na(v) | (v<0 | v>1) ) {f_errorWindow("Need values in the range 0.0 and 1.0"); stop()}
      else values[[i]] <- v
    }
    names(values) <- names(objList)
    assign(varName,values,envir=mmTK)
    gWidgets2::dispose(window)
  }

  f_contributors <- function(h,...){
    #Get pedigree data
    if(!exists("persons_ped1",envir=mmTK)){
      f_errorWindow("Specify pedigrees first")
    }
    persons1 <- get("persons_ped1",envir=mmTK)
    persons2 <- get("persons_ped2",envir=mmTK)
    if(!identical(persons1,persons2)) {
      f_errorWindow("Persons in pedigree 1 and 2 must match!")
    } else{
      contWindow <- gWidgets2::gwindow("Contributors")
      contGroup <- gWidgets2::ggroup(horizontal = FALSE, container = contWindow)
      contFrame <- gWidgets2::gframe("Specify contributors in mixture",container=contGroup)
      #Check if contributors have already been specified. If so, set these as checked boxes
      if(!exists("idxC1",envir=mmTK)) idxC1 <- NULL
      if(!exists("idxC2",envir=mmTK)) idxC2 <- NULL
      else {
        idxC1 <- get("idxC1",envir=mmTK)
        idxC2 <- get("idxC2",envir=mmTK)
      }
      persons <- union(persons1,persons2)

      #Different contributors under each hypothesis
      contFrame <- gWidgets2::gframe("Contributors in pedigree 1",container=contGroup)
      gWidgets2::gcheckboxgroup(persons, checked=persons1%in%idxC1,horizontal=FALSE,container=contFrame,
                     handler=f_values,action="idxC1")
      contFrame2 <- gWidgets2::gframe("Contributors in pedigree 2",container=contGroup)
      gWidgets2::gcheckboxgroup(persons, checked=persons2%in%idxC2,horizontal=FALSE,container=contFrame2,
                     handler=f_values,action="idxC2")

      gWidgets2::gbutton("ok", container=contGroup,handler=function(h,...){
        gWidgets2::dispose(h$obj)
      })
    }
  }

  #Get mixture in the right format
  f_mixture <- function(E){
    m <- length(unique(E$Marker)) #Number of markers
    n <- length(unique(E$SampleName)) #Number of samples (if replicates)
    #Split according to markers and then sample
    mix <- split(E,E$Marker)
    mix2 <- lapply(mix,function(x) split(x[,-c(1,2)],x$SampleName))
    #Remove NA's
    lapply(mix2,function(z) lapply(z, function(z2) z2[!is.na(z2)]))
  }
  #Reads genotypes of known contributors and stores in list
  f_genotypes <- function(GT){
    m <- length(unique(GT$Marker)) #Number of markers
    n <- length(unique(GT$SampleName)) #Number of contributors
    gt <- split(GT,GT$Marker) #Split according to markers
    lapply(gt,function(x) split(x[,3:4],x$SampleName)) #Split according to individual
  }

  #Add silent allele to database
  f_silent <- function(freqs,ps){

    if(ps>0) {
      freqsS <- rbind(freqs,rep(as.numeric(ps),ncol(freqs)))
      rownames(freqsS) <- c(rownames(freqs),'Silent')
    } else freqsS <- freqs
    assign('freqsS',freqsS,envir=mmTK)
  }

  #Add alleles not in database with frequency MAF
  #Checks if any alleles are below MAF (by running f_MAF)
  f_unobserved <- function(freqs,M,G,MAF){

    alleleNames <- rownames(freqs)
    markerNames <- colnames(freqs)

    n <- ncol(M)
    keepIx <- alNotDB <- numeric()
    #Go through each marker in mixture profile to look for alleles not in database
    #Simultaneously checking in reference profile (assumes that markers in reference and
    #mixture are the same)
    for(i in 1:nrow(M)){
      #Get all database frequencies for the marker
      mark <- M[i,2]
      ix <- which(mark==markerNames)
      #Check if marker name is the same in mixture file and database, otherwise give error
      #(This has already been checked once the data was imported in)
      if(length(ix)==0) f_errorWindow(paste("Marker",mark,"not found in database"))
      al <- alleleNames[!is.na(freqs[,ix])]
      #Alleles in mixture
      am <- M[i,3:n][!is.na(M[i,3:n])]
      #Alleles in genotypes
      ag <- unlist(G[G[,2]==mark,3:4])
      aa <- unique(c(am,ag))
      #Check if all alleles in mixture and reference exist in database
      idx <- c(aa%in%al)

      if(any(!idx)){
        cat(i,"\n")
        keepIx <- c(keepIx,rep(ix,sum(!idx))) #Index of marker
        alNotDB <- c(alNotDB,aa[!idx]) #Alleles not found in db
      }
    }

    #Change format of database before adding new alleles
    db <- numeric()
    for(i in 1:ncol(freqs)){
      ix <- which(!is.na(freqs[,i]))
      a <- alleleNames[ix]
      f <- freqs[ix,i]
      Marker <- rep(markerNames[i],length(ix))
      dbNew <- data.frame(Marker,a,f)
      colnames(dbNew) <- c("Marker","Allele","Frequency")
      db <- rbind(db,dbNew)
    }

    if(length(alNotDB)>0) { #There are alleles not in database

      mess <- paste("Allele",alNotDB, "will be added to marker",markerNames[keepIx], "with frequency",MAF)
      gWidgets2::gmessage(mess, title="Note",icon = "info")
      #Add new allele at the end of database
      newData <- data.frame(markerNames[keepIx],alNotDB,MAF)
      colnames(newData) <- c('Marker','Allele','Frequency')
      db <- data.frame(Marker=c(as.character(db$Marker),as.character(newData$Marker)),
                       Allele=c(as.character(db$Allele),as.character(newData$Allele)),
                       Frequency=c(db$Frequency,newData$Frequency))
      assign('db',db,mmTK)
      #Check for MAF, scale and sort
      f_MAF(db,MAF)
      # }) #end handler add alleles

    } else{ #No alleles added to database
      #Check for MAF, scale and sort
      assign('db',db,mmTK)
      f_MAF(db,MAF)
    }
  }# end f_unobserved

  #Check if any allele frequencies are below MAF and sets frequency to MAF
  #Scales frequencies if necessary (with f_MAF)
  #Function used by f_unobserved
  f_MAF <- function(db,MAF){
    db <- get('db',mmTK)
    if(any(db$Frequency<MAF)){ #Frequencies below MAF
      gWidgets2::gconfirm("Some frequencies are below the min. allele frequency.
               Change the indicated frequencies?", title="Note",icon = "question",handler = function(h,...){
                 db$Frequency[db$Frequency<MAF] <- MAF
                 assign('db',db,mmTK)
    })

    }
    #Check if scaling is necessary
    db <- get('db',mmTK)
    f_scale(db)
  }


  #Check that all frequencies sum to 1, otherwise scale
  #Sort database and assign final database to environment
  #Function used by f_MAF
  f_scale <- function(db){
    db <- get('db',mmTK)
    markerNames <- unique(db[,1])

    #Sort database according to marker, then allele. Silent allele last
    #First reorder levels of Allele
    aL <- levels(db[,2])
    if(any(aL=='Silent')) { db$Allele <- factor(db$Allele,c(sort(as.numeric(aL[which(!aL=='Silent')])),'Silent'))
    } else db$Allele <- factor(db$Allele,sort(as.numeric(aL)))

    db <- db[order(db$Marker,db$Allele),]
    #Assign database that will be used if scaling is not done
    assign('dbF',db,envir=mmTK) #Final database

    #Check that frequencies sum to 1, otherwise scale
    sums <- sapply(1:length(markerNames),function(i) sum(db[db[,1]==markerNames[i],3]))
    ix <- which(sums!=1)
    if(length(ix)>0) {

      w <- gWidgets2::gconfirm(format("Frequencies do not sum to 1. Do you want to scale? If not, a rest allele will be added.",jusity="centre"), title="Note",icon = "info")
      if(w){ #Scale
        for(m in markerNames[ix]){
            db[db[,1]==m,3] <- db[db[,1]==m,3]/sum(db[db[,1]==m,3])
        }
      } else{ #Rest allele, or scale if frequencies sum > 1
            for(i in ix){
                if(sums[i]>1) { #Enforce scaling
                  db[db[,1]==markerNames[i],3] <- db[db[,1]==markerNames[i],3]/sum(db[db[,1]==markerNames[i],3])
                } else { #Rest allele
                    db <- rbind(db,data.frame(Marker=markerNames[i],Allele='r',Frequency=1-sums[i]))
                }
            }
        if(any(sums[ix]>1)) {
          mess <- paste("Enforced scaling for markers:",paste(markerNames[sums[ix]>1],collapse=", "))
          gWidgets2::gmessage(mess, title="Note",icon = "info")
        }
      }
      assign('dbF',db,envir=mmTK) #Final database
    }
  }



  f_LR <- function(){

    ####### Get input ########
    #Mutation parameters
    mutPar <- get("mutPar",envir=mmTK)
    r1 <- mutPar$fMutRate
    r2 <- mutPar$mMutRate
    range1 <- mutPar$mutRange
    mutModel <- get("Mut",envir=mmTK)
    optPar <- get('optPar',mmTK)
    theta <- optPar$theta

    #Profiles
    if(!exists('mixture',envir=mmTK)){ f_errorWindow("Mixture not imported"); stop()}
    if(!exists('reference',envir=mmTK)){ f_errorWindow("Reference profiles not imported"); stop()}
    if(!exists('frequencies',envir=mmTK)){ f_errorWindow("Missing allele frequencies"); stop()}
    if(!exists("ped1",envir=mmTK)){ f_errorWindow("Missing pedigree information"); stop()}
    if(!exists("ped2",envir=mmTK)){ f_errorWindow("Missing pedigree information"); stop()}
    if(!exists("idxC1",envir=mmTK)){ f_errorWindow("Contributors not specified"); stop()}


    E <- get("mixture",envir=mmTK) #get object
    G <- get("reference",envir=mmTK) #get object
    G$SampleName <- sapply(G$SampleName,titleize)
    #Remove AMEL marker? Or not allow for it?
    R <- f_mixture(E)
    knownGenos <- f_genotypes(G)
    #Frequencies
    db2 <- get("dbF",envir=mmTK)
    alleleNames <- as.character(db2[,2])


    #Check if there are non-numeric allele names and stepwise mutation model
    #If so, give a warning
    isNumeric <- suppressWarnings(as.numeric(alleleNames[alleleNames!="Silent"]))
    if(any(is.na(isNumeric)) && mutModel=='Stepwise'){
      f_errorWindow(format("Stepwise mutation model requires numeric allele names.
                    Change allele names or choose a different mutation model.",justify="centre"))
    }

    #Make lists of alleles and frequencies per marker
    allelesAll <- split(db2$Allele,db2$Marker)
    afreqAll <- split(db2$Frequency,db2$Marker)

    #Pedigree data
    ped1 <- get("ped1",envir=mmTK)
    ped2 <- get("ped2",envir=mmTK)
    persons1 <- get("persons_ped1",envir=mmTK)
    persons2 <- get("persons_ped2",envir=mmTK)
    pedigrees <- list(ped1,ped2)
    #idxC <- get("idxC",envir=mmTK) #Individuals who are contributors
    idxC1 <- get("idxC1",envir=mmTK) #Contributors under H1
    idxC2 <- get("idxC2",envir=mmTK) #Contributors under H2
    idxK <- unique(G[,1]) #Individuals with known genotypes
    idxC <- union(idxC1,idxC2)
    idxU <- idxC[!idxC%in%idxK] #Contributors with uknown genotypes. Assuming that all individuals are represented in both pedigrees!

    # #Check if the names of individuals in the reference file correspond to names in pedigree
    # notFound <- unique(G$SampleName)[!unique(G$SampleName)%in%persons1]
    # if(length(notFound)>0) { f_errorWindow(c("Following individuals in reference file not found in pedigree:",paste(notFound,collapse=", "))); stop() }

    #Dropout/drop-in
    #Set default dropout 0 for all contributors if none is set
    # if(!exists('dropliste',envir=mmTK)){
    #   dDef <- get('drop',envir=mmTK)
    #   D <- rep(list(dDef$d),length(idxC))
    #   names(D) <- idxC
    #   di <- dDef$di
    # } else{
    drop <- get('dropliste',mmTK)
    if(!all(idxC%in%names(drop))){ f_errorWindow("Specify dropout and drop-in"); stop()
    } else {
      D <- drop[idxC]
      di <- drop[["dropin"]]
    }


    ############# Computations ###########
    #infoWindowLR <- gWidgets2::gwindow("",visible=TRUE)
    #gWidgets2::glabel("Computing LR...",container=infoWindowLR)

    markers <- names(R)
    LRmarker <- numeric(length(markers))
    lik1 <- lik2 <- numeric(length(markers))
    for(i in 1:length(R)){

      #Create locus object for each marker
      alleles <- as.character(allelesAll[[which(names(allelesAll)==markers[i])]])
      afreq <- afreqAll[[which(names(afreqAll)==markers[i])]]
      #Set MutationRate2 to a very small number to avoid probability 0 of mutation to and from microinvariants
      #This is same as is done in Familias (MutationRate2 will be ignored unless mutation model is stepwise)
      locus <- Familias::FamiliasLocus(frequencies=afreq,name=markers[i],
                             allelenames=alleles, MutationModel=mutModel, femaleMutationRate=r1,
                             maleMutationRate=r2, MutationRange=range1, MutationRate2 = 1e-06)

      #If there is a silent allele we need to modify the mutation matrix.
      #Silent allele must be given as 's' (not 'silent' as in Familias)
      #That way Familias will treat it like a regular allele,
      #while relMix will treat is specially
      if('Silent'%in%alleles){
        newAlleles <- c(alleles[-length(alleles)],'s')
        mm <- locus$femaleMutationMatrix #Assuming same mutation matrix for male and female
        colnames(mm) <- rownames(mm) <- newAlleles
        locus <- Familias::FamiliasLocus(frequencies=afreq,name=markers[i],
                               allelenames= newAlleles, MutationModel='Custom', MutationMatrix=mm)
      }

      #datamatrix <- createDatamatrix(locus,knownGenos[[which(names(knownGenos)==markers[i])]],idsU=idxU)
      names(pedigrees) <- c("H1","H2")

      if(identical(idxC1,idxC2)){
        #Find genotypes for known and unknown individuals involved in the hypothesis
        datamatrix <- createDatamatrix(locus,knownGenos[[which(names(knownGenos)==markers[i])]],idsU=idxU)
        res <- relMix(pedigrees, locus, R=R[[i]], datamatrix, ids=idxC, D=lapply(D[idxC],function(x) c(x,x^2)),di=di, kinship=theta)
        lik1[i] <- res$H1
        lik2[i] <- res$H2

      } else {
        #relMix must be run twice with different contributors specified
        #Find genotypes for known and unknown individuals involved in each hypothesis
        datamatrix1 <- createDatamatrix(locus,knownGenos[[which(names(knownGenos)==markers[i])]],idsU=intersect(idxU,idxC1))
        datamatrix2 <- createDatamatrix(locus,knownGenos[[which(names(knownGenos)==markers[i])]],idsU=intersect(idxU,idxC2))
        res1 <- relMix(pedigrees$H1, locus, R=R[[i]], datamatrix1, ids=idxC1, D=lapply(D[idxC1],function(x) c(x,x^2)),di=di, kinship=theta)
        res2 <- relMix(pedigrees$H2, locus, R=R[[i]], datamatrix2, ids=idxC2, D=lapply(D[idxC2],function(x) c(x,x^2)),di=di, kinship=theta)
        lik1[i] <- res1[[1]]
        lik2[i] <- res2[[1]]
      }
    }
    LRmarker <- lik1/lik2
    #Set NaN's to 0 (likelihood 0 also under H2)
    LRmarker[is.na(LRmarker)] <- 0

    Data <- data.frame(Marker=markers,LR=LRmarker,LikH1=lik1,LikH2=lik2,stringsAsFactors=FALSE)
    #Make result window
    LRwindow <- gWidgets2::gwindow("Results",height=800,width=800,visible=FALSE)
    LRwindowGroup <- gWidgets2::ggroup(horizontal = TRUE, container = LRwindow)

    #### left col ####
    paramGroup <- gWidgets2::ggroup(horizontal=FALSE,container=LRwindowGroup)
    #gWidgets2::size(paramGroup) <- c(250,450)
    #Database options
    databaseData <- data.frame(Parameter=c("theta","Silent allele freq.","Min. allele freq."),Value=unlist(optPar))
    databaseData[,1] <- as.character(databaseData[,1])
    databaseData[,2] <- as.character(databaseData[,2])

    databaseGroup <- gWidgets2::ggroup(horizontal=FALSE,container=paramGroup)
    databaseFrame <- gWidgets2::gframe("Database",container=databaseGroup,horizontal=FALSE)
    #dataTab <- gWidgets2::gtable(databaseData,container=databaseFrame)
    dataTab <- gWidgets2::gtable(format(databaseData, justify = "centre"),container=databaseFrame)
    gWidgets2::size(dataTab) <- list(height=110,width=240,column.widths=c(120,50))

    #Mutations
    if(mutModel=="Stepwise"){
      mutData <- data.frame(Parameter=c("Model","Range"),Value=c(mutModel,mutPar$mutRange))
    } else mutData <-  data.frame(Parameter="Model",Value=mutModel)
    mutData2 <- data.frame(Parameter=c("Female mut. rate","Male mut. rate"),Value=c(mutPar$fMutRate,mutPar$mMutRate))
    mutData[,1] <- as.character(mutData[,1])
    mutData2[,1] <- as.character(mutData2[,1])
    mutData[,2] <- as.character(mutData[,2])
    mutData2[,2] <- as.character(mutData2[,2])
    mutData <- rbind(mutData,mutData2)

    mutGroup <- gWidgets2::ggroup(horizontal=FALSE,container=paramGroup)
    mutFrame <- gWidgets2::gframe("Mutations",container=mutGroup,horizontal=FALSE)
    #mutTab <- gWidgets2::gtable(mutData,container=mutFrame)
    mutTab <- gWidgets2::gtable(format(mutData, justify = "centre"),container=mutFrame)
    gWidgets2::size(mutTab) <- list(height=110,width=240,column.widths=c(120,50))

    contGroup <- gWidgets2::ggroup(horizontal=FALSE,container=paramGroup)
    contFrame <- gWidgets2::gframe("Contributors",container=contGroup,horizontal=FALSE)
    if(length(idxC1)<length(idxC2)){
      contData <- cbind("Pedigree 1"=idxC1[seq(idxC2)],"Pedigree 2"=idxC2)
    } else{
      contData <- cbind("Pedigree 1"=idxC1,"Pedigree 2"=idxC2[seq(idxC1)])
    }
    contData[is.na(contData)] <- ""
    contTab <- gWidgets2::gtable(format(as.data.frame(contData), justify = "centre"),container=contFrame)
    gWidgets2::size(contTab) <- list(height=90,width=240,column.widths=c(85,85))

    #Dropout/drop-in
    dropData <- data.frame(Parameter=c(paste("Dropout",names(D)),"Drop-in"),Value=c(unlist(D),di))
    dropData[,1] <- as.character(dropData[,1])
    dropData[,2] <- as.character(dropData[,2])
    colnames(dropData) <- c("Parameter","Value")
    dropGroup <-  gWidgets2::ggroup(horizontal=FALSE,container=paramGroup,spacing=7)
    dropFrame <- gWidgets2::gframe("Dropout and drop-in",container=dropGroup,horizontal=FALSE)
    dropTab <- gWidgets2::gtable(format(dropData,justify="centre"),container=dropFrame)
    gWidgets2::size(dropTab) <- list(height=110,width=240,column.widths=c(120,50))


    #### right col ####
    LRgroup2 <- gWidgets2::ggroup(container=LRwindowGroup,horizontal=FALSE,spacing=7)
    #Plot pedigrees
    pedGroup <- gWidgets2::ggroup(horizontal=TRUE,container=LRgroup2,spacing=5)
    #pedFrame <- gWidgets2::gframe("Pedigrees",container=pedGroup,horizontal=TRUE,expand=FALSE)
    pedFrame1 <- gWidgets2::gframe("Pedigree 1",container=pedGroup,horizontal=TRUE,expand=TRUE,fill=TRUE)
    #gWidgets2::size(pedFrame1) <- list(height=110,width=240)
    pedFrame2 <- gWidgets2::gframe("Pedigree 2",container=pedGroup,horizontal=TRUE,expand=TRUE,fill=TRUE)
    #gWidgets2::size(pedFrame2) <- list(height=110,width=240)
    if(all(ped1$findex==0 & ped1$mindex==0)) {
      gWidgets2::glabel("No relations",container=pedFrame1)
    }else{
      img1 <- tkrplot::tkrplot(gWidgets2::getToolkitWidget(pedFrame1),
                      fun = function() {
                        #err <- try(plot(ped1),silent=TRUE)
                        plot(ped1,cex=0.8)},hscale=0.4,vscale=0.8)
      gWidgets2::add(pedGroup,img1)
    }
    if(all(ped2$findex==0 & ped2$mindex==0)) {
      gWidgets2::glabel("No relations",container=pedFrame2)
    }else{
      img2 <- tkrplot::tkrplot(gWidgets2::getToolkitWidget(pedFrame2),
                    fun = function() {
                      #err <- try(plot(ped2),silent=TRUE)
                    plot(ped2,cex=0.8)},hscale=0.4,vscale=0.8)
    gWidgets2::add(pedGroup,img2)
    }
    #List of LRs
    print(format(Data[,1:2],digits=4,scientific=TRUE))
    gWidgets2::glabel(paste("Total LR:",format(prod(LRmarker),digits=4,scientific=TRUE)),container=LRgroup2)
    tab <- gWidgets2::gtable(format(Data[,1:2],digits=4,scientific=TRUE),container=LRgroup2, expand = TRUE)
    gWidgets2::size(tab) <- list(height = 300, width=200, column.widths = c(100, 100))


    #Prepare data to write to file
    allParam <- rbind(databaseData,mutData,dropData)
    colnames(contData) <- paste("Contributors",colnames(contData))
    allParam <- cbind(allParam," "=rep("",nrow(allParam)),rbind(contData,matrix("",nrow=nrow(allParam)-nrow(contData),ncol=ncol(contData))))
    #exportData <- data.frame(c(allParam[,1],"","Marker",Data[,1],"Total"),c(allParam[,2],"","LR",Data[,2],prod(LRmarker)),
    #                         c(rep("",nrow(allParam)+1),"Lik H1",Data[,3],prod(Data[,3])),c(rep("",nrow(allParam)+1),"Lik H2",Data[,4],prod(Data[,4])))
    #Parameters and LR
    #paramLR <- data.frame(c(allParam[,1],"","Marker",Data[,1],"Total"),c(allParam[,2],"","LR",Data[,2],prod(LRmarker)),
    #                         c(rep("",nrow(allParam)+1),"Lik H1",Data[,3],prod(Data[,3])),c(rep("",nrow(allParam)+1),"Lik H2",Data[,4],prod(Data[,4])))

    #Reference profiles
    gtlist <- split(G[,3:4],G[,1])
    GT <- do.call(cbind,gtlist)
    #Mixture profile
    mix <- E[-1]
    colnames(mix) <- c("Marker",paste("Mix.",colnames(mix[,-1]),sep=""))
    GT <- GT[rep(1:nrow(GT),length(unique(E$SampleName))),]
    rownames(GT) <- rownames(E)
    prof <- cbind(E,GT)
    prof[is.na(prof)] <- ""

    #LRgt <- cbind(rbind(Data,c("Total",prod(LRmarker),prod(Data[,3]),prod(Data[,4]))),rbind(prof[,-1],rep("",ncol(prof)-1)))
    DataTot <- rbind(Data,c("Total",prod(LRmarker),prod(Data[,3]),prod(Data[,4])))
    fill <- matrix("",nrow=nrow(DataTot),ncol=ncol(prof)-ncol(DataTot))
    DataTot <- cbind(DataTot,fill)
    colNames <- c(colnames(Data),rep(" ",ncol(fill)))
    colNamesFill <- rbind(rep("",length(colNames)),colNames)
    colnames(colNamesFill) <- colnames(DataTot)
    DataTot <- rbind(colNamesFill,DataTot)

    colnames(DataTot) <- colnames(prof)
    LRgt <- rbind(prof,DataTot)
    #LRgt <- rbind(rbind(Data,c("Total",prod(LRmarker),prod(Data[,3]),prod(Data[,4]))),rbind(prof[,-1],rep("",ncol(prof)-1))))
    LRgt <- rbind(colnames(LRgt),LRgt)
    fill <- matrix("",nrow=nrow(allParam),ncol=ncol(LRgt)-ncol(allParam))
    allParam2 <- rbind(cbind(allParam,fill),rep("",ncol(LRgt)))
    colnames(LRgt) <- colnames(allParam2)
    exportData <- rbind(allParam2,LRgt)
    exportData[is.na(exportData)] <- ""
    rownames(exportData) <- NULL
    colnames(exportData) <- c("Parameter","Value","",colnames(contData),rep("",ncol(exportData)-5))

    LRbut <- gWidgets2::gbutton("Save results to file", container=LRgroup2,handler=function(h,...){
      f_export(exportData)
      gWidgets2::dispose(h$obj)
    })

    gWidgets2::visible(LRwindow) <- TRUE
    }

  ############################################################

  mmTK = new.env() #create new environment object

  #Assign some default parameters
  defaultList <- list(Mut='Equal',
                      mutPar=list(mutRange=0.5,
                                  fMutRate=0,
                                  mMutRate=0),
                      optPar=list(theta=0,
                                  silent=0,
                                  MAF=0.001),
                      dropliste=list(d=0.00,
                                     di=0.00)
  )

  for(i in 1:length(defaultList)){
    assign(names(defaultList)[i],defaultList[[i]],envir=mmTK)
  }


  win <- gWidgets2::gwindow(paste0("RelMix v",packageVersion("relMix")),height=500,width=500,visible=FALSE)
  top <- gWidgets2::ggroup(horizontal = FALSE, container = win)
  group1 <- gWidgets2::ggroup(horizontal = TRUE, container=top)

  ###### Import data #####
  #Buttons for importing files
  dataGroup <- gWidgets2::ggroup(horizontal = FALSE, container=group1,spacing=7)
  impFrame <- gWidgets2::gframe("Import data",container=dataGroup,horizontal = FALSE)
  objMix <- gWidgets2::gbutton(text="Import mixture profile",container=impFrame,handler=f_importprof,action="mixture")
  objRef <- gWidgets2::gbutton(text="Import reference profiles",container=impFrame,handler=f_importprof,action="reference")
  objDB <- gWidgets2::gbutton(text="Database",container=impFrame,handler=f_database)

  ###### Mutations #####
  #Mutation frame
  mutGroup <- gWidgets2::ggroup(horizontal = FALSE, container=group1,spacing=7)
  mutFrame <- gWidgets2::gframe("Mutations",container=mutGroup,horizontal=FALSE)
  mutButton <- gWidgets2::gbutton("Mutations", handler=f_mutations,container=mutFrame)

  group2 <- gWidgets2::ggroup(horizontal = TRUE, container=top)

  ####### Pedigrees ######
  #Pedigree frame
  pedGroup <- gWidgets2::ggroup(horizontal = FALSE, container=group2,spacing=7)
  pedFrame <- gWidgets2::gframe("Pedigrees",container=pedGroup)
  #Should we have the option to include different contributors under each hypothesis???
  #Pedigree 1 button
  gWidgets2::glabel("Pedigree 1",container=pedFrame)
  objPed <- gcombobox(c("Paternity","Unrelated","Custom"), selected=0,container=pedFrame,handler=f_pedigree,action="1")
  #Pedigree 2 button
  gWidgets2::glabel("Pedigree 2",container=pedFrame)
  objPed2 <- gcombobox(c("Paternity","Unrelated","Custom"), selected=0,container=pedFrame,handler=f_pedigree,action="2")

  ######## Contributors #######
  contFrame <- gWidgets2::gframe("Contributors",container=group2,horizontal=TRUE)
  objCont <- gWidgets2::gbutton(text="Specify contributors",container=contFrame,handler=f_contributors)

  ####### Dropout/drop-in #########
  group3 <- gWidgets2::ggroup(horizontal = TRUE, container=top,spacing=10)
  #Dropout button
  dropFrame <- gWidgets2::gframe("Dropout and drop-in",container=group3,horizontal=TRUE)
  dropButton <- gWidgets2::gbutton("Specify dropout and drop-in", handler=f_dropout, container=dropFrame)

  ###### Compute LR ######
  group4 <- gWidgets2::ggroup(horizontal = TRUE, container=top,spacing=10)
  #Compute LR button
  LRbutton <- gWidgets2::gbutton("Compute LR", container=group4, handler=function(h,...){
    f_LR()
  })

  ####### Restart ######
  gWidgets2::addSpace(group4, 10)
  restartButton = gWidgets2::gbutton(text="RESTART",container=group4,handler=function(h,...) {
    gWidgets2::dispose(win) #remove main window
    relMixGUI() #reopen relMix
  })

  gWidgets2::visible(win) <- TRUE
}


