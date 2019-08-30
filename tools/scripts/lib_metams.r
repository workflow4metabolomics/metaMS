# lib_metams.r version 2.1.1
# R function for metaMS runGC under W4M
# author Yann GUITTON CNRS IRISA/LINA Idealg project 2014-2015
# author Yann GUITTON Oniris Laberca 2015-2017


#@author G. Le Corguille
# This function will
# - load the packages
# - display the sessionInfo
loadAndDisplayPackages <- function(pkgs) {
    for(pkg in pkgs) suppressPackageStartupMessages( stopifnot( library(pkg, quietly=TRUE, logical.return=TRUE, character.only=TRUE)))

    sessioninfo = sessionInfo()
    cat(sessioninfo$R.version$version.string,"\n")
    cat("Main packages:\n")
    for (pkg in names(sessioninfo$otherPkgs)) { cat(paste(pkg,packageVersion(pkg)),"\t") }; cat("\n")
    cat("Other loaded packages:\n")
    for (pkg in names(sessioninfo$loadedOnly)) { cat(paste(pkg,packageVersion(pkg)),"\t") }; cat("\n")
}

#This function list the compatible files within the directory as xcms did
#@author Gildas Le Corguille lecorguille@sb-roscoff.fr ABiMS TEAM
getMSFiles <- function (directory) {
    filepattern <- c("[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]","[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]")
    filepattern <- paste(paste("\\.", filepattern, "$", sep=""),collapse="|")
    info <- file.info(directory)
    listed <- list.files(directory[info$isdir], pattern=filepattern,recursive=TRUE, full.names=TRUE)
    files <- c(directory[!info$isdir], listed)
    exists <- file.exists(files)
    files <- files[exists]
    return(files)
}

# This function retrieve a xset like object
#@author Gildas Le Corguille lecorguille@sb-roscoff.fr
getxcmsSetObject <- function(xobject) {
    # XCMS 1.x
    if (class(xobject) == "xcmsSet")
        return (xobject)
    # XCMS 3.x
    if (class(xobject) == "XCMSnExp") {
        # Get the legacy xcmsSet object
        suppressWarnings(xset <- as(xobject, 'xcmsSet'))
        if (!is.null(xset@phenoData$sample_group))
            sampclass(xset) <- xset@phenoData$sample_group
        else
            sampclass(xset) <- "."
        return (xset)
    }
}

#J.Saint-Vanne
#Function to correct the file names which can be written like "..alg8.mzData" and we just want "alg8"
getCorrectFileName <- function(peaktable,sampleMetadata){

    #Try to start for the first sample, avoid description of line with colnamesdontwant
    i <- 1
    while(!(sampleMetadata[1,1] %in% strsplit(colnames(peaktable)[i],"\\.")[[1]])) {
        if(i < length(peaktable)) {
            i <- i + 1
        } else {
            break
        }
    }
    #i now correspond to the first column with a sample
    for(j in 1:(nrow(sampleMetadata))) {
        col <- j + i-1 #minus 1 cause i is the good column to start and j start at 1
        if(col <= length(colnames(peaktable))) {
            newname <- gsub("(^.*)(\\..*$)","\\1",colnames(peaktable)[col])
            if(newname != sampleMetadata[j,1]){
                #Correction for 2 points starting the name (I don't know why they are here...)
                if(".." == gsub("(^\\.+)(.*)","\\1",newname)){
                    newname <- sub("(^\\.+)(.*)","\\2",newname)
                }
            }
            colnames(peaktable)[col] <- newname
        }
    }
    return(peaktable)
}


#@author Gildas Le Corguille lecorguille@sb-roscoff.fr
# This function get the raw file path from the arguments
getRawfilePathFromArguments <- function(singlefile, zipfile, listArguments) {
    if (!is.null(listArguments[["zipfile"]]))           zipfile = listArguments[["zipfile"]]
    if (!is.null(listArguments[["zipfilePositive"]]))   zipfile = listArguments[["zipfilePositive"]]
    if (!is.null(listArguments[["zipfileNegative"]]))   zipfile = listArguments[["zipfileNegative"]]

    if (!is.null(listArguments[["singlefile_galaxyPath"]])) {
        singlefile_galaxyPaths = listArguments[["singlefile_galaxyPath"]];
        singlefile_sampleNames = listArguments[["singlefile_sampleName"]]
    }
    if (!is.null(listArguments[["singlefile_galaxyPathPositive"]])) {
        singlefile_galaxyPaths = listArguments[["singlefile_galaxyPathPositive"]];
        singlefile_sampleNames = listArguments[["singlefile_sampleNamePositive"]]
    }
    if (!is.null(listArguments[["singlefile_galaxyPathNegative"]])) {
        singlefile_galaxyPaths = listArguments[["singlefile_galaxyPathNegative"]];
        singlefile_sampleNames = listArguments[["singlefile_sampleNameNegative"]]
    }
    if (exists("singlefile_galaxyPaths")){
        singlefile_galaxyPaths = unlist(strsplit(singlefile_galaxyPaths,","))
        singlefile_sampleNames = unlist(strsplit(singlefile_sampleNames,","))

        singlefile=NULL
        for (singlefile_galaxyPath_i in seq(1:length(singlefile_galaxyPaths))) {
            singlefile_galaxyPath=singlefile_galaxyPaths[singlefile_galaxyPath_i]
            singlefile_sampleName=singlefile_sampleNames[singlefile_galaxyPath_i]
            singlefile[[singlefile_sampleName]] = singlefile_galaxyPath
        }
    }
    for (argument in c("zipfile","zipfilePositive","zipfileNegative","singlefile_galaxyPath","singlefile_sampleName","singlefile_galaxyPathPositive","singlefile_sampleNamePositive","singlefile_galaxyPathNegative","singlefile_sampleNameNegative")) {
        listArguments[[argument]]=NULL
    }
    return(list(zipfile=zipfile, singlefile=singlefile, listArguments=listArguments))
}


#@author Gildas Le Corguille lecorguille@sb-roscoff.fr
# This function retrieve the raw file in the working directory
#   - if zipfile: unzip the file with its directory tree
#   - if singlefiles: set symlink with the good filename
retrieveRawfileInTheWorkingDirectory <- function(singlefile, zipfile) {
    if(!is.null(singlefile) && (length("singlefile")>0)) {
        for (singlefile_sampleName in names(singlefile)) {
            singlefile_galaxyPath = singlefile[[singlefile_sampleName]]
            if(!file.exists(singlefile_galaxyPath)){
                error_message=paste("Cannot access the sample:",singlefile_sampleName,"located:",singlefile_galaxyPath,". Please, contact your administrator ... if you have one!")
                print(error_message); stop(error_message)
            }
            file.symlink(singlefile_galaxyPath,singlefile_sampleName)
        }
        directory = "."
    }
    if(!is.null(zipfile) && (zipfile!="")) {
        if(!file.exists(zipfile)){
            error_message=paste("Cannot access the Zip file:",zipfile,". Please, contact your administrator ... if you have one!")
            print(error_message)
            stop(error_message)
        }

        #list all file in the zip file
        #zip_files=unzip(zipfile,list=T)[,"Name"]

        #unzip
        suppressWarnings(unzip(zipfile, unzip="unzip"))

        #get the directory name
        filesInZip=unzip(zipfile, list=T);
        directories=unique(unlist(lapply(strsplit(filesInZip$Name,"/"), function(x) x[1])));
        directories=directories[!(directories %in% c("__MACOSX")) & file.info(directories)$isdir]
        directory = "."
        if (length(directories) == 1) directory = directories

        cat("files_root_directory\t",directory,"\n")
    }
    return (directory)
}

##ADDITIONS FROM Y. Guitton
getBPC <- function(file,rtcor=NULL, ...) {
    #Error with chromatogram function with CDF file (only on Galaxy, not in R command...)
    if(MSnbase:::isCdfFile(file)){
        xdata <- readMSData(file,mode="onDisk")
        chromBPC <- bpi(xdata,initial=FALSE)
        chromRT <- rtime(xdata)
        cbind(chromRT,chromBPC)
    } else {
        xdata <- readMSData(file,mode="onDisk")
        c <- chromatogram(xdata, aggregationFun = "max")
        cbind(c@.Data[[1]]@rtime,c@.Data[[1]]@intensity)
    }
}

getBPC2s <- function (files, xset = NULL, pdfname="BPCs.pdf", rt = c("raw","corrected"), scanrange=NULL) {
    require(xcms)
  
    #create sampleMetadata, get sampleMetadata and class
    if(!is.null(xset)) {
    	#When files come from XCMS3 directly before metaMS
    	sampleMetadata <- xset@phenoData
    } else {
    	#When files come from a zip directory with raw files before metaMS
	    sampleMetadata<-xcms:::phenoDataFromPaths(files)
    }
    class<-unique(sampleMetadata[,"class"]) #create phenoData like table
    classnames<-vector("list",length(class))
    for (i in 1:length(class)){
        classnames[[i]]<-which( sampleMetadata[,"class"]==class[i])
    }

    N <- dim(sampleMetadata)[1]
    TIC <- vector("list",N)

    for (j in 1:N) {
        cat(files[j],"\n")
        TIC[[j]] <- getBPC(files[j])
        #good for raw 
        # seems strange for corrected
        #errors if scanrange used in xcmsSetgeneration
        if (!is.null(xcmsSet) && rt == "corrected"){
            rtcor <- xcmsSet@rt$corrected[[j]] 
        }else{
            rtcor <- NULL
        }
        TIC[[j]] <- getBPC(files[j],rtcor=rtcor)
    }

    pdf(pdfname,w=16,h=10)
    cols <- rainbow(N)
    lty = 1:N
    pch = 1:N
    #search for max x and max y in BPCs
    xlim = range(sapply(TIC, function(x) range(x[,1],na.rm=TRUE)),na.rm=TRUE)
    ylim = range(sapply(TIC, function(x) range(x[,2],na.rm=TRUE)),na.rm=TRUE)
    ylim = c(-ylim[2], ylim[2])

    ##plot start
    if (length(class)>2){
        for (k in 1:(length(class)-1)){
            for (l in (k+1):length(class)){
                cat(paste(class[k],"vs",class[l],sep=" ","\n")) 
                plot(0, 0, type="n", xlim = xlim/60, ylim = ylim, main = paste("Base Peak Chromatograms \n","BPCs_",class[k]," vs ",class[l], sep=""), xlab = "Retention Time (min)", ylab = "BPC")
                colvect<-NULL
                for (j in 1:length(classnames[[k]])) {
                    tic <- TIC[[classnames[[k]][j]]]
                    # points(tic[,1]/60, tic[,2], col = cols[i], pch = pch[i], type="l")
                    points(tic[,1]/60, tic[,2], col = cols[classnames[[k]][j]], pch = pch[classnames[[k]][j]], type="l")
                    colvect<-append(colvect,cols[classnames[[k]][j]])
                }
                for (j in 1:length(classnames[[l]])) {
                    # i=class2names[j]
                    tic <- TIC[[classnames[[l]][j]]]
                    points(tic[,1]/60, -tic[,2], col = cols[classnames[[l]][j]], pch = pch[classnames[[l]][j]], type="l")
                    colvect<-append(colvect,cols[classnames[[l]][j]])
                }
                legend("topright",paste(gsub("(^.+)\\..*$","\\1",basename(files[c(classnames[[k]],classnames[[l]])]))), col = colvect, lty = lty, pch = pch)
            }
        }
    }#end if length >2

    if (length(class)==2){
        k=1
		l=2
        colvect<-NULL
        plot(0, 0, type="n", xlim = xlim/60, ylim = ylim, main = paste("Base Peak Chromatograms \n","BPCs_",class[k],"vs",class[l], sep=""), xlab = "Retention Time (min)", ylab = "BPC")

        for (j in 1:length(classnames[[k]])) {
            tic <- TIC[[classnames[[k]][j]]]
            # points(tic[,1]/60, tic[,2], col = cols[i], pch = pch[i], type="l")
            points(tic[,1]/60, tic[,2], col = cols[classnames[[k]][j]], pch = pch[classnames[[k]][j]], type="l")
            colvect<-append(colvect,cols[classnames[[k]][j]])
        }
        for (j in 1:length(classnames[[l]])) {
            # i=class2names[j]
            tic <- TIC[[classnames[[l]][j]]]
            points(tic[,1]/60, -tic[,2], col = cols[classnames[[l]][j]], pch = pch[classnames[[l]][j]], type="l")
            colvect<-append(colvect,cols[classnames[[l]][j]])
        }
        legend("topright",paste(gsub("(^.+)\\..*$","\\1",basename(files[c(classnames[[k]],classnames[[l]])]))), col = colvect, lty = lty, pch = pch)
    }#end length ==2
    
    if (length(class)==1){
        k=1
		ylim = range(sapply(TIC, function(x) range(x[,2],na.rm=TRUE)),na.rm=TRUE)
        colvect<-NULL
        plot(0, 0, type="n", xlim = xlim/60, ylim = ylim, main = paste("Base Peak Chromatograms \n","BPCs_",class[k], sep=""), xlab = "Retention Time (min)", ylab = "BPC")

        for (j in 1:length(classnames[[k]])) {
            tic <- TIC[[classnames[[k]][j]]]
            # points(tic[,1]/60, tic[,2], col = cols[i], pch = pch[i], type="l")
            points(tic[,1]/60, tic[,2], col = cols[classnames[[k]][j]], pch = pch[classnames[[k]][j]], type="l")
            colvect<-append(colvect,cols[classnames[[k]][j]])
        }
        legend("topright",paste(gsub("(^.+)\\..*$","\\1",basename(files[c(classnames[[k]])]))), col = colvect, lty = lty, pch = pch)
    }#end length ==1
    dev.off()
}

getTIC <- function(file,rtcor=NULL) {
    #Error with chromatogram function with CDF file (only on Galaxy, not in R command...)
    if(MSnbase:::isCdfFile(file)){
        xdata <- readMSData(file,mode="onDisk")
        chromTIC <- tic(xdata,initial=FALSE)
        chromRT <- rtime(xdata)
        cbind(chromRT,chromTIC)
    } else {
        xdata <- readMSData(file,mode="onDisk")
        c <- chromatogram(xdata, aggregationFun = "sum")
        cbind(c@.Data[[1]]@rtime,c@.Data[[1]]@intensity)
    }
}

##  overlay TIC from all files in current folder or from xcmsSet, create pdf
getTIC2s <- function(files, xset=NULL, pdfname="TICs.pdf", rt=c("raw","corrected")) {
    require(xcms)

    #create sampleMetadata, get sampleMetadata and class
    if(!is.null(xset)){
		#When files come from XCMS3 before metaMS treatment
    	sampleMetadata<-xset@phenoData
    } else {
		#When files come from a zip directory with raw files before metaMS
    	sampleMetadata<-xcms:::phenoDataFromPaths(files)
    }
    class<-as.vector(levels(sampleMetadata[,"class"])) #create phenoData like table
    classnames<-vector("list",length(class))
    for (i in 1:length(class)){
        classnames[[i]]<-which( sampleMetadata[,"class"]==class[i])
    }
        
    N <- dim(sampleMetadata)[1]
    TIC <- vector("list",N)

    for (i in 1:N) {
        cat(files[i],"\n")
        if (!is.null(xcmsSet) && rt == "corrected")
            rtcor <- xcmsSet@rt$corrected[[i]]
        else
            rtcor <- NULL
        TIC[[i]] <- getTIC(files[i],rtcor=rtcor)
    }
 
    pdf(pdfname,w=16,h=10)
    cols <- rainbow(N)
    lty = 1:N
    pch = 1:N
    #search for max x and max y in TICs
    xlim = range(sapply(TIC, function(x) range(x[,1],na.rm=TRUE)),na.rm=TRUE)
    ylim = range(sapply(TIC, function(x) range(x[,2],na.rm=TRUE)),na.rm=TRUE)
    ylim = c(-ylim[2], ylim[2])  
	  
    ##plot start
    if (length(class)>2){
        for (k in 1:(length(class)-1)){
            for (l in (k+1):length(class)){
                print(paste(class[k],"vs",class[l],sep=" ")) 
                plot(0, 0, type="n", xlim = xlim/60, ylim = ylim, main = paste("Total Ion Chromatograms \n","TICs_",class[k]," vs ",class[l], sep=""), xlab = "Retention Time (min)", ylab = "TIC")
                colvect<-NULL
                for (j in 1:length(classnames[[k]])) {
                    tic <- TIC[[classnames[[k]][j]]]
                    # points(tic[,1]/60, tic[,2], col = cols[i], pch = pch[i], type="l")
                    points(tic[,1]/60, tic[,2], col = cols[classnames[[k]][j]], pch = pch[classnames[[k]][j]], type="l")
                    colvect<-append(colvect,cols[classnames[[k]][j]])
                }
                for (j in 1:length(classnames[[l]])) {
                    # i=class2names[j]
                    tic <- TIC[[classnames[[l]][j]]]
                    points(tic[,1]/60, -tic[,2], col = cols[classnames[[l]][j]], pch = pch[classnames[[l]][j]], type="l")
                    colvect<-append(colvect,cols[classnames[[l]][j]])
                }
                legend("topright",paste(gsub("(^.+)\\..*$","\\1",basename(files[c(classnames[[k]],classnames[[l]])]))), col = colvect, lty = lty, pch = pch)
            }
        }
    }#end if length >2

    if (length(class)==2){
        k=1
        l=2
        plot(0, 0, type="n", xlim = xlim/60, ylim = ylim, main = paste("Total Ion Chromatograms \n","TICs_",class[k],"vs",class[l], sep=""), xlab = "Retention Time (min)", ylab = "TIC")
        colvect<-NULL
        for (j in 1:length(classnames[[k]])) {
            tic <- TIC[[classnames[[k]][j]]]
            # points(tic[,1]/60, tic[,2], col = cols[i], pch = pch[i], type="l")
            points(tic[,1]/60, tic[,2], col = cols[classnames[[k]][j]], pch = pch[classnames[[k]][j]], type="l")
            colvect<-append(colvect,cols[classnames[[k]][j]])
        }
        for (j in 1:length(classnames[[l]])) {
            # i=class2names[j]
            tic <- TIC[[classnames[[l]][j]]]
            points(tic[,1]/60, -tic[,2], col = cols[classnames[[l]][j]], pch = pch[classnames[[l]][j]], type="l")
            colvect<-append(colvect,cols[classnames[[l]][j]])
        }
        legend("topright",paste(gsub("(^.+)\\..*$","\\1",basename(files[c(classnames[[k]],classnames[[l]])]))), col = colvect, lty = lty, pch = pch)
    }#end length ==2

    if (length(class)==1){
        k=1
        ylim = range(sapply(TIC, function(x) range(x[,2],na.rm=TRUE)),na.rm=TRUE)
        plot(0, 0, type="n", xlim = xlim/60, ylim = ylim, main = paste("Total Ion Chromatograms \n","TICs_",class[k], sep=""), xlab = "Retention Time (min)", ylab = "TIC")
        colvect<-NULL
        for (j in 1:length(classnames[[k]])) {
            tic <- TIC[[classnames[[k]][j]]]
            # points(tic[,1]/60, tic[,2], col = cols[i], pch = pch[i], type="l")
            points(tic[,1]/60, tic[,2], col = cols[classnames[[k]][j]], pch = pch[classnames[[k]][j]], type="l")
            colvect<-append(colvect,cols[classnames[[k]][j]])
        }
        legend("topright",paste(gsub("(^.+)\\..*$","\\1",basename(files[c(classnames[[k]])]))), col = colvect, lty = lty, pch = pch)
    }#end length ==1
    dev.off()
}


#Update by J.Saint-Vanne
##addition for quality control of peak picking
#metaMS EIC and pspectra plotting option
#version 20190520
#only for Galaxy 
plotUnknowns<-function(resGC, unkn="", DB=NULL, fileFrom=NULL){

    ##Annotation table each value is a pcgrp associated to the unknown 
    ##NOTE pcgrp index are different between xcmsSet and resGC due to filtering steps in metaMS
    ##R. Wehrens give me some clues on that and we found a correction
    
	#correction of annotation matrix due to pcgrp removal by quality check in runGCresult
	#matrix of correspondance between an@pspectra and filtered pspectra from runGC
    #Select only pspectra which correpond to them select in metaMS
    # col1 = filtered spectra from runGC and col2 = an@spectra
	allPCGRPs <-lapply(1:length(resGC$xset),
				function(i) {
					an <- resGC$xset[[i]]
					huhn <- an@pspectra[which(sapply(an@pspectra, length) >=
					metaSetting(resGC$settings,
					"DBconstruction.minfeat"))]
					matCORR<-cbind(1:length(huhn), match(huhn, an@pspectra))
				})
    #Build a new annotation list with sampnames and pseudospectra number from xset
    helpannotation <- list()
    for(j in 1:length(resGC$xset)){
        helpannotation[[j]] <- resGC$annotation[[j]][1:2]
        pspvector <- vector()
        
        if(nrow(helpannotation[[j]]) > 0){
            for(i in 1: nrow(helpannotation[[j]])){
                #Find corresponding pspec
                psplink <- allPCGRPs[[j]][match(helpannotation[[j]][i,1],allPCGRPs[[j]]),2]
                pspvector <- c(pspvector,psplink)
                #Change the annotation column into sampname column
                if(helpannotation[[j]][i,2] < 0){
                    #It's an unknown
                    new_name <- paste("Unknown",abs(as.integer(helpannotation[[j]][i,2])))
                    helpannotation[[j]][i,2] <- new_name
                }else{
                    #It has been found in local database
                    for(k in 1:length(DB)){
                        if(helpannotation[[j]][i,2] == k){
                            helpannotation[[j]][i,2] <- DB[[k]]$Name
                            break
                        }
                    }
                }
            }
        }
        helpannotation[[j]] <- cbind(helpannotation[[j]],pspvector)
        names(helpannotation)[j] <- names(resGC$annotation[j])
    }
    peaktable <- resGC$PeakTable
		
	par (mar=c(5, 4, 4, 2) + 0.1)
	#For each unknown
	for (l in 1:length(unkn)){
		#recordPlot
		perpage=3 #if change change layout also!
		dev.new(width=21/2.54, height=29.7/2.54, file=paste("Unknown_",unkn[l],".pdf", sep="")) #A4 pdf
		# par(mfrow=c(perpage,2))
		layout(matrix(c(1,1,2,3,4,4,5,6,7,7,8,9), 6, 2, byrow = TRUE), widths=rep(c(1,1),perpage), heights=rep(c(1,5),perpage))
		# layout.show(6)
		oma.saved <- par("oma")
		par(oma = rep.int(0, 4))
		par(oma = oma.saved)
		o.par <- par(mar = rep.int(0, 4))
		on.exit(par(o.par))
		#For each sample
		for (c in 1:length(resGC$xset)) {
            #get sample name
            sampname<-basename(resGC$xset[[c]]@xcmsSet@filepaths)
            #remove .cdf, .mzXML filepattern
            filepattern <- c("[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]", 
                                    "[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]")
            filepattern <- paste(paste("\\.", filepattern, "$", sep = ""), collapse = "|")
            sampname<-gsub(filepattern, "",sampname)                     
            
            an<-resGC$xset[[c]]
            if(fileFrom == "zipfile") {
                an@xcmsSet@filepaths <- paste0("./",an@xcmsSet@phenoData[,"class"],"/",basename(an@xcmsSet@filepaths))
            }#else {
                #print(an@xcmsSet@filepaths)
                #an@xcmsSet@filepaths <- paste0("./",basename(an@xcmsSet@filepaths))
            #}
            #Make different plot according to type of file (CDF file)
            if(MSnbase:::isCdfFile(resGC$xset[[c]]@xcmsSet@filepaths)) {
                plotEICsCDF(an, helpannotation, filepattern, peaktable[unkn[l],1], sampname)
            }else{
                plotEICsOthers(an, helpannotation, filepattern, peaktable[unkn[l],1], sampname)
            }
        }
        graphics.off()
    }#end  for unkn[l]
}#end function

plotEICsOthers <- function(an, helpannotation, filepattern, name, sampname){
    title1<-paste(name,"from",sampname, sep = " ")
    #Find the good annotation for this sample
    for(a in 1:length(helpannotation)){
        if(gsub(filepattern, "", names(helpannotation)[a]) == paste0("./",sampname)){
            #Find the unkn or the matched std in this sample
            findunkn <- FALSE
            if(nrow(helpannotation[[a]]) > 0){
                for(r in 1:nrow(helpannotation[[a]])){
                    #Recherche du bon inconnu
                    if(helpannotation[[a]][r,"annotation"] == name){
                        findunkn <- TRUE
                        pcgrp <- helpannotation[[a]][r,"pspvector"]
                        par (mar=c(0, 0, 0, 0) + 0.1)
                        #Write title
                        plot.new()
                        box()
                        text(0.5, 0.5, title1, cex=2)                       
                        par (mar=c(3, 2.5, 3, 1.5) + 0.1)
                        #Window for EIC
                        plotEICs(an, pspec=pcgrp, maxlabel=2)
                        #Window for pseudospectra
                        plotPsSpectrum(an, pspec=pcgrp, maxlabel=2)
                    }
                }
            }
            
            if(!findunkn) {
                par (mar=c(0, 0, 0, 0) + 0.1)
                #Write title
                plot.new()
                box()
                text(0.5, 0.5, title1, cex=2)
                #Window for EIC
                plot.new()
                box()
                text(0.5, 0.5, "NOT FOUND", cex=2)
                #Window for pseudospectra
                plot.new()
                box()
                text(0.5, 0.5, "NOT FOUND", cex=2)
            }
            break
        }
    } 
}
			

#Inspire from plotEICs function "https://github.com/sneumann/CAMERA/blob/master/R/xsVisualise.R"
#Modify by J.Saint-Vanne to be able to run it on CDF files
plotEICsCDF <- function(an, helpannotation, filepattern, name, sampname){
    title1<-paste(name,"from",sampname, sep = " ")
    #Find the good annotation for this sample
    for(a in 1:length(helpannotation)){
        if(gsub(filepattern, "", names(helpannotation)[a]) == paste0("./",sampname)){
            #Find the unkn or the matched std in this sample
            findunkn <- FALSE
            if(nrow(helpannotation[[a]]) > 0){
                for(r in 1:nrow(helpannotation[[a]])){
                    #Recherche du bon inconnu
                    if(helpannotation[[a]][r,"annotation"] == name){
                        findunkn <- TRUE
                        pcgrp <- helpannotation[[a]][r,"pspvector"]
                        par (mar=c(0, 0, 0, 0) + 0.1)
                        #Write title
                        plot.new()
                        box()
                        text(0.5, 0.5, title1, cex=2)                       
                        par (mar=c(3, 2.5, 3, 1.5) + 0.1)
                  
                        #####################
                        #Remake plotEIC function for CDF
                        #####################
                        maxlabel=2
                        sleep=0
                        method="bin"
                        #Which different samples need to be addressed for extracting raw data
                        smpls <- unique(an@psSamples[pcgrp])

                        xeic  <- new("xcmsEIC");
                        xeic@rtrange <- matrix(nrow=length(pcgrp), ncol=2)
                        xeic@mzrange <- matrix(nrow=length(pcgrp), ncol=2)
                        #iterator for ps-grp
                        pcpos <- 1
                        #one second overlap
                        rtmargin <- 1

                        for (d in seq(along=smpls)) { ## sample-wise EIC collection
                            #read rawData into one xcmsRaw
                            xraw <- xcmsRaw(an@xcmsSet@filepaths[smpls[d]], profmethod=method, profstep=0)
                            pspecS <- pcgrp[which(an@psSamples[pcgrp] == smpls[d])]
                            ## getting ALL peaks from the current sample (not that bad)
                            peaks <- CAMERA:::getPeaks(an@xcmsSet, smpls[d])  
                            eic <- lapply (pspecS, function(pc) {
                                pidx <- an@pspectra[[pc]]
                                pks <- peaks[pidx, , drop=FALSE]
                                gks <- an@groupInfo[pidx, , drop=FALSE]
                                nap <- which(is.na(pks[, 1]))
                                pks[nap, ] <- cbind(gks[nap, c(1:6), drop=FALSE], matrix(nrow=length(nap), ncol=5, 0))
                                bbox <- c(rtmin = min(pks[, "rtmin"])-rtmargin,
                                          rtmax = max(pks[, "rtmax"])+rtmargin,
                                          mzmin = min(pks[, "mzmin"]),
                                          mzmax = max(pks[, "mzmax"]))
                                eic <- lapply(1:nrow(pks), function(x){
                                    reic <- xcms::rawEIC(xraw,rtrange=pks[x,c("rtmin","rtmax"),drop=FALSE],
                                                         mzrange=pks[x,c("mzmin","mzmax"),drop=FALSE])
                                    new <- cbind("rt"=xraw@scantime[reic$scan],"intensity"=reic$intensity)  
                                })
                                mzrange <- apply(pks, 1,function(x){
                                    x[c("mzmin","mzmax")]
                                })
                                #Transposition
                                mzrange <- t(mzrange)
                                rtrange <- apply(pks, 1, function(x){
                                    x[c("rtmin","rtmax")]
                                })
                                #Transposition
                                rtrange <- t(rtrange)
                                #Create EIC object
                                eic <- new("xcmsEIC", eic = list(xcmsRaw=eic), rt = "raw", 
                                           rtrange=rtrange, mzrange=mzrange, groupnames = character(0))

                                #write resulting bounding box into xcmsEIC
                                xeic@rtrange[pcpos, ] <<- bbox[c("rtmin","rtmax")]
                                xeic@mzrange[pcpos, ] <<- bbox[c("mzmin","mzmax")]

                                cat("-->", pcpos, "\n")
                                pcpos <<- pcpos+1
                                eic@eic[[1]]
                            })
                            xeic@eic <- c(xeic@eic, eic)
                        }
                        names(xeic@eic) <- paste("Pseudospectrum ", pcgrp, sep="")

                        if(maxlabel > 0){
                            col <- rainbow(maxlabel);
                        } else {
                            col <- c();
                        }
  
                        ##
                        ## Loop through all pspectra
                        ##

                        for (ps in seq(along=pcgrp)) {
                            EIC <- xeic@eic[[ps]];
                            pidx <- an@pspectra[[pcgrp[ps]]];
                            peaks <- CAMERA:::getPeaks(an@xcmsSet,an@psSamples[pcgrp[ps]])[pidx,,drop=FALSE]
                            grps <- an@groupInfo[pidx, ]
                            nap <- which(is.na(peaks[, 1]))
                            naps <- rep(FALSE, nrow(peaks))

                            if(length(nap) > 0){  
                                naps[nap] <- TRUE;
                                peaks[nap,] <- cbind(grps[nap,c(1:6), drop=FALSE], matrix(nrow=length(nap), ncol=5,0))
                            }
    
                            main <- paste("Pseudospectrum ", pcgrp[ps], sep="")
    
                            ## Calculate EICs and plot ranges
                            neics     <- length(pidx)
                            lmaxlabel <- min(maxlabel, neics)
                            eicidx <- 1:neics
                            maxint <- numeric(length(eicidx))
                            for (j in eicidx) {
                                maxint[j] <- max(EIC[[j]][, "intensity"])
                            }
                            o  <- order(maxint, decreasing = TRUE)
                            rt <- xeic@rtrange[ps, ]
                            rt.min <- round(mean(peaks[, "rtmin"]), digits=3)
                            rt.med <- round(peaks[o[1], "rt"], digits=3)
                            rt.max <- round(mean(peaks[, "rtmax"]), digits=3)
    
                            ## Open Plot
                            plot(0, 0, type = "n", xlim = rt, ylim = c(0, max(maxint)),xaxs='i',
                                 xlab = "Retention Time (seconds)", ylab = "Intensity",
                                 main = paste("Extracted Ion Chromatograms for ", main,"\nTime: From",rt.min,"to",rt.max,", mean",rt.med))
  
                            ## Plot Peak and surrounding
                            lcol <- rgb(0.6, 0.6, 0.6)
                            lcol <- c(col, rep(lcol, max(nrow(peaks) - maxlabel, 0)))
                            cnt  <- 1
                            for (j in eicidx[o]) {
                                pts <- xeic@eic[[ps]][[j]]
                                points(pts, type = "l", col = lcol[cnt]);
                                cnt <- cnt + 1;
                                peakrange <- peaks[,c("rtmin","rtmax"), drop=FALSE]
                                ptsidx <- pts[,"rt"] >= peakrange[j,1] & pts[,"rt"] <= peakrange[j,2]
                                if (naps[j]){ 
                                    points(pts[ptsidx, ], type = "l", col = col[j], lwd=1.3, lty=3)
                                } else {
                                    points(pts[ptsidx, ], type = "l", col = col[j], lwd=1.3)
                                }
                            }
                            ## Plot Annotation Legend
                            pspectrum <- getpspectra(an, grp=pcgrp[ps])
                            mz <- pspectrum[o, "mz"];
                            #Check adduct annotation
                            if (lmaxlabel > 0 & "adduct" %in% colnames(pspectrum)) {
                                adduct <- sub("^ ", "", pspectrum[o, "adduct"]) #Remove Fronting Whitespaces
                                mass   <- sapply(strsplit(adduct, " "), function(x) {x[2]})
                                adduct <- sapply(strsplit(adduct, " "), function(x) {x[1]})
                                umass <- unique(na.omit(mass[1:maxlabel]));
                                adduct[is.na(adduct)] <- "";
                                test <- vector("list",length=length(mz));
                                mz <- format(pspectrum[o[1:lmaxlabel], "mz"], digits=5);
        
                                if(length(umass) > 0){
                                    for(i in 1:length(umass)){
                                        ini <- which(mass==umass[i]);
                                        for(ii in 1:length(ini)){
                                            firstpart  <- strsplit(adduct[ini[ii]], "M")[[1]][1]
                                            secondpart <- strsplit(adduct[ini[ii]], "M")[[1]][2]
                                            masspart   <- mz[ini[ii]];
                                            test[[ini[ii]]] <- substitute(paste(masspart, " ", firstpart, M[i], secondpart), 
                                                                list(firstpart=firstpart, i=i, secondpart=secondpart,
                                                                     masspart = masspart))
                                        }
                                    }
                                }
        
                                for(i in seq(along=lmaxlabel)){
                                    if(is.null(test[[i]])){
                                        test[[i]] <- mz[i];
                                    }
                                }
        
                                leg <- as.expression(test[1:lmaxlabel]);  
                                legend("topright", legend=leg, col=lcol, lty=1)
                            }
            
                            if (sleep > 0) {
                                Sys.sleep(sleep)
                            }
                        }

                        #Window for pseudospectra
                        plotPsSpectrum(an, pspec=pcgrp, maxlabel=2)
                    }
                } 
            }
            
            if(!findunkn) {
                par (mar=c(0, 0, 0, 0) + 0.1)
                #Write title
                plot.new()
                box()
                text(0.5, 0.5, title1, cex=2)
                #Window for EIC
                plot.new()
                box()
                text(0.5, 0.5, "NOT FOUND", cex=2)
                #Window for pseudospectra
                plot.new()
                box()
                text(0.5, 0.5, "NOT FOUND", cex=2)
            }
            break
        }
    }
}