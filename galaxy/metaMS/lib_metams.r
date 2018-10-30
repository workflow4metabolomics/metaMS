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

##ADDITIONS FROM Y. Guitton
getBPC <- function(file,rtcor = NULL, ...) {
    print("getBPC")
    print(file)
    object <- xcmsRaw(file)
	sel <- profRange(object, ...)
	cbind(if (is.null(rtcor)) object@scantime[sel$scanidx] else rtcor ,xcms:::colMax(object@env$profile[sel$massidx,sel$scanidx,drop=FALSE]))    
}

getBPC2s <- function (files, pdfname = "BPCs.pdf", rt = c("raw","corrected"), scanrange = NULL) {
    print("getBPC2s")
    require(xcms)
    print(files)        
    #create sampleMetadata, get sampleMetadata and class
    sampleMetadata <- xcms:::phenoDataFromPaths(unname(files))
    class <- class <- as.vector(levels(sampleMetadata[,"class"])) #create phenoData like table
    classnames <- vector("list",length(class))
    for (i in 1:length(class)){
        classnames[[i]] <- which( sampleMetadata[,1] == class[i])
    }

    N <- dim(sampleMetadata)[1]
      
    TIC <- vector("list",N)

    for (j in 1:N) {

        cat("Processing file",names(files[j]),"\n")
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

    pdf(pdfname, w = 16, h = 10)
    cols <- rainbow(N)
    lty = 1:N
    pch = 1:N
    #search for max x and max y in BPCs
    xlim = range(sapply(TIC, function(x) range(x[,1])))
    ylim = range(sapply(TIC, function(x) range(x[,2])))
    ylim = c(-ylim[2], ylim[2])


    ##plot start
    
    if (length(class) > 2){
      	for (k in 1:(length(class)-1)){
        	for (l in (k+1):length(class)){
            	print(paste(class[k],"vs",class[l],sep=" ")) 
            	plot(0, 0, type = "n", xlim = xlim/60, ylim = ylim, main = paste("Base Peak Chromatograms \n","BPCs_",class[k]," vs ",class[l], sep=""), xlab = "Retention Time (min)", ylab = "BPC")
            	colvect <- NULL
           		for (j in 1:length(classnames[[k]])) {
      
              		tic <- TIC[[classnames[[k]][j]]]
              		# points(tic[,1]/60, tic[,2], col = cols[i], pch = pch[i], type="l")
              		points(tic[,1]/60, tic[,2], col = cols[classnames[[k]][j]], pch = pch[classnames[[k]][j]], type = "l")
              		colvect <- append(colvect,cols[classnames[[k]][j]])
            	}
          		for (j in 1:length(classnames[[l]])) {
          			# i=class2names[j]
          			tic <- TIC[[classnames[[l]][j]]]
          			points(tic[,1]/60, -tic[,2], col = cols[classnames[[l]][j]], pch = pch[classnames[[l]][j]], type = "l")
          			colvect <- append(colvect,cols[classnames[[l]][j]])
          		}
          		legend("topright",paste(names(files[c(classnames[[k]],classnames[[l]])])), col = colvect, lty = lty, pch = pch)
        	}
      	}
    }#end if length >2
    
    if (length(class) == 2){
        k = 1
		l = 2
        colvect <- NULL
        plot(0, 0, type = "n", xlim = xlim/60, ylim = ylim, main = paste("Base Peak Chromatograms \n","BPCs_",class[k],"vs",class[l], sep=""), xlab = "Retention Time (min)", ylab = "BPC")

       	for (j in 1:length(classnames[[k]])) {

          	tic <- TIC[[classnames[[k]][j]]]
          	# points(tic[,1]/60, tic[,2], col = cols[i], pch = pch[i], type="l")
          	points(tic[,1]/60, tic[,2], col = cols[classnames[[k]][j]], pch = pch[classnames[[k]][j]], type = "l")
          	colvect <- append(colvect,cols[classnames[[k]][j]])
        }
      	for (j in 1:length(classnames[[l]])) {
          	# i=class2names[j]
          	tic <- TIC[[classnames[[l]][j]]]
          	points(tic[,1]/60, -tic[,2], col = cols[classnames[[l]][j]], pch = pch[classnames[[l]][j]], type = "l")
          	colvect <- append(colvect,cols[classnames[[l]][j]])
      	}
      	legend("topright",paste(names(files[c(classnames[[k]],classnames[[l]])])), col = colvect, lty = lty, pch = pch)

    }#end length ==2
    
    if (length(class) == 1){
        k = 1
		ylim = range(sapply(TIC, function(x) range(x[,2])))
        colvect <- NULL
        plot(0, 0, type = "n", xlim = xlim/60, ylim = ylim, main = paste("Base Peak Chromatograms \n","BPCs_",class[k], sep=""), xlab = "Retention Time (min)", ylab = "BPC")

       	for (j in 1:length(classnames[[k]])) {

          	tic <- TIC[[classnames[[k]][j]]]
          	# points(tic[,1]/60, tic[,2], col = cols[i], pch = pch[i], type="l")
          	points(tic[,1]/60, tic[,2], col = cols[classnames[[k]][j]], pch = pch[classnames[[k]][j]], type = "l")
          	colvect <- append(colvect,cols[classnames[[k]][j]])
        }
      
        legend("topright",paste(names(files[c(classnames[[k]])])), col = colvect, lty = lty, pch = pch)

    }#end length ==1
    dev.off()

    # invisible(TIC)
}

getTIC <- function(file,rtcor=NULL) {
    object <- xcmsRaw(file)
    cbind(if (is.null(rtcor)) object@scantime else rtcor, rawEIC(object,mzrange=range(object@env$mz))$intensity)
}

##
##  overlay TIC from all files in current folder or from xcmsSet, create pdf
##
getTIC2s <- function(files, pdfname="TICs.pdf", rt=c("raw","corrected")) {
         
    #create sampleMetadata, get sampleMetadata and class
    sampleMetadata<-xcms:::phenoDataFromPaths(files)
    class<-class<-as.vector(levels(sampleMetadata[,"class"])) #create phenoData like table
    classnames<-vector("list",length(class))
    for (i in 1:length(class)){
        classnames[[i]]<-which( sampleMetadata[,1]==class[i])
    }
        
    N <- dim(sampleMetadata)[1]
    TIC <- vector("list",N)

    for (i in 1:N) {
        cat(files[i],"\n")
        if (!is.null(xcmsSet) && rt == "corrected"){
        	rtcor <- xcmsSet@rt$corrected[[i]]
        }else{
        	rtcor <- NULL
        }
        TIC[[i]] <- getTIC(files[i],rtcor=rtcor)
    }
 
    pdf(pdfname,w=16,h=10)
    cols <- rainbow(N)
    lty = 1:N
    pch = 1:N
    #search for max x and max y in TICs
    xlim = range(sapply(TIC, function(x) range(x[,1])))
    ylim = range(sapply(TIC, function(x) range(x[,2])))
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
                legend("topright",paste(names(files[c(classnames[[k]],classnames[[l]])])), col = colvect, lty = lty, pch = pch)
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
        legend("topright",paste(names(files[c(classnames[[k]],classnames[[l]])])), col = colvect, lty = lty, pch = pch)

    }#end length ==2

    if (length(class)==1){
        k=1
        ylim = range(sapply(TIC, function(x) range(x[,2])))
		
        plot(0, 0, type="n", xlim = xlim/60, ylim = ylim, main = paste("Total Ion Chromatograms \n","TICs_",class[k], sep=""), xlab = "Retention Time (min)", ylab = "TIC")
        colvect<-NULL
        for (j in 1:length(classnames[[k]])) {

            tic <- TIC[[classnames[[k]][j]]]
            # points(tic[,1]/60, tic[,2], col = cols[i], pch = pch[i], type="l")
            points(tic[,1]/60, tic[,2], col = cols[classnames[[k]][j]], pch = pch[classnames[[k]][j]], type="l")
            colvect<-append(colvect,cols[classnames[[k]][j]])
        }

        legend("topright",paste(names(files[c(classnames[[k]])])), col = colvect, lty = lty, pch = pch)

    }#end length ==1

    dev.off()

    # invisible(TIC)
}


##addition for quality control of peak picking
#metaMS EIC and pspectra plotting option
#version 20150512
#only for Galaxy 

plotUnknowns<-function(resGC, unkn=""){

    ##Annotation table each value is a pcgrp associated to the unknown 
    ##NOTE pcgrp index are different between xcmsSet and resGC due to filtering steps in metaMS
    ##R. Wehrens give me some clues on that and we found a correction
    #if unkn="none"
	
	if(unkn=="none") {
	   	pdf("Unknown_Empty.pdf")
	   	plot.new()
	   	text(x=0.5,y=1,pos=1, labels="No EIC ploting required")
	   	dev.off()
	}else{

		mat<-matrix(ncol=length(resGC$xset), nrow=dim(resGC$PeakTable)[1])
		 
		for (j in 1: length(resGC$xset)){
			test<-resGC$annotation[[j]]
			print(paste("j=",j))
			for (i in 1:dim(test)[1]){
				if (as.numeric(row.names(test)[i])>dim(mat)[1]){
					next
				} else {
					mat[as.numeric(row.names(test)[i]),j]<-test[i,1]
				}
			}
		}
		colnames(mat)<-colnames(resGC$PeakTable[,c((which(colnames(resGC$PeakTable)=="rt"|colnames(resGC$PeakTable)=="RI")[length(which(colnames(resGC$PeakTable)=="rt"|colnames(resGC$PeakTable)=="RI"))]+1):dim(resGC$PeakTable)[2])])
		
		#debug

		# print(dim(mat))
		# print(mat[1:3,]) 
		# write.table(mat, file="myannotationtable.tsv", sep="\t", row.names=FALSE)
		#correction of annotation matrix due to pcgrp removal by quality check in runGCresult
		#matrix of correspondance between an@pspectra and filtered pspectra from runGC

		allPCGRPs <- lapply(1:length(resGC$xset),
						function(i) {
							an <- resGC$xset[[i]]
							huhn <- an@pspectra[which(sapply(an@pspectra, length) >= metaSetting(resGC$settings,"DBconstruction.minfeat"))]
							matCORR<-cbind(1:length(huhn), match(huhn, an@pspectra))
						})

		if (unkn[1]==""){    
		#plot EIC and spectra for all unknown for comparative purpose
	   

			par (mar=c(5, 4, 4, 2) + 0.1)
			for (l in 1:dim(resGC$PeakTable)[1]){ #l=2
				#recordPlot
				perpage=3 #if change change layout also!
				num.plots <- ceiling(dim(mat)[2]/perpage) #three pcgroup per page
				my.plots <- vector(num.plots, mode='list')
				dev.new(width=21/2.54, height=29.7/2.54, file=paste("Unknown_",l,".pdf", sep="")) #A4 pdf
				# par(mfrow=c(perpage,2))
				layout(matrix(c(1,1,2,3,4,4,5,6,7,7,8,9), 6, 2, byrow = TRUE), widths=rep(c(1,1),perpage), heights=rep(c(1,5),perpage))
				# layout.show(6)
				oma.saved <- par("oma")
				par(oma = rep.int(0, 4))
				par(oma = oma.saved)
				o.par <- par(mar = rep.int(0, 4))
				on.exit(par(o.par))
				stop=0 #initialize
				for (i in 1:num.plots) {
					start=stop+1
					stop=start+perpage-1 #
					for (c in start:stop){
						if (c <=dim(mat)[2]){
								
							#get sample name
							sampname<-basename(resGC$xset[[c]]@xcmsSet@filepaths)

							#remove .cdf, .mzXML filepattern
							filepattern <- c("[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]", 
									"[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]")
							filepattern <- paste(paste("\\.", filepattern, "$", sep = ""), 
									collapse = "|")
							sampname<-gsub(filepattern, "",sampname)
							 
							title1<-paste("unknown", l,"from",sampname, sep=" ")
							an<-resGC$xset[[c]]
							 
							par (mar=c(0, 0, 0, 0) + 0.1)
							plot.new()
							box()
							text(0.5, 0.5, title1, cex=2)
							if (!is.na(mat[l,c])){
								pcgrp=allPCGRPs[[c]][which(allPCGRPs[[c]][,1]==mat[l,c]),2]
								if (pcgrp!=mat[l,c]) print ("pcgrp changed")
								par (mar=c(3, 2.5, 3, 1.5) + 0.1)
								plotEICs(an, pspec=pcgrp, maxlabel=2)
								plotPsSpectrum(an, pspec=pcgrp, maxlabel=2)
							} else {
								plot.new()
								box()
								text(0.5, 0.5, "NOT FOUND", cex=2)
								plot.new()
								box()
								text(0.5, 0.5, "NOT FOUND", cex=2)
							}
						}
					}
					# my.plots[[i]] <- recordPlot()
				}
				graphics.off()

					# pdf(file=paste("Unknown_",l,".pdf", sep=""), onefile=TRUE)
					# for (my.plot in my.plots) {
						# replayPlot(my.plot)
					# }
					# my.plots
					# graphics.off()

			}#end  for l
		}#end if unkn=""
		else{
			par (mar=c(5, 4, 4, 2) + 0.1)
			l=unkn
			if (length(l)==1){
				#recordPlot
				perpage=3 #if change change layout also!
				num.plots <- ceiling(dim(mat)[2]/perpage) #three pcgroup per page
				my.plots <- vector(num.plots, mode='list')
				
				dev.new(width=21/2.54, height=29.7/2.54, file=paste("Unknown_",l,".pdf", sep="")) #A4 pdf
				# par(mfrow=c(perpage,2))
				layout(matrix(c(1,1,2,3,4,4,5,6,7,7,8,9), 6, 2, byrow = TRUE), widths=rep(c(1,1),perpage), heights=rep(c(1,5),perpage))
				# layout.show(6)
				oma.saved <- par("oma")
				par(oma = rep.int(0, 4))
				par(oma = oma.saved)
				o.par <- par(mar = rep.int(0, 4))
				on.exit(par(o.par))
				stop=0 #initialize
				for (i in 1:num.plots) {
					start=stop+1
					stop=start+perpage-1 #
					for (c in start:stop){
						if (c <=dim(mat)[2]){
								
							#get sample name
							sampname<-basename(resGC$xset[[c]]@xcmsSet@filepaths)

							#remove .cdf, .mzXML filepattern
							filepattern <- c("[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]", "[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]")
							filepattern <- paste(paste("\\.", filepattern, "$", sep = ""), collapse = "|")
							sampname<-gsub(filepattern, "",sampname)
							 
							title1<-paste("unknown", l,"from",sampname, sep=" ")
							an<-resGC$xset[[c]]
							 
							par (mar=c(0, 0, 0, 0) + 0.1)
							plot.new()
							box()
							text(0.5, 0.5, title1, cex=2)
							if (!is.na(mat[l,c])){
								pcgrp=allPCGRPs[[c]][which(allPCGRPs[[c]][,1]==mat[l,c]),2]
								if (pcgrp!=mat[l,c]) print ("pcgrp changed")
								par (mar=c(3, 2.5, 3, 1.5) + 0.1)
								plotEICs(an, pspec=pcgrp, maxlabel=2)
								plotPsSpectrum(an, pspec=pcgrp, maxlabel=2)
							} else {
								plot.new()
								box()
								text(0.5, 0.5, "NOT FOUND", cex=2)
								plot.new()
								box()
								text(0.5, 0.5, "NOT FOUND", cex=2)
							}
						}
					}
					# my.plots[[i]] <- recordPlot()
				}
				graphics.off()

				# pdf(file=paste("Unknown_",l,".pdf", sep=""), onefile=TRUE)
				# for (my.plot in my.plots) {
					# replayPlot(my.plot)
				# }
				# my.plots
				# graphics.off()
			} else {
				par (mar=c(5, 4, 4, 2) + 0.1)
				for (l in 1:length(unkn)){ #l=2
					#recordPlot
					perpage=3 #if change change layout also!
					num.plots <- ceiling(dim(mat)[2]/perpage) #three pcgroup per page
					my.plots <- vector(num.plots, mode='list')
					dev.new(width=21/2.54, height=29.7/2.54, file=paste("Unknown_",unkn[l],".pdf", sep="")) #A4 pdf
					# par(mfrow=c(perpage,2))
					layout(matrix(c(1,1,2,3,4,4,5,6,7,7,8,9), 6, 2, byrow = TRUE), widths=rep(c(1,1),perpage), heights=rep(c(1,5),perpage))
					# layout.show(6)
					oma.saved <- par("oma")
					par(oma = rep.int(0, 4))
					par(oma = oma.saved)
					o.par <- par(mar = rep.int(0, 4))
					on.exit(par(o.par))
					stop=0 #initialize
					for (i in 1:num.plots) {
						start=stop+1
						stop=start+perpage-1 #
						for (c in start:stop){
							if (c <=dim(mat)[2]){
							
								#get sample name
								sampname<-basename(resGC$xset[[c]]@xcmsSet@filepaths)

								#remove .cdf, .mzXML filepattern
								filepattern <- c("[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]", "[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]")
								filepattern <- paste(paste("\\.", filepattern, "$", sep = ""), collapse = "|")
								sampname<-gsub(filepattern, "",sampname)
								 
								title1<-paste("unknown",unkn[l],"from",sampname, sep=" ")
								an<-resGC$xset[[c]]
								 
								par (mar=c(0, 0, 0, 0) + 0.1)
								plot.new()
								box()
								text(0.5, 0.5, title1, cex=2)
								if (!is.na(mat[unkn[l],c])){
									pcgrp=allPCGRPs[[c]][which(allPCGRPs[[c]][,1]==mat[unkn[l],c]),2]
									if (pcgrp!=mat[unkn[l],c]) print ("pcgrp changed")
									par (mar=c(3, 2.5, 3, 1.5) + 0.1)
									plotEICs(an, pspec=pcgrp, maxlabel=2)
									plotPsSpectrum(an, pspec=pcgrp, maxlabel=2)
								} else {
									plot.new()
									box()
									text(0.5, 0.5, "NOT FOUND", cex=2)
									plot.new()
									box()
									text(0.5, 0.5, "NOT FOUND", cex=2)
								}
							}
						}
						# my.plots[[i]] <- recordPlot()
					}
					graphics.off()

					# pdf(file=paste("Unknown_",unkn[l],".pdf", sep=""), onefile=TRUE)
					# for (my.plot in my.plots) {
						# replayPlot(my.plot)
					# }
					# my.plots
					# graphics.off()

				}#end  for unkn[l]
			
			}
		
		}
	}
} #end function 

# This function get the raw file path from the arguments
#@author Gildas Le Corguille lecorguille@sb-roscoff.fr
getRawfilePathFromArguments <- function(singlefile, zipfile, args, prefix="") {

  if (!(prefix %in% c("","Positive","Negative","MS1","MS2", "XSET"))) stop("prefix must be either '', 'Positive', 'Negative', 'MS1' ,'MS2' or 'XSET'")

  if (!is.null(args[[paste0("zipfile",prefix)]])) zipfile <- args[[paste0("zipfile",prefix)]]

  if (!is.null(args[[paste0("singlefile_galaxyPath",prefix)]])) {
    singlefile_galaxyPaths <- args[[paste0("singlefile_galaxyPath",prefix)]]
    singlefile_sampleNames <- args[[paste0("singlefile_sampleName",prefix)]]
  }

  if (exists("singlefile_galaxyPaths")){
    singlefile_galaxyPaths <- unlist(strsplit(singlefile_galaxyPaths,"\\|"))
    singlefile_sampleNames <- unlist(strsplit(singlefile_sampleNames,"\\|"))

    singlefile <- NULL
    for (singlefile_galaxyPath_i in seq(1:length(singlefile_galaxyPaths))) {
      singlefile_galaxyPath <- singlefile_galaxyPaths[singlefile_galaxyPath_i]
      singlefile_sampleName <- singlefile_sampleNames[singlefile_galaxyPath_i]
      # In case, an url is used to import data within Galaxy
      singlefile_sampleName <- tail(unlist(strsplit(singlefile_sampleName,"/")), n=1)
      singlefile[[singlefile_sampleName]] <- singlefile_galaxyPath
    }
  }
  return(list(zipfile=zipfile, singlefile=singlefile))
}


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