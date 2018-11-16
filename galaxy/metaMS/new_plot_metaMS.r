#!/usr/local/public/bin/Rscript --vanilla --slave --no-site-file
# metams.r version="2.1.2"
#created by Yann GUITTON 
#use RI options + add try on plotUnknown add session Info
#use make.names in sampleMetadata to avoid issues with files names 


# ----- LOG FILE -----
log_file <- file("plot_metaMS.log", open = "wt")
sink(log_file)
sink(log_file, type = "output")


# ----- PACKAGE -----
cat("\tSESSION INFO\n")

#Import the different functions
source_local <- function(fname) {
    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    source(paste(base_dir, fname, sep="/"))
}
source_local("lib_metams.r")

pkgs <- c("metaMS","batch") #"batch" necessary for parseCommandArgs function
loadAndDisplayPackages(pkgs)

cat("\n")

modNamC <- "plot_metaMS" ## module name

cat("\nStart of the '", modNamC, "' Galaxy module call: ", 
    format(Sys.time(), "%a %d %b %Y %X"), "\n\n", sep="")


# ----- ARGUMENTS -----
cat("\tARGUMENTS INFO\n")
args = parseCommandArgs(evaluate=FALSE) #interpretation of arguments given in command line as an R list of objects
write.table(as.matrix(args), col.names=F, quote=F, sep='\t\t')


# ----- PROCESSING INFILE -----
cat("\n\n\tARGUMENTS PROCESSING INFO\n")

cat("\n\n")


# ----- INFILE PROCESSING -----
cat("\tINFILE PROCESSING INFO\n\n")

# Loading RData file
load(args[["metaMS"]])
if (!exists("resGC")) stop("\n\nERROR: The RData doesn't contain any object called 'resGC' which is provided by the tool: new_metaMS.runGC")

# Handle infiles
if (!exists("singlefile")) singlefile <- NULL
if (!exists("zipfile")) zipfile <- NULL
rawFilePath <- getRawfilePathFromArguments(singlefile, zipfile, args)
zipfile <- rawFilePath$zipfile
singlefile <- rawFilePath$singlefile
directory <- retrieveRawfileInTheWorkingDirectory(singlefile, zipfile)


# ----- MAIN PROCESSING INFO -----
cat("\n\n\tMAIN PROCESSING INFO\n")


cat("\t\tCOMPUTE\n")

#TIC/BPC picture generation

#Use getTIC2s and getBPC2s because getTICs and getBPCs can exists due to transfert of function in Rdata
cat("\nProcessing BPCs...\n")
if(!is.null(singlefile)){
    print("singlefile BPC")
    files <- paste("./",names(singlefile),sep="")
    if(!is.null(files)){
        c <- getBPC2s(files = files, xset = xset, rt="raw", pdfname="BPCs_raw.pdf")
    }else{
        print("Error files is empty")
    }
}
if(!is.null(zipfile)){
    print("zipfile BPC")
    files <- getMSFiles(directory)
    if(!is.null(files)){
        c <- getBPC2s(files = files, rt="raw", pdfname="BPCs_raw.pdf") 
    }else{
        print("Error files is empty")
    } 
}
cat("BPC created...\n")

cat("\nProcessing TICs...\n")
if(!is.null(singlefile)){
    print("singlefile TIC")
    files <- paste("./",names(singlefile),sep="")
    if(!is.null(files)){
        b <- getTIC2s(files = files, xset = xset, rt="raw", pdfname="TICs_raw.pdf")
    }else{
        print("Error files is empty")
    }
}
if(!is.null(zipfile)){
    print("zipfile TIC")
    files <- getMSFiles(directory)
    if(!is.null(files)){
        b <- getTIC2s(files = files, rt="raw", pdfname="TICs_raw.pdf")  
    }else{
        print("Error files is empty")
    }
}
cat("TIC created...\n")

#test
system(paste('gs  -o TICsBPCs_merged.pdf  -sDEVICE=pdfwrite  -dPDFSETTINGS=/prepress  *Cs_raw.pdf'))
system(paste('gs  -o GCMS_EIC.pdf  -sDEVICE=pdfwrite  -dPDFSETTINGS=/prepress  Unknown_*'))

cat("\nEnd of '", modNamC, "' Galaxy module call: ", as.character(Sys.time()), "\n", sep = "")