#!/usr/local/public/bin/Rscript --vanilla --slave --no-site-file
# metams.r version="2.1.2"
#created by Yann GUITTON 
#use RI options + add try on plotUnknown add session Info
#use make.names in sampleMetadata to avoid issues with files names 

# ----- LOG FILE -----
#log_file <- file("plot_metams.log", open = "wt")
#sink(log_file)
#sink(log_file, type = "output")

# ----- PACKAGE -----
cat("\tSESSION INFO\n")

#Import the different functions
source_local <- function(fname) {
    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    print(paste(base_dir, fname, sep="/"))
    source(paste(base_dir, fname, sep="/"))
}
source_local("lib_metams.r")

pkgs <- c("metaMS","batch") #"batch" necessary for parseCommandArgs function
loadAndDisplayPackages(pkgs)

cat("\n\n")

modNamC <- "plot_metaMS" ## module name

## log file
##---------
cat("\nStart of the '", modNamC, "' Galaxy module call: ", 
    format(Sys.time(), "%a %d %b %Y %X"), "\n\n", sep="")

# ----- ARGUMENTS -----
cat("\tARGUMENTS INFO\n")
listArguments = parseCommandArgs(evaluate=FALSE) 
#interpretation of arguments given in command line as an R list of objects
cat("\nParameters:\n")
print(listArguments)
cat("\n\n")

# ----- PROCESSING INFILE -----
cat("\tARGUMENTS PROCESSING INFO\n")

# Loading RData file
load(listArguments[["metaMS"]])

# Saving the specific parameters

#TIC/BPC picture generation
 # use files as entry not xset that do not exist   

#Use getTIC2s and getBPC2s because getTICs and getBPCs can exists due to transfert of function in Rdata
print("Processing BPCs...")
print(filepathsXset)
c <- getBPC2s(files = filepathsXset, rt="raw", pdfname="BPCs_raw.pdf")  

print("Processing TICs...")
b <- getTIC2s(files = filepathsXset, rt="raw", pdfname="TICs_raw.pdf")

print("Step QC plot")

#to do check if no peaks found
#Quality controls plots but only working in R (don't know why)
a <- try(plotUnknowns(resGC=resGC, unkn=unknarg)); #use unknparam value
if(class(a) == "try-error") {
    pdf("Unknown_Error.pdf")
    plot.new()
    text(x=0.5,y=1,pos=1, labels="Error generating EICs\n please use none instead of a vector in plotUnknown")
    dev.off()
}
# create a mergpdf

#test
system(paste('gs  -o TICsBPCs_merged.pdf  -sDEVICE=pdfwrite  -dPDFSETTINGS=/prepress  *Cs_raw.pdf'))
system(paste('gs  -o GCMS_EIC.pdf  -sDEVICE=pdfwrite  -dPDFSETTINGS=/prepress  Unknown_*'))

#delete the parameters to avoid the passage to the next tool in .RData image
rm(listArguments)

## Closing
##--------

cat("\nEnd of '", modNamC, "' Galaxy module call: ", as.character(Sys.time()), "\n", sep = "")

rm(list = ls())
#sink()