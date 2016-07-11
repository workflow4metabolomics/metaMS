####----------------------------------------------------------------------####
####                  GC-MS data processing version 0.99.5                ####
####----------------------------------------------------------------------####

I. Date: 2015-05-26

II. Author and maintainer:

   Yann GUITTON
   LINA - équipe Combi (CNRS, Univ-Nantes, EMN, INRIA)
   IRISA -  équipe Dyliss (CNRS, Univ-Rennes 1, INRIA) 
   Phone: +33 2 51 12 53 90
   E-mail: yann.guitton@irisa.fr / yann.guitton@gmail.com

III. Funding

   Developed within IDEALG project (http://www.idealg.ueb.eu/versionAnglaise/)

IV. Usage restrictions

   Use of this tool is restricted to the service conditions of the MetaboHUB-IFB infrastructures.
   For any question regarding the use of these services, please contact: yann.guitton@univ-nantes.fr

V. Installation

   4 files are required for installation:

   1) 'README.txt'
         Instructions for installation
   
   2) 'idealg_metams_runGC.xml'
         Configuration file; to be put into the './galaxy-dist/tools/' directory
		+ 2.png files for illustration

   3) 'metams.R'
         Wrapper code written in R aimed at launching the runGC function from the metaMS package given the arguments entered by the user through the Galaxy interface
   
   4) 'metaMS R package '
         The 'metaMS' package requires dependencies and can be installed with 
                source("http://bioconductor.org/biocLite.R")
                biocLite("metaMS")
 
         This code is for installation of the Galaxy module on the Workflow4metabolomics.org MetaboHUB-IFB platform only and must not be distributed without the author agreement

   
        
