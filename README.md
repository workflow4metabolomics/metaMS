metaMS for Galaxy
=================

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/bioconductor-metams/README.html) [![Build Status](https://travis-ci.org/workflow4metabolomics/metaMS.svg?branch=master)](https://travis-ci.org/workflow4metabolomics/metaMS)

Our project
-----------
The [Workflow4Metabolomics](http://workflow4metabolomics.org), W4M in short, is a French infrastructure offering software tool processing, analyzing and annotating metabolomics data. It is based on the Galaxy platform.


metaMS
------
MS-based metabolomics data processing and compound annotation pipeline.

Author: Ron Wehrens [aut, cre] (author of GC-MS part), Pietro Franceschi [aut] (author of LC-MS part), Nir Shahaf [ctb], Matthias Scholz [ctb], Georg Weingart [ctb] (development of GC-MS approach), Elisabete Carvalho [ctb] (testing and feedback of GC-MS pipeline)

Maintainer: Ron Wehrens <ron.wehrens at gmail.com>

Citation (from within R, enter citation("metaMS")):

Wehrens R, Weingart G and Mattivi F (2014). “metaMS: An open-source pipeline for GC-MS-based untargeted metabolomics.” J. Chrom. B, 966, pp. 109-116.

Homepage: [https://github.com/Bioconductor-mirror/metaMS/](https://github.com/Bioconductor-mirror/metaMS/tree/release-3.3)


Galaxy
------
Galaxy is an open, web-based platform for data intensive biomedical research. Whether on the free public server or your own instance, you can perform, reproduce, and share complete analyses. 

Homepage: [https://galaxyproject.org/](https://galaxyproject.org/)

Dependencies using Conda
------------------------
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/bioconductor-metams/README.html)

[Conda](http://conda.pydata.org/) is package manager that among many other things can be used to manage Python packages.


```
#To install miniconda2
#http://conda.pydata.org/miniconda.html
#To install the metaMS R library using conda:
conda install bioconductor-metams bioconductor-xcms r-batch
#To set an environment:
conda create -n bioconductor-metams bioconductor-metams bioconductor-xcms r-batch
#To activate the environment:
. activate bioconductor-metams
```

Travis
------
[![Build Status](https://travis-ci.org/workflow4metabolomics/metaMS.svg?branch=master)](https://travis-ci.org/workflow4metabolomics/metaMS)

Test and Deploy with Confidence. Easily sync your GitHub projects with Travis CI and you'll be testing your code in minutes!


Test Status
-----------

Planemo test using conda: passed


Historic contributors
---------------------
 - Yann Guitton @yguitton - [LABERCA - Laboratory of Food Contaminants and Residue Analysis](http://www.laberca.org/) - Ecole Nationale Vétérinaire, Agroalimentaire et de l'Alimentation Nantes-Atlantique - France
 - Gildas Le Corguillé @lecorguille - [ABiMS](http://abims.sb-roscoff.fr/) / [IFB](http://www.france-bioinformatique.fr/) - [UPMC](www.upmc.fr)/[CNRS](www.cnrs.fr) - [Station Biologique de Roscoff](http://www.sb-roscoff.fr/) - France

