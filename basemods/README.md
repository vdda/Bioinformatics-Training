Bacterial_Basemod_Analysis
==========================

This repository contains R functions and several walk-through examples of how to use those functions to do custom analysis of the base modification data output by SMRTPortal.  The R functions are focused on studying methylation patterns created by the methyltransferase / restriction enzyme (MT-RE) systems found in bacterial genomes.  Additionally, the exercise in the kinetics folder reviews the use of functions available in the R library pbh5 (also available from github) that would likely be useful in analyzing other types of methylation.

To complete the exercises, you will need access to R, free software available here: http://www.r-project.org/.  New users may find RGui easier to work with than command line versions.  RGui can be used on Windows, Mac, and Unix platforms, though to use it on a unix platform your sysadmin will have to install it for you. You will also need to install the following libraries in R: Hmisc, ggplot2, Biostrings, plyr, seqinr and pbh5.  All of these are available from the R cran mirrors, with the exception of pbh5, a PacBio library that is available on github. The exercises were developed with R 2.15.2, and library versions ggplot2_0.9.2.1, plyr_1.7.1, Hmisc_3.9-3, and Biostrings_2.24.1. Use of other versions may result in apparent bugs.

Having secured access to R and installed the required libraries, download this repository. For the kinetics exercise, you will also need to download the aligned_reads.cmp.h5 example data available here:  [basemod data](https://datasets.pacb.com.s3.amazonaws.com/2013/basemod_workshop/basemod_dataset.tgz).  Put the aligned_reads.cmp.h5 file into the kinetics folder in your local copy of this repository.

