Bacterial_Basemod_Analysis
==========================

This repository contains R functions and a walk-through example of using those functions to do additional custom analysis of the base modification data output by SMRTPortal.  The functions are focused on studying methylation patterns within methyltransferase recognition motifs commonly found in bacterial genomes.  However, the pbh5tools.R exercise reviews the use of functions that would likely be useful in analyzing other types of methylation.

To complete the exercises, you will need access to R, free software available here: http://www.r-project.org/.  New users may find RGui to be easier to work with than command line versions.  RGui can be used on Windows, Mac, and Unix platforms, though for unix platforms your sysadmin will have to install it for you. 

Having secured access to R, download the files in this directory to one folder and create three additional empty folders in that same folder called 'analysis', 'circos', and 'cmph5_demo'.  Do not rename any of the files you have downloaded.  The R files you downloaded contain links to the filenames, and will use the empty folders you just made to write the pdfs and csvs you will generate as you work through the exercises.  

To being, open basemodAnalysis.R and work your way through the exercises therein by executing sections of code and reading the commentary.  Note that to make circos plots, you will need to exit R and run circos at the command line of a unix system where circos is installed, within the directory containing all circos text files produced during the exercise ( /circos).  Once complete, do the same with pbh5tools.R
