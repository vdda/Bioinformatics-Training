HGAP 1.4 wrapper script
-----------------------
This is a convenience wrapper script to execute the HGAP assembly process using  
the 1.4 SMRT Pipe software release.  Its designed to execute the entire pipeline  
in one go with minimal input required from the user.

Files
-----
* hgap-1.4.sh - The main wrapper script
* params_preasm.xml - Default parameters for PreAssembler
* params_reseq.xml - Default parameters for Resequencing

Download these to a directory and include it in your PATH environment variable.   


Command Line
------------
Be sure to have your SMRT Pipe installation and environment ready prior to  
execution.  You first need to set your SEYMOUR_HOME env variable to the root   
installation of SMRT Analysis, e.g., /opt/smrtanalysis.  Then you must source the  
environment setup script.  This is necessary because the hgap14.sh script   
make use of programs and tools installed with SMRT Analysis.

    > export SEYMOUR_HOME=/opt/smrtanalysis
    > source $SEYMOUR_HOME/etc/setup.sh

You can then proceed to execute the hgap14.sh script:

    > hgap14.sh --help
    USAGE: hgap14.sh [params] <input.xml>
     -p    Optional path to a preassembler params file
     -r    Optional path to a resequencing params file
     -s    Optional path to a celera-assembler spec file
     -x    Override default options to smrtpipe
     -l    Run everything locally (e.g., no cluster)
     -d    Keep all the intermediate files and output
           some statistics at the end


You must provide an input.xml which lists the base files to include in the    
analysis.

If the -p and/or -r options are provided, the default parameter files are   
overridden.  The default files are 'params_preasm.xml' and 'params_reseq.xml'.
These must be present alongside the hgap14.sh for the defaults to work.

A recommended default spec file 'ca_default.spec' for runCA.  If present in   
the same directory as hgap14.sh, it will be used.

If you do not have access to a cluster or the host is not a submit host, then  
you may supply the -l flag to run everything locally.

The -d flag will keep some extra files around, but also generate some useful   
statistics at the end to help you evaluate the finished assembly.

Some success has been achieved by using a newer development versions of runCA   
than the one we currently ship with (7.0).  The runCA executable called by the   
script can be overriden on the command line like this:

    > RunCA=/path/to/custom/runCA hgap14.sh 


Output
------
The final output of HGAP will be located in the data/ directory:
* consensus.fasta.gz
* consensus.fastq.gz

Notes
-----

The provided xml parameter files are suitable for C2 data.  For XL data,
you may find more success by adding --allowPartialAlignments to the 
layoutOps section of the params_preasm.xml file.

Setting the overlap mer size smaller has also helped in some cases with XL  
data.
