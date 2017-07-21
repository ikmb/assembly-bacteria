![](images/ikmb_bfx_logo.png)
# IKMB Genomes - Bacterial assembly and annotation from PE reads

## Overview

This pipeline takes PE genomic reads, performs trimming using trim_galore,  assembles draft genomes using Spades and annotates the assembly using Prokka. 

In addition, reads will be mapped back to the genome sequence to obtain stats about mapping rate, insert size, duplication rates, etc. 

## Running the pipeline

The pipeline is currently configured to run on RZCluster. Other execute environment may be added in the future. 

To run the pipeline, you first have to clone the repository from git: 

`git clone git@git.ikmb.uni-kiel.de:bfx-core/NF-genomes-bacteria.git`

If you already have a local copy, make sure it is up-to-date:

`git pull`

With a local installation of the pipeline, you can next proceed to run it in a folder where you would like to generate the various result files:

Load the environment modules for Nextflow:

`module load IKMB Java Nextflow`

Next, move to the folder where the pipeline will be run from:

`cd /path/to/run_folder`

Finally, run the pipeline by providing a path to the sample info:

`nextflow -c /path/to/nextflow.config run /path/to/main.nf --samples Samples.csv`

The file Samples is a CSV formatted list of data sets to process. An example format is included in the template subfolder. Required information include:

* Sample name
* Library name
* Path to forward read
* path to reverse read

## Outputs

The pipeline generates all outputs in the folder "outputs". For each PE data set (a library), a folder will be created containing the assembly, the annotation and the re-mapped reads for computing statistics on coverage etc. 

Furthermore, trimmed reads will be jointly placed in the folder "trimgalore". A compilation of statistics for the various analysis is included in the older MultiQC - for library-level analytics (library_multiqc.html) and the assembly, annotation etc (sample_multiqc.html).




