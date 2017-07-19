# IKMB Genomes - Bacterial assembly and annotation from PE reads

## Overview

This pipeline takes PE genomic reads, performs trimming using trim_galore,  assembles draft genomes using Spades and annotates the assembly using Prokka. 

In addition, reads will be mapped back to the genome sequence to obtain stats about mapping rate, insert size, duplication rates, etc. 

## Running the pipeline

The pipeline is currently configured to run on RZCluster. Other execute environment may be added in the future. 

To run the pipeline, you first have to clone the repository from git: 

`git clone git@git.ikmb.uni-kiel.de:genomes/genomes-bacteria.git`

If you already have a local copy, make sure it is up-to-date:

`git pull`

With a local installation of the pipeline, you can next proceed to run it in a folder where you would like to generate the various result files:

Load the environment modules for Nextflow:

`module load IKMB Java Nextflow`

Next, move to the folder where the pipeline will be run from:

`cd /path/to/run_folder`

Finally, run the pipeline:

`nextflow -c /path/to/nextflow.config run /path/to/main.nf --folder /path/to/reads`

Where "path/to/reads" points to the folder in which the PE files are located. 

## Outputs

The pipeline generates all outputs in the folder "outputs". For each PE data set (a library), a folder will be created containing the assembly, the annotation and the re-mapped reads for computing statistics on coverage etc. 

Furthermore, trimmed reads will be jointly placed in the folder "trimgalore". A compilation of statistics for the various analysis is included in the older MultiQC - for library-level analytics (library_multiqc.html) and the assembly, annotation etc (sample_multiqc.html).




